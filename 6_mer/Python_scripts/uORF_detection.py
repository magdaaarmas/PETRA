from itertools import product
import pandas as pd
import os
import sys

# ----------------------------- args -----------------------------
gene = sys.argv[1]
UTR_info_file = sys.argv[2]   # '5UTR_information.csv'
save_folder = sys.argv[3]     # ensure it ends with '/' in the caller
# ----------------------------------------------------------------

def generate_dna_sequences(length):
    bases = ['A', 'C', 'G', 'T']
    return [''.join(seq) for seq in product(bases, repeat=length)]

def detect_stop_codons(seq_from_start):
    """
    seq_from_start: DNA sequence starting at the ATG (3' direction).
    Returns index (0-based) of first stop codon end *relative to seq_from_start start*,
    or False if no in-frame stop found.
    """
    seq = seq_from_start.upper()
    stop_codons = {'TGA', 'TAG', 'TAA'}
    # step in codons of 3
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3]
        if len(codon) < 3:
            break
        if codon in stop_codons:
            return i  # position of first base of stop codon
    return False

def get_insertion_context(df):
    """Assumes df has columns '5UTR_sequence' and 'insert_position' (1-based)."""
    if df.empty:
        raise ValueError("No rows for the requested gene in UTR info.")
    UTR_sequence = str(df['5UTR_sequence'].values[0]).upper()
    insert_position = int(df['insert_position'].values[0]) - 1
    pre_insertion = UTR_sequence[:insert_position+1][-5:]
    post_insertion = UTR_sequence[insert_position+1:][:3]
    return pre_insertion, post_insertion

def analyse_inserts(df, pre_sequence, post_sequence):
    df = df.copy()
    UTR_sequence = str(df['5UTR_sequence'].values[0]).upper()
    insert_position = int(df['insert_position'].values[0]) - 1
    post_sequence_in_context = UTR_sequence[insert_position + 4:]

    inserts = generate_dna_sequences(6)
    rows = []
    for insert in inserts:
        insert_in_context = pre_sequence + insert + post_sequence
        start_pos = insert_in_context.find('ATG')
        if start_pos != -1:
            # sequence starting at ATG, extended with downstream UTR
            seq_from_start = (insert_in_context + post_sequence_in_context)[start_pos:]
            stop_rel = detect_stop_codons(seq_from_start)
            if stop_rel is not False:
                type_ = 'uORF'
                length_uORF = stop_rel + 3
                uORF_start = insert_position + start_pos - 5
                uORF_end = uORF_start + length_uORF
            else:
                # no stop; classify oORF by frame
                type_ = 'oORF_inframe' if (len(seq_from_start) % 3 == 0) else 'oORF_outframe'
                length_uORF = 'NaN'
                uORF_start = 'NaN'
                uORF_end = 'NaN'

            # Kozak check with bounds
            up3_ok = (start_pos - 3) >= 0 and insert_in_context[start_pos-3] in ('A', 'G')
            dn3_ok = (start_pos + 3) < len(insert_in_context) and insert_in_context[start_pos+3] == 'G'
            if up3_ok and dn3_ok:
                kozak = 'strong'
            elif up3_ok or dn3_ok:
                kozak = 'moderate'
            else:
                kozak = 'weak'
        else:
            type_ = 'NaN'
            kozak = 'NaN'
            length_uORF = 'NaN'
            uORF_start = 'NaN'
            uORF_end = 'NaN'

        rows.append({
            'sequence': insert,
            'uORF': type_,
            'uORF_length': length_uORF,
            'uORF_start': uORF_start,
            'uORF_end': uORF_end,
            'kozak_strength': kozak,
        })
    return pd.DataFrame(rows)

def analyse_UTR(gene, UTR_info_file):
    # read and normalize
    UTR_info = pd.read_csv(UTR_info_file)
    if 'gene' not in UTR_info.columns:
        raise ValueError(f"'gene' column not found in {UTR_info_file}")
    UTR_info['gene'] = UTR_info['gene'].astype(str).str.strip()
    sub = UTR_info[UTR_info['gene'] == gene]
    if sub.empty:
        raise ValueError(f"Gene '{gene}' not found in {UTR_info_file}")
    pre_sequence, post_sequence = get_insertion_context(sub)
    return analyse_inserts(sub, pre_sequence, post_sequence)

# ----------------------------- run ------------------------------
uORF_analysis = analyse_UTR(gene, UTR_info_file)

# incorporate uORF classification
scores = pd.read_csv(f"{save_folder}{gene}_scores_variant_features.csv")
merged_scores = pd.merge(scores, uORF_analysis, on='sequence', how='left')
merged_scores.to_csv(f"{save_folder}{gene}_scores_variant_features.csv", index=False)
