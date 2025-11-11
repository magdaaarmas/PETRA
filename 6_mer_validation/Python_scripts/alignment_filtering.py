import os
import sys
import pandas as pd


working_dir = os.getcwd()
file_locations=sys.argv[1]
min_alignment_score = sys.argv[2]
min_alignment_score=int(min_alignment_score)
experiment_folder=sys.argv[3]


def sam_to_variant_cigar_files(sam_file, min_alignment_score, gene_directory):
    line_count = -2
    reads_not_aligning = 0
    reads_aligning = 0
    index_of_first_dot = sam_file.name.find('.')
    last_slash = sam_file.name.rfind("/")
    sample = sam_file.name[last_slash+1:index_of_first_dot]
    seq_output_dir=gene_directory+str(min_alignment_score)+'_aligned_sequences/'
    cigar_output_dir=gene_directory+str(min_alignment_score)+'_aligned_cigars/'
    os.makedirs(seq_output_dir, exist_ok=True)
    os.makedirs(cigar_output_dir, exist_ok=True)
    seq_output_file = seq_output_dir+sample + '_aligned_sequences.txt'
    cigar_output_file = cigar_output_dir+sample + '_aligned_cigar.txt'
    alignment_info_file = gene_directory+sample + '_'+str(min_alignment_score)+'_alignment_info.txt'
    

    with open(seq_output_file, 'w') as seq_outfile, open(cigar_output_file, 'w') as cigar_outfile, open(alignment_info_file,'w') as alignment_info_outfile:
        for line in sam_file:
            # skip SAM header lines (start with '@') and empty lines
            if not line.strip() or line.startswith('@'):
                continue

            sam_data = line.rstrip('\n').split('\t')

            # need at least the 11 mandatory SAM fields (0..10) + typically optional fields
            if len(sam_data) < 11:
                # log malformed/short line and skip it
                # you could write it to a debug file instead of printing
                print(f"[WARN] skipping short/malformed SAM line: {line[:80]!r}")
                continue

            # safe access to standard fields
            # QNAME(0), FLAG(1), RNAME(2), POS(3), MAPQ(4), CIGAR(5), ...
            cigar = sam_data[5]

            # OPTIONAL FIELDS: find the AS tag anywhere in sam_data[11:]
            AS_field = None
            for f in sam_data[11:]:
                # accept formats like "AS:i:123" or "AS:i:-45"
                if f.startswith('AS:'):
                    AS_field = f
                    break

            if AS_field is None:
                # no AS tag found — either skip, or treat as low-score
                print(f"[INFO] no AS tag for read (treating as score 0): {sam_data[0]}")
                score = 0.0
            else:
                # parse the numeric score robustly
                # AS_field is like 'AS:i:123' -> take last part after ':'
                try:
                    score = float(AS_field.split(':')[-1])
                except ValueError:
                    print(f"[WARN] could not parse AS field {AS_field!r} for read {sam_data[0]}; skipping")
                    continue

            seq = sam_data[9]

            if score > min_alignment_score:
                seq_outfile.write(seq + '\n')
                cigar_outfile.write(cigar + '\n')
                reads_aligning += 1
            else:
                reads_not_aligning += 1


        alignment_info_outfile.write(f"Reads aligning: {reads_aligning}\n")
        alignment_info_outfile.write(f"Reads not aligning: {reads_not_aligning}\n")

file_paths=pd.read_csv(file_locations)
genes=file_paths['gene'].tolist()

for gene in genes:
    specific_directory=f"{working_dir}/{experiment_folder}/aligned_reads/{gene}_aligned/"
    for file in os.listdir(specific_directory):
        if file.endswith(".sam"):
            with open(specific_directory+file, 'r') as sam_file:
                sam_to_variant_cigar_files(sam_file, min_alignment_score, specific_directory)
