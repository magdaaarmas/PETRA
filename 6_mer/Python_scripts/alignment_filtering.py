import os
import sys
import pandas as pd


working_dir = os.getcwd()
file_locations=sys.argv[1]
min_alignment_score = sys.argv[2]
min_alignment_score=int(min_alignment_score)

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
            if line_count < 0:
                line_count += 1
                continue
            else:
                line_count += 1
                sam_data = line.strip().split('\t')
                # cigar string
                cigar = sam_data[5]
                # this is the alignment score -- not necessarily useful unless mispriming common
                AS = sam_data[11]
                # sequence of read
                seq = sam_data[9]
                # filter based on the alignment score -- good alignments should be over ~500
                index_of_score = AS.find('i:') + 2
                # this is the actual score
                score = float(AS[index_of_score:])
                # hard coded alignment score cuttoff (using needleall with 10 gap open and 0.5 gap extension -- 100 is fairly arbitrary -- very low to be inclusive. could raise to only look at better reads)
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
    specific_directory=f"{working_dir}/aligned_reads/{gene}_aligned/"
    for file in os.listdir(specific_directory):
        if file.endswith(".sam"):
            with open(specific_directory+file, 'r') as sam_file:
                sam_to_variant_cigar_files(sam_file, min_alignment_score, specific_directory)
