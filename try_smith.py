from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment

# Load reference genome
ref_genome_file = "lambda.fasta"
ref_genome = next(SeqIO.parse(ref_genome_file, "fasta"))

# Load query sequences
query_file = "lambda_simulated_reads.fasta"
query_seqs = list(SeqIO.parse(query_file, "fasta"))

# Set alignment parameters
match_score = 2
mismatch_score = -1
gap_open_penalty = -2
gap_extend_penalty = -1

# Perform sequence alignment for each query sequence
for query_seq in query_seqs:
    alignments = pairwise2.align.localms(ref_genome.seq, query_seq.seq, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty)
    best_alignment = alignments[0]
    print(format_alignment(*best_alignment))  # print alignment for debugging purposes

    # Create new FASTA record for aligned sequence
    aligned_seq_record = query_seq.copy()
    aligned_seq_record.seq = best_alignment[1]
    aligned_seq_record.description = "Aligned to reference genome (score={})".format(best_alignment[2])

    # Write aligned sequence to new FASTA file
    with open("aligned_sequences.fasta", "a") as outfile:
        SeqIO.write(aligned_seq_record, outfile, "fasta")