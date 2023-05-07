import Bio
from Bio import SeqIO
from collections import defaultdict

from Bio.Align import PairwiseAligner

# Load the reference genome
print("reading ecoli")
ref_seq = SeqIO.read("lambda.fasta", "fasta")

# Build a k-mer index for the reference genome
k = 40
kmer_index = defaultdict(list)
for i in range(len(ref_seq) - k + 1):
    kmer = str(ref_seq[i:i + k].seq)
    kmer_index[kmer].append(i)

# Load the reads
print("reading reads")
reads = list(SeqIO.parse("lambda_simulated_reads.fasta", "fasta"))

# # Map the reads to the reference genome using the k-mer index and Needleman-Wunsch algorithm
alignments = []
aligner = PairwiseAligner()
aligner.mode = 'local'
print("being alignment")
counter = 0
for read in reads:
    for i in range(len(read) - k + 1):
        print("start {} of read {}".format(i, len(read) - k + 1))
        print("we are doing that for {}".format(len(reads)))
        kmer = str(read[i:i + k].seq)
        if kmer in kmer_index:
            for ref_pos in kmer_index[kmer]:
                ref_subseq = str(ref_seq[ref_pos:ref_pos + len(read)].seq)
                alignment = aligner.align(ref_subseq, str(read.seq))[0]
                if alignment.score > len(read) * 0.9:
                    alignments.append((read.id, ref_pos, alignment))

# Identify mutations
print("Identify mutations")
mutations = []
for read_id, ref_pos, score in alignments:
    read_seq = str(reads[read_id].seq)
    ref_seq = str(ref_seq[ref_pos:ref_pos + len(read_seq)].seq)
    mutation_type = ""
    mutation_pos = -1
    mutation_length = 0
    for i in range(len(read_seq)):
        if read_seq[i] != ref_seq[i]:
            if read_seq[i:i + 3] == ref_seq[i:i + 3]:
                mutation_type = "substitution"
                mutation_pos = ref_pos + i
                mutation_length = 1
            elif read_seq[i:i + 3] == "-" * 3:
                mutation_type = "deletion"
                mutation_pos = ref_pos + i
                while i < len(read_seq) and read_seq[i] == "-":
                    mutation_length += 1
                    i += 1
            elif ref_seq[i:i + 3] == "-" * 3:
                mutation_type = "insertion"
                mutation_pos = ref_pos + i
                while i < len(read_seq) and ref_seq[i] == "-":
                    mutation_length += 1
                    i += 1
    if mutation_type:
        mutations.append((mutation_type, mutation_pos, mutation_length))

# Write the mutations to a CSV file
print("write")
print(mutations)
with open("mutations.csv", "w") as f:
    f.write("Mutation Type,Position,Length\n")
    for mutation in mutations:
        f.write("{},{},{}\n".format(mutation[0], mutation[1], mutation[2]))