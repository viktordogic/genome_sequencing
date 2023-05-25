from Bio import SeqIO
import json

with open('minimized_lambda_w75_k15.json') as f:
    referent_genome_json = json.load(f)

with open('minimized_lambda_reads_w75_k15.json') as g:
    genome_sequence_json_forward = json.load(g)

with open('minimized_lambda_reads_w75_k15_reversed.json') as h:
    genome_sequence_json_reversed = json.load(h)

referent_genome_io = SeqIO.read("lambda.fasta", "fasta")


def convert_to_string(char_list):
    string = ''.join(char_list)
    return string

def count_chars(char_list):
    count = 0
    for element in char_list:
        if element != '-':
            count += 1
    return count


ref_genome = ['-' for _ in range(len(referent_genome_io))]

str_ref_genome = convert_to_string(ref_genome)

def recreate_genome(ref_genome_minimizers, seq_minimizers, kmer_length):
    for seq_idx, seq_minimizers in seq_minimizers.items():
        for kmer, indices in seq_minimizers.items():
            if kmer in ref_genome_minimizers:
                ref_idx = ref_genome_minimizers[kmer][0]
                for seq_idx in indices:
                    ref_genome[ref_idx:ref_idx + kmer_length] = kmer
    return ''.join(ref_genome)

# recreate_genome(referent_genome_json, genome_sequence_json_forward, 15)
# recreate_genome(referent_genome_json, genome_sequence_json_reversed, 15)

print(count_chars(ref_genome))
