from Bio import SeqIO
import json
from mapping_verification import global_alignment
from mapping_verification_2 import create_cigar_string


def replace_substring(old_string, substring, index):
    return old_string[:index] + substring + old_string[index + len(substring):]


def sort_by_index(lst, index):
    return sorted(lst, key=lambda x: x[index])


def count_same_chars(str1, str2):
    count = 0
    for char1, char2 in zip(str1, str2):
        if char1 == char2:
            count += 1
    return count


with open('minimized_lambda_w75_k15.json') as f:
    referent_genome_json = json.load(f)

with open('minimized_lambda_reads_w75_k15.json') as g:
    genome_sequence_json_forward = json.load(g)

with open('minimized_lambda_reads_w75_k15_reversed.json') as h:
    genome_sequence_json_reversed = json.load(h)

referent_genome_io = SeqIO.read("lambda.fasta", "fasta")
sequences_io_forward = list(SeqIO.parse("lambda_simulated_reads.fasta", "fasta"))
sequences_io_reversed = list()
for sequence in sequences_io_forward:
    sequences_io_reversed.append(sequence.reverse_complement())

char = "-"
new_referent_genome = char.rjust(len(referent_genome_io), "-")

list_of_best_finds = []
for seq_index, seq_seq in genome_sequence_json_forward.items():
    ndata = []
    for k, v in referent_genome_json.items():
        if k in seq_seq.keys():
            is_contained = False
            if len(ndata) == 0:
                ndata = [[v[0], seq_seq[k][0], 1, seq_index, 'f']]
                continue
            for element in ndata:
                diff_ref = v[0] - element[0]
                diff_seq = seq_seq[k][0] - element[1]
                if diff_ref == diff_seq:
                    element[2] += 1
                    is_contained = True
                    break
            if not is_contained:
                ndata.append([v[0], seq_seq[k][0], 1, seq_index, 'f'])

    if len(ndata) != 0:
        best_find = ndata[0]

        for other_elements in range(1, len(ndata)):
            if ndata[other_elements][2] > best_find[2]:
                best_find = ndata[other_elements]

        list_of_best_finds.append(best_find)

for seq_index, seq_seq in genome_sequence_json_reversed.items():
    ndata = []
    for k, v in referent_genome_json.items():
        if k in seq_seq.keys():
            is_contained = False
            if len(ndata) == 0:
                ndata = [[v[0], seq_seq[k][0], 1, seq_index, 'r']]
                continue
            for element in ndata:
                diff_ref = v[0] - element[0]
                diff_seq = seq_seq[k][0] - element[1]
                if diff_ref == diff_seq:
                    element[2] += 1
                    is_contained = True
                    break
            if not is_contained:
                ndata.append([v[0], seq_seq[k][0], 1, seq_index, 'r'])

    if len(ndata) != 0:
        best_find = ndata[0]

        for other_elements in range(1, len(ndata)):
            if ndata[other_elements][2] > best_find[2]:
                best_find = ndata[other_elements]

        list_of_best_finds.append(best_find)

for index in range(len(list_of_best_finds)):
    list_of_best_finds[index][3] = int(list_of_best_finds[index][3])

list_of_best_finds = sort_by_index(list_of_best_finds, 1)

indeksi = []
for index in range(len(list_of_best_finds)):
    if list_of_best_finds[index][1] > list_of_best_finds[index][0]:
        indeksi.append(index)

for index in indeksi:
    del list_of_best_finds[index]

del list_of_best_finds[48]
del list_of_best_finds[224]

print(list_of_best_finds[-1])
print(sequences_io_forward[int(list_of_best_finds[-1][3])][13406:].seq)

# for bla in list_of_best_finds:
#     print(bla)

print(referent_genome_io.seq[43213:])

# for indeks, tapl in enumerate(list_of_best_finds):
#     if tapl[4] == 'f':
#         new_referent_genome = replace_substring(new_referent_genome, sequences_io_forward[int(tapl[3])].seq,
#                                                 tapl[0] - tapl[1])
#     else:
#         new_referent_genome = replace_substring(new_referent_genome, sequences_io_reversed[int(tapl[3])].seq,
#                                                 tapl[0] - tapl[1])

