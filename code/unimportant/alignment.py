from Bio import SeqIO
import json
from utils import *


list_of_best_finds = []
for seq_index, seq_seq in genome_sequence_json_forward.items():
    ndata = []
    for k, v in referent_genome_json.items():
        if k in seq_seq.keys():
            is_contained = False
            if len(ndata) == 0:
                ndata = [[v[0], seq_seq[k][0], 1, seq_index, 'f', seq_seq[k][0]]]
                continue
            for element in ndata:
                diff_ref = v[0] - element[0]
                diff_seq = seq_seq[k][0] - element[1]
                if diff_ref == diff_seq:
                    element[2] += 1
                    element[5] = seq_seq[k][0]
                    is_contained = True
                    break
            if not is_contained:
                ndata.append([v[0], seq_seq[k][0], 1, seq_index, 'f', seq_seq[k][0]])

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
                ndata = [[v[0], seq_seq[k][0], 1, seq_index, 'r', seq_seq[k][0]]]
                continue
            for element in ndata:
                diff_ref = v[0] - element[0]
                diff_seq = seq_seq[k][0] - element[1]
                if diff_ref == diff_seq:
                    element[2] += 1
                    element[5] = seq_seq[k][0]
                    is_contained = True
                    break
            if not is_contained:
                ndata.append([v[0], seq_seq[k][0], 1, seq_index, 'r', seq_seq[k][0]])

    if len(ndata) != 0:
        best_find = ndata[0]

        for other_elements in range(1, len(ndata)):
            if ndata[other_elements][2] > best_find[2]:
                best_find = ndata[other_elements]

        list_of_best_finds.append(best_find)

for index in range(len(list_of_best_finds)):
    list_of_best_finds[index][3] = int(list_of_best_finds[index][3])

list_of_best_finds = sort_by_index(list_of_best_finds, 2)

list_of_best_finds_forward = []
list_of_best_finds_reversed = []

for lista in list_of_best_finds:
    if lista[4] == 'f':
        list_of_best_finds_forward.append(lista)
    else:
        list_of_best_finds_reversed.append(lista)
