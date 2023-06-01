import csv
import json
from Bio import SeqIO

referent_genome_io = SeqIO.read("ecoli.fasta", "fasta")
sequences_io_forward = list(SeqIO.parse("ecoli_simulated_reads.fasta", "fasta"))
sequences_io_reversed = list()
for sequence in sequences_io_forward:
    sequences_io_reversed.append(sequence.reverse_complement())


with open('minimized_ecoli_w75_k15.json') as f:
    referent_genome_json = json.load(f)

with open('double.json') as g: # 'minimized_ecoli_reads_w75_k15.json'
    genome_sequence_json_forward = json.load(g) # 'filtered_ecoli_reads_w75_k15.json'

with open('double_reversed.json') as h: # minimized_lambda_reads_w75_k15_reversed.json
    genome_sequence_json_reversed = json.load(h) # 'filtered_lambda_reads_w75_k15_reversed.json'


def export_minimizers_to_json(ref, forward, reverse):
    with open('minimized_ecoli_w75_k15.json', 'w') as f:
        json.dump(ref, f)
    with open('minimized_lambda_reads_w75_k15.json', 'w') as g:
        json.dump(forward, g)
    with open('minimized_lambda_reads_w75_k15.json', 'w') as h:
        json.dump(reverse, h)


def generate_chars(length):
    return ["."] * length


index_checking = generate_chars(len(referent_genome_io.seq))


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


def convert_to_string(char_list):
    string = ''.join(char_list)
    return string


def count_chars(char_list):
    count = 0
    for element in char_list:
        if element != '-':
            count += 1
    return count


def create_simple_cigar_string(referent, sequence):
    i = 0
    j = 0
    counter = 0

    while i < len(referent) and j < len(sequence):
        if referent[i] == sequence[j]:
            counter += 1
        i += 1
        j += 1

    return counter


def create_cigar_string(referent, sequence):
    cigar_string = ''
    i = 0
    j = 0
    counter = 0

    while i < len(referent) and j < len(sequence):
        if referent[i] == sequence[j]:
            cigar_string += 'M'
            i += 1
            j += 1
            counter += 1
        elif i < len(referent) - 1 and referent[i+1] == sequence[j]:
            cigar_string += 'D'
            i += 2
            j += 1
        elif j < len(sequence) - 1 and sequence[j+1] == referent[i]:
            cigar_string += 'I'
            i += 1
            j += 2
        else:
            cigar_string += 'X'
            i += 1
            j += 1

    # In case sequence or referent still has unprocessed elements
    while i < len(referent):
        cigar_string += 'D'
        i += 1
    while j < len(sequence):
        cigar_string += 'I'
        j += 1

    return cigar_string


def generate_mutation_report(alignment, report_lista, t_begin):
    aligned_ref = alignment[0]
    aligned_mut = alignment[1]

    for i in range(len(aligned_ref)):
        if index_checking[i + t_begin] == '.':
            if aligned_ref[i] == '-':  # This is an insertion
                report_lista.append(['I', i + t_begin, aligned_mut[i]])
            elif aligned_mut[i] == '-':  # This is a deletion
                report_lista.append(['D', i + t_begin, '-'])
            elif aligned_ref[i] != aligned_mut[i]:  # This is a substitution
                report_lista.append(['X', i + t_begin, aligned_mut[i]])


def export_to_csv(data, filename):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(data)


def sort_csv_by_second_column(filename):
    with open(filename, 'r') as file:
        csv_reader = csv.reader(file)
        headers = next(csv_reader)  # Read the header row
        sorted_rows = sorted(csv_reader, key=lambda row: int(row[1]))  # Sort rows by the second column

    # Write the sorted rows to a new file
    with open('sorted_' + filename, 'w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow(headers)  # Write the header row
        csv_writer.writerows(sorted_rows)  # Write the sorted rows