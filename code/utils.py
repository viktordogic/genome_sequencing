import csv
import json

with open('minimized_ecoli_w75_k15.json') as f:
    referent_genome_json = json.load(f)

with open('double.json') as g:
    genome_sequence_json_forward = json.load(g)

with open('double_reversed.json') as h:
    genome_sequence_json_reversed = json.load(h)


def export_minimizers_to_json(ref, forward, reverse):
    with open('minimized_ecoli_w75_k15.json', 'w') as f:
        json.dump(ref, f)
    with open('minimized_lambda_reads_w75_k15.json', 'w') as g:
        json.dump(forward, g)
    with open('minimized_lambda_reads_w75_k15.json', 'w') as h:
        json.dump(reverse, h)


def generate_chars(length):
    return ["."] * length


def generate_mutation_report(alignment, report_lista, t_begin, string_checker):
    aligned_ref = alignment[0]
    aligned_mut = alignment[1]

    for i in range(len(aligned_ref)):
        if string_checker[i + t_begin] == '.':
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
        headers = next(csv_reader)
        sorted_rows = sorted(csv_reader, key=lambda row: int(row[1]))

    with open('sorted_' + filename, 'w', newline='') as file:
        csv_writer = csv.writer(file)
        csv_writer.writerow(headers)
        csv_writer.writerows(sorted_rows)