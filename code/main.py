from minimizing import create_sequence_minimizer
from minimizing_reads import create_minimizer_sequences
from mapping import create_alignments
from filtering_minimizer import remove_recurring_keys
from utils import *


if __name__ == '__main__':
    print("Starting the program...")
    print("Reading the genomes...")
    referent_genome_io = SeqIO.read("ecoli.fasta", "fasta")
    sequences_io_forward = list(SeqIO.parse("ecoli_simulated_reads.fasta", "fasta"))
    sequences_io_reversed = list()
    for sequence in sequences_io_forward:
        sequences_io_reversed.append(sequence.reverse_complement())

    print("Creating reference genome minimizers...")
    referent_json = create_sequence_minimizer(referent_genome_io)

    print("Creating minimizers for forward reads...")
    forward_json = create_minimizer_sequences(sequences_io_forward)

    print("Creating minimizers for reversed reads...")
    reversed_json = create_minimizer_sequences(sequences_io_reversed)

    print("Filtering the minimizers...")
    remove_recurring_keys(forward_json)
    remove_recurring_keys(reversed_json)

    report_list = [["Mutacija", "Pozicija", "Vrijednost"]]

    print("Creating alignments...")
    alignment_list = create_alignments(referent_json, forward_json, reversed_json)

    for alignment in alignment_list:
        if alignment[0] < 50:
            break
        generate_mutation_report(alignment[1], report_list, alignment[2])
        index_checking[alignment[2]:alignment[3]] = ['x'] * (alignment[3] - alignment[2])
        # print(sum([1 for x in index_checking if x == 'x']))

    export_to_csv(report_list, 'output_file.csv')
    sort_csv_by_second_column('output_file.csv')
