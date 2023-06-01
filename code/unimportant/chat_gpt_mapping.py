# from alignment import generate_mutation_report
# from utils import *
#
# ref_genome = ['-' for _ in range(len(referent_genome_io))]
# str_ref_genome = convert_to_string(ref_genome)
#
# for k, v in genome_sequence_json_forward.items():
#     for k2, v2 in v.items():
#         if k2 in referent_genome_json.keys():
#             print(k, k2, v2, referent_genome_json[k2])
#
# from Bio import SeqIO, Align
# from utils import *
#
# aligner = Align.PairwiseAligner()
# aligner.match_score = 1
# aligner.mode = 'global'
# alignments = aligner.align(sequences_io_forward[212].seq[1200:2200], referent_genome_io.seq[2100:3100])
# alignment = alignments[0]
# generate_mutation_report(alignments, 'cum.csv')
# print(alignment)
# print("Score = %.1f:" % alignment.score)
import csv
import json
import math
from collections import defaultdict
from utils import *

from typing import Tuple, List

from Bio.Align import PairwiseAligner

from copyminimizers import Minimize

kmer_len = 15
EPSILON_BAND = 10
cigar_flag = 1
Match = Tuple[bool, int, int]

# def generate_mutation_report(alignment, filename, t_begin, repeat):
#     with open(filename, 'a', newline='') as file:
#         writer = csv.writer(file)
#         if repeat:
#             writer.writerow(["Mutacija", "Pozicija", "Vrijednost"])
#
#         aligned_ref = alignment[0]
#         aligned_mut = alignment[1]
#
#         for i in range(len(aligned_ref)):
#             if aligned_ref[i] == '-':  # This is an insertion
#                 writer.writerow(['I', i + t_begin, aligned_mut[i]])
#             elif aligned_mut[i] == '-':  # This is a deletion
#                 writer.writerow(['D', i + t_begin, '-'])
#             elif aligned_ref[i] != aligned_mut[i]:  # This is a substitution
#                 writer.writerow(['X', i + t_begin, aligned_mut[i]])


report_list = [["Mutacija", "Pozicija", "Vrijednost"]]


def generate_chars(length):
    return ["."] * length


index_checking = generate_chars(len(referent_genome_io.seq))


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


def longest_increasing_subsequence(entries):
    # Sort entries based on position in reference genome
    entries.sort(key=lambda x: x[0])

    n = len(entries)

    # Memoization arrays
    memo = [1] * n  # Lengths of longest increasing subsequences
    pred = [None] * n  # Predecessors in longest increasing subsequences

    # Fill memo and pred
    for i in range(n):
        for j in range(i):
            if entries[j][0] < entries[i][0] and memo[i] < memo[j] + 1:
                memo[i] = memo[j] + 1
                pred[i] = j

    # Find end of longest increasing subsequence
    end = memo.index(max(memo))

    # Reconstruct longest increasing subsequence
    lis = []
    while end is not None:
        lis.append(entries[end])
        end = pred[end]
    lis.reverse()

    return lis


def GetCeilIndex(curr_cluster, T, l, r, key):
    while r - l > 1:
        m = l + (r - l) // 2
        if curr_cluster[T[m]][2] >= key:
            r = m
        else:
            l = m
    return r


def longestIncreasingSubsequence(curr_cluster):
    n = len(curr_cluster)
    tailIndices = [0] * n
    prevIndices = [-1] * n

    len_ = 1
    for i in range(1, n):
        if curr_cluster[i][2] < curr_cluster[tailIndices[0]][2]:
            tailIndices[0] = i
        elif curr_cluster[i][2] > curr_cluster[tailIndices[len_ - 1]][2]:
            prevIndices[i] = tailIndices[len_ - 1]
            tailIndices[len_] = i
            len_ += 1
        else:
            pos = GetCeilIndex(curr_cluster, tailIndices, -1, len_ - 1, curr_cluster[i][2])
            if pos == 0:
                prevIndices[i] = -1
            else:
                prevIndices[i] = tailIndices[pos - 1]
            tailIndices[pos] = i

    lisIndexes = []
    i = tailIndices[len_ - 1]
    while i >= 0:
        lisIndexes.append(i)
        i = prevIndices[i]

    lisIndexes.reverse()
    result = [curr_cluster[i] for i in lisIndexes]

    return result


def getQueryLoc(x: Match) -> int:
    if x[0]:
        return x[1] - x[2]
    else:
        return x[1] + x[2]


def getRefLoc(x: Match) -> int:
    return x[2]


def kmerOverlap(prev: int, curr: int) -> int:
    if abs(curr - prev) < kmer_len:
        return abs(curr - prev)
    else:
        return kmer_len


def col10Approx(lis_of_cluster: List[Match]) -> int:
    approx_val_q = kmer_len
    approx_val_r = kmer_len
    prev_q = getQueryLoc(lis_of_cluster[0])
    prev_r = getRefLoc(lis_of_cluster[0])
    first = True
    for match in lis_of_cluster:
        if first:
            first = False
            continue
        curr_q = getQueryLoc(match)
        curr_r = getRefLoc(match)
        approx_val_q += kmerOverlap(prev_q, curr_q)
        approx_val_r += kmerOverlap(prev_r, curr_r)
        prev_q = curr_q
        prev_r = curr_r
    return min(approx_val_r, approx_val_q)


def splitIntoClustersAndFindBest(matches):
    curr_cluster = [matches[0]]
    max_col10 = 0
    match_cluster = []

    for i in range(1, len(matches) + 1):
        if (i == len(matches) or
                matches[i][0] != matches[i - 1][0] or
                matches[i][1] - matches[i - 1][1] >= EPSILON_BAND):

            lis_of_cluster = longestIncreasingSubsequence(curr_cluster)
            curr_col10 = col10Approx(lis_of_cluster)
            if curr_col10 > max_col10 and len(lis_of_cluster) >= 4:
                max_col10 = curr_col10
                match_cluster = lis_of_cluster.copy()

            curr_cluster = []

        if i < len(matches):  # to avoid index out of range error
            curr_cluster.append(matches[i])

    return match_cluster


def bestMatchCluster(fragment_index, reference_index):
    matches = []
    for f_entry in fragment_index:
        if f_entry in reference_index:
            for f_location in fragment_index[f_entry]:
                for r_location in reference_index[f_entry]:
                    diff_strand = f_location[1] ^ r_location[1]
                    if diff_strand:
                        relative_position = f_location[0] + r_location[0]
                    else:
                        relative_position = f_location[0] - r_location[0]
                    match = (diff_strand, relative_position, r_location[0])
                    matches.append(match)

    if matches:
        # print(matches)
        matches.sort(key=lambda x: (not x[0], x[1], x[2]))
        return splitIntoClustersAndFindBest(matches)


def makeIndex(sequence):
    index = defaultdict(list)
    minimizers = Minimize(sequence.seq, len(sequence.seq))
    for minimizer in minimizers:
        index[minimizer[0]].append((minimizer[1], minimizer[2]))
    print(index)
    return index


def make_cleaned_reference_index(sequence, frequency):
    print("Making reference index...")
    index = makeIndex(sequence)
    minimizers_occurance = []
    num_of_singletons = 0

    for key, value in index.items():
        if len(value) == 1:
            num_of_singletons += 1
        minimizers_occurance.append((len(value), key))

    print("Sorting reference minimizers by occurance...")
    minimizers_occurance.sort(key=lambda x: x[0], reverse=True)
    minimizers_to_skip = math.ceil(len(index) * frequency)

    if minimizers_to_skip >= len(index):
        minimizers_to_skip = len(index) - 1

    print(f"Number of distinct minimizers: {len(index)}")
    print(f"Fraction of singletons: {num_of_singletons / len(index)}")
    print(
        f"Number of occurrences of the most frequent minimizer with top {frequency * 100}% most frequent ignored: {minimizers_occurance[minimizers_to_skip][0]}")

    print("Removing too frequent minimizers from reference index...")
    for i in range(minimizers_to_skip):
        del index[minimizers_occurance[i][1]]

    print("Done.\n")

    return index


def getRelStrandString(x):
    return '+' if x[0] else '-'


def reverseComplement(seq):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    result = ""
    for base in reversed(seq):
        result += complement_dict.get(base, base)
    return result


def getPaf(fragment, reference, match_cluster):
    if match_cluster[0][0]:
        q_begin = int(getQueryLoc(match_cluster[-1]))
        q_end = int(getQueryLoc(match_cluster[0])) + kmer_len - 1
        is_reverse_complement = True
    else:
        q_begin = int(getQueryLoc(match_cluster[0]))
        q_end = int(getQueryLoc(match_cluster[-1])) + kmer_len - 1
        is_reverse_complement = False
    t_begin = int(getRefLoc(match_cluster[0]))
    t_end = int(getRefLoc(match_cluster[-1])) + kmer_len - 1

    if cigar_flag:
        # try:
        if is_reverse_complement:
            rev_compl = reverseComplement(fragment.seq[q_begin:q_end])
            fragment_sequence = rev_compl
        else:
            fragment_sequence = fragment.seq[q_begin:q_end]

        aligner = PairwiseAligner()
        aligner.match_score = 1
        aligner.mismatch_score = -1
        aligner.mode = 'local'
        aligner.gap_score = -1
        # print("tbeing", t_begin)
        # print("tend", t_end)
        alignment = aligner.align(reference.seq[t_begin:t_end], fragment_sequence)[0]
        if alignment.score / (t_end - t_begin) > 0.75:
            print(
                f"FOR AN INTERVAL OF {t_begin} TO {t_end}, THE ALIGNMENT SCORE IS {(alignment.score / (t_end - t_begin)):0.2f}%")
            generate_mutation_report(alignment, report_list, t_begin)
            index_checking[t_begin:t_end] = ['x'] * (t_end - t_begin)
        else:
            print("NEIN")
        return
        # except Exception:
        #     return ""


ref_index = make_cleaned_reference_index(referent_genome_io, 0.001)
count = 0
count_2 = 0
for fragment in sequences_io_forward:
    frag_index = makeIndex(fragment)
    match_cluster = bestMatchCluster(frag_index, ref_index)
    if match_cluster is not None and len(match_cluster) != 0:
        getPaf(fragment, referent_genome_io, match_cluster)
        count += 1

export_to_csv(report_list, 'mutation_report.csv')
sort_csv_by_second_column('mutation_report.csv')

print(sum(1 for i in index_checking if i == 'x'))

print(index_checking)
