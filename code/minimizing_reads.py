from Bio import SeqIO
from collections import defaultdict, deque
import json
import time

from utils import sequences_io_forward, sequences_io_reversed


def create_minimizer_sequences(fasta_file):
    w = 75
    k = 15
    kmer_index = defaultdict(lambda: defaultdict(list))

    for seq_index, ref_seq in enumerate(fasta_file):
        L = len(ref_seq)
        min_kmer_queue = deque()  # Queue to store the minimum kmers in the current window
        min_kmer = None  # Current minimum kmer

        for i in range(L - k + 1):  # Loop over the kmers in the sequence
            while min_kmer_queue and min_kmer_queue[0][1] < i - w + k:  # Remove kmers outside of the current window
                min_kmer_queue.popleft()

            curr_kmer = ref_seq[i:i + k]
            while min_kmer_queue and min_kmer_queue[-1][0].seq >= curr_kmer.seq:  # Remove kmers larger than the current
                min_kmer_queue.pop()

            min_kmer_queue.append((curr_kmer, i))  # Add the current kmer to the queue
            if i >= w - k:  # Starting from the first complete window
                min_kmer = min_kmer_queue[0]  # The tuple at the front of the queue is the minimum in the current window
                kmer_index[seq_index][str(min_kmer[0].seq)].append(min_kmer[1])

        print("The current read is", seq_index)

    return {int(k): {k2: list(set(v2)) for k2, v2 in v.items()} for k, v in kmer_index.items()}


# forward_minimizers = create_minimizer_sequences(sequences_io_forward)
backward_minimizers = create_minimizer_sequences(sequences_io_reversed)

# with open('minimized_lambda_reads_w75_k15.json', 'w') as f:
#     json.dump(forward_minimizers, f)

with open('minimized_ecoli_reads_w75_k15_reversed.json', 'w') as f:
    json.dump(backward_minimizers, f)
