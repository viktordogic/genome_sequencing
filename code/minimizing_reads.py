from Bio import SeqIO
from collections import defaultdict, deque
import json
import time

def create_minimizer_sequences(fasta_file):
    w = 75
    k = 15
    kmer_index = defaultdict(lambda: defaultdict(list))

    for seq_index, ref_seq in enumerate(fasta_file):
        L = len(ref_seq)
        min_kmer_queue = deque()
        min_kmer = None

        for i in range(L - k + 1):
            while min_kmer_queue and min_kmer_queue[0][1] < i - w + k:
                min_kmer_queue.popleft()

            curr_kmer = ref_seq[i:i + k]
            while min_kmer_queue and min_kmer_queue[-1][0].seq >= curr_kmer.seq:
                min_kmer_queue.pop()

            min_kmer_queue.append((curr_kmer, i))
            if i >= w - k:
                min_kmer = min_kmer_queue[0]
                kmer_index[seq_index][str(min_kmer[0].seq)].append(min_kmer[1])

        print("The current read is", seq_index)

    return {int(k): {k2: list(set(v2)) for k2, v2 in v.items()} for k, v in kmer_index.items()}
