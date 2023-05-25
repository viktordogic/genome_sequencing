from Bio import SeqIO
from collections import defaultdict, deque
import json

ref_seq = SeqIO.read("lambda.fasta", "fasta")
w = 75
k = 15
L = len(ref_seq)
print("L is:", L)

kmer_index = defaultdict(list)
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
		kmer_index[str(min_kmer[0].seq)].append(min_kmer[1])

kmer_index = {k: list(set(v)) for k, v in kmer_index.items()}

with open('minimized_lambda_w75_k15.json', 'w') as f:
	json.dump(kmer_index, f)
