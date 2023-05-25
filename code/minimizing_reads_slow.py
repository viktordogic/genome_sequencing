from Bio import SeqIO
from collections import defaultdict
import json

list_of_reads = list(SeqIO.parse("lambda_simulated_reads.fasta", "fasta"))
w = 75
k = 15

'''
USING SETS
'''

kmer_index = defaultdict(lambda: defaultdict(set))

for seq_index, ref_seq in enumerate(list_of_reads):
	L = len(ref_seq)
	for i in range(0, L - w + 1):
		char = "T"
		curr_min = char.rjust(k, "T")
		sub_f = ref_seq[i:i + w]

		index = 0
		for j in range(0, w - k + 1):
			sub2 = sub_f[j:j + k]
			if type(curr_min) == str:
				if sub2.seq < curr_min:
					curr_min = sub2
					index = j
			else:
				if sub2.seq < curr_min.seq:
					curr_min = sub2
					index = j

		kmer_index[seq_index][str(curr_min.seq)].add(i + index)
	print("The current read is", seq_index)

kmer_index = {k: {k2: list(v2) for k2, v2 in v.items()} for k, v in kmer_index.items()}

with open('minimized_lambda_reads_w75_k15.json', 'w') as f:
	json.dump(kmer_index, f)
