from Bio import SeqIO
from collections import defaultdict
import json

ref_seq = SeqIO.read("ecoli.fasta", "fasta")
w = 31
k = 7
L = len(ref_seq)
print("L is:", L)

'''
USING LISTS
'''

kmer_index = defaultdict(list)

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

	if i % 100000 == 0:
		print("The current sequence is: " + sub_f.seq)
		print("\tIts minimizer is: " + curr_min.seq)
		print("the index is:", index + i)
		print("Percentage of whole code: ", (i / (L - w + 1)) * 100, "%")

	kmer_index[str(curr_min.seq)].append(i + index)

kmer_index = {k: list(set(v)) for k, v in kmer_index.items()}

with open('minimized_ecoli_w31_k7.json', 'w') as f:
	json.dump(kmer_index, f)

'''
USING SETS
'''

# kmer_index = defaultdict(set)
#
# for i in range(0, L - w + 1):
# 	char = "T"
# 	curr_min = char.rjust(k, "T")
# 	sub_f = ref_seq[i:i + w]
#
# 	index = 0
# 	for j in range(0, w - k + 1):
# 		sub2 = sub_f[j:j + k]
# 		if type(curr_min) == str:
# 			if sub2.seq < curr_min:
# 				curr_min = sub2
# 				index = j
# 		else:
# 			if sub2.seq < curr_min.seq:
# 				curr_min = sub2
# 				index = j
#
# 	kmer_index[str(curr_min.seq)].add(i + index)
# 	print("The current sequence is: " + sub_f.seq)
# 	print("\tIts minimizer is: " + curr_min.seq)
# 	print("the index is:", index + i)
#
# kmer_index = {k: list(v) for k, v in kmer_index.items()}
#
# with open('minimized_lambda_w31_k7.json', 'w') as f:
# 	json.dump(kmer_index, f)
