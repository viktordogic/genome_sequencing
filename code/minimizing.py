from collections import defaultdict, deque


def create_sequence_minimizer(ref_gen, w=75, k=15):
    L = len(ref_gen)
    print("L is:", L)

    kmer_index = defaultdict(list)
    min_kmer_queue = deque()
    min_kmer = None

    for i in range(L - k + 1):
        while min_kmer_queue and min_kmer_queue[0][1] < i - w + k:
            min_kmer_queue.popleft()

        curr_kmer = ref_gen[i:i + k]
        while min_kmer_queue and min_kmer_queue[-1][0].seq >= curr_kmer.seq:
            min_kmer_queue.pop()

        min_kmer_queue.append((curr_kmer, i))
        if i >= w - k:
            min_kmer = min_kmer_queue[0]
            kmer_index[str(min_kmer[0].seq)].append(min_kmer[1])

    return {k: list(set(v)) for k, v in kmer_index.items()}
