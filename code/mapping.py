from Bio import SeqIO, Align


def create_alignments(ref_json, forward_json, reversed_json, ref_seq, forward_seq, reversed_seq):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'

    last_index = 0
    counter = 0
    list_of_alignments = list()
    for k_ref, v_ref in ref_json.items():
        counter += 1
        # if counter % 1000 == 0:
            # print("FEATURED NUMBER OF KEYS:", counter)
        find_info = list()
        found = False
        for k, v in forward_json.items():
            if k_ref in v.keys():
                find_info = v_ref[0], k, v[k_ref][0], 'f'
                found = True
        if not found:
            for k, v in reversed_json.items():
                if k_ref in v.keys():
                    find_info = v_ref[0], k, v[k_ref][0], 'r'

        if len(find_info) != 0:
            ref_pos, seq_nr, seq_pos, strand = find_info
            length_of_string = ref_pos - last_index
            seq_l = max(0, seq_pos - length_of_string)
            if strand == 'f':
                try:
                    alignment = aligner.align(ref_seq.seq[last_index:ref_pos],
                                              forward_seq[int(seq_nr)].seq[seq_l:seq_pos])[0]
                    list_of_alignments.append(
                        [100 * alignment.score / (ref_pos - last_index), alignment, last_index, ref_pos])
                except (IndexError, ValueError):
                    print("")
            else:
                try:
                    alignment = aligner.align(ref_seq.seq[last_index:ref_pos],
                                              reversed_seq[int(seq_nr)].seq[seq_l:seq_pos])[0]
                    list_of_alignments.append(
                        [100 * alignment.score / (ref_pos - last_index), alignment, last_index, ref_pos])
                except (IndexError, ValueError):
                    print("")
            last_index = ref_pos
    list_of_alignments.sort(key=lambda x: x[0], reverse=True)
    return list_of_alignments
