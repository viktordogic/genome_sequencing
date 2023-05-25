def create_cigar_string(referent, sequence):
    cigar_string = ''
    i = 0
    j = 0
    counter = 0

    while i < len(referent) and j < len(sequence):
        if referent[i] == sequence[j]:
            cigar_string += 'M'
            i += 1
            j += 1
            counter += 1
        elif i < len(referent) - 1 and referent[i+1] == sequence[j]:
            cigar_string += 'D'
            i += 2
            j += 1
        elif j < len(sequence) - 1 and sequence[j+1] == referent[i]:
            cigar_string += 'I'
            i += 1
            j += 2
        else:
            cigar_string += 'X'
            i += 1
            j += 1

    # In case sequence or referent still has unprocessed elements
    while i < len(referent):
        cigar_string += 'D'
        i += 1
    while j < len(sequence):
        cigar_string += 'I'
        j += 1

    return cigar_string
