from collections import defaultdict


def count_minimizers(nested_dict):
    minimizer_counts = defaultdict(int)

    for sequence_index in nested_dict:
        print("COUNTING:", sequence_index)
        for minimizer in nested_dict[sequence_index]:
            minimizer_counts[minimizer] += 1

    single_occurrence_minimizers = {k: v for k, v in minimizer_counts.items() if v <= 2}
    return single_occurrence_minimizers


def remove_recurring_keys(nested_dict):
    single_occurrence_minimizers = count_minimizers(nested_dict)

    for sequence_index in nested_dict:
        print("The current sequence is:", sequence_index)
        nested_dict[sequence_index] = {k: v for k, v in nested_dict[sequence_index].items() if
                                       k in single_occurrence_minimizers}
