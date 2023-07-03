

def generate_kmers(sequence, k):
    """
    Given a sequence of characters and a length k, generate all possible k-mers
    (substrings) of length k from the sequence.
    """
    kmers = []
    n = len(sequence)
    for i in range(n - k + 1):
        kmer = sequence[i:i+k]
        kmers.append(kmer)
    return kmers



def union_vectors(vector1, vector2):
    """
    Given two vectors, return a new vector that is the union of the two input vectors.
    """
    union = list(set(vector1).union(set(vector2)))
    return union


def count_items(vector, sequence):
    """
    Given a vector and a sequence, count how many times each item in the vector appears in the sequence.
    """
    counts = []
    for item in vector:
        counts.append(sequence.count(item))
    return counts

def count_items_norm(vector, sequence):
    """
    Given a vector and a sequence, count how many times each item in the vector appears in the sequence.
    """
    counts = []
    for item in vector:
        counts.append(sequence.count(item)/len(sequence))
    return counts