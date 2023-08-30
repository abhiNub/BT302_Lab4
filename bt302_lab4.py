def pairwise_alignment(seq1, seq2, gap_penalty=-2, match_score=1, mismatch_penalty=-1):
    # Initialize the scoring matrix and traceback matrix
    n = len(seq1) + 1
    m = len(seq2) + 1
    dp = [[0] * m for _ in range(n)]

    for i in range(n):
        dp[i][0] = gap_penalty * i
    for j in range(m):
        dp[0][j] = gap_penalty * j

    for i in range(1, n):
        for j in range(1, m):
            match = dp[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = dp[i - 1][j] + gap_penalty
            insert = dp[i][j - 1] + gap_penalty
            dp[i][j] = max(match, delete, insert)

    aligned_seq1, aligned_seq2 = "", ""
    i, j = n - 1, m - 1

    while i > 0 and j > 0:
        if dp[i][j] == dp[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty):
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i - 1][j] + gap_penalty:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    while i > 0:
        aligned_seq1 = seq1[i - 1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        i -= 1
    while j > 0:
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[j - 1] + aligned_seq2
        j -= 1

    return aligned_seq1, aligned_seq2

def progressive_alignment(sequences, gap_penalty=-2, match_score=1, mismatch_penalty=-1):
    num_sequences = len(sequences)
    
    # Perform pair-wise alignment of all sequences
    alignments = [sequences[0]]
    for i in range(1, num_sequences):
        aligned_seq, _ = pairwise_alignment(alignments[i - 1], sequences[i], gap_penalty, match_score, mismatch_penalty)
        alignments.append(aligned_seq)
    
    # Create a profile alignment using the pair-wise alignments
    profile_alignment = alignments[0]
    for i in range(1, num_sequences):
        profile_alignment, _ = pairwise_alignment(profile_alignment, alignments[i], gap_penalty, match_score, mismatch_penalty)
    
    return profile_alignment

# Collect user input for protein sequences
num_sequences = int(input("Enter the number of protein sequences: "))
sequences = []

for i in range(num_sequences):
    sequence = input(f"Enter protein sequence {i + 1}: ")
    sequences.append(sequence)

# Perform progressive multiple sequence alignment
aligned_sequence = progressive_alignment(sequences)

# Print the aligned sequence
print("Progressive Multiple Sequence Alignment:")
print(aligned_sequence)