# Smith-Waterman Algorithm Implementation in Python

def read_input(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    seq1 = lines[0].strip()
    seq2 = lines[1].strip()
    gap_penalty = int(lines[2].strip())
    mismatch_penalty = int(lines[3].strip())
    match_score = int(lines[4].strip())
    return seq1, seq2, gap_penalty, mismatch_penalty, match_score

def initialize_matrix(rows, cols):
    return [[0 for _ in range(cols)] for _ in range(rows)]

def smith_waterman(seq1, seq2, gap_penalty, mismatch_penalty, match_score):
    m, n = len(seq1), len(seq2)
    H = initialize_matrix(m+1, n+1)
    max_score = 0
    max_pos = None

    # Scoring matrix calculation
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = H[i-1][j-1] + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty)
            delete = H[i-1][j] + gap_penalty
            insert = H[i][j-1] + gap_penalty
            H[i][j] = max(0, match, delete, insert)
            if H[i][j] >= max_score:
                max_score = H[i][j]
                max_pos = (i, j)

    return H, max_score, max_pos

def traceback(H, seq1, seq2, max_pos, gap_penalty, mismatch_penalty, match_score):
    aligned_seq1 = ''
    aligned_seq2 = ''
    i, j = max_pos

    while H[i][j] != 0:
        current_score = H[i][j]
        diagonal_score = H[i-1][j-1]
        up_score = H[i-1][j]
        left_score = H[i][j-1]

        if current_score == diagonal_score + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif current_score == up_score + gap_penalty:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        elif current_score == left_score + gap_penalty:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        else:
            break

    return aligned_seq1, aligned_seq2

def print_matrix(H, seq1, seq2):
    print("-----------------------------------------------------------** matrix **")
    print("===========================================================")
    seq2 = ' ' + seq2
    header = '     ' + '  '.join(seq2)
    print(header)
    for i, row in enumerate(H):
        if i == 0:
            print('  ', '  '.join(map(str, row)))
        else:
            print(seq1[i-1], '  '.join(map(str, row)))
    print("===========================================================")

def main():
    # Read input
    seq1, seq2, gap_penalty, mismatch_penalty, match_score = read_input('input.txt')

    # Compute Smith-Waterman algorithm
    H, max_score, max_pos = smith_waterman(seq1, seq2, gap_penalty, mismatch_penalty, match_score)

    # Traceback to get alignment
    aligned_seq1, aligned_seq2 = traceback(H, seq1, seq2, max_pos, gap_penalty, mismatch_penalty, match_score)

    # Print results
    print_matrix(H, seq1, seq2)
    print(f"Score = {max_score}")
    print(f"** Match = {match_score} | mismatch = {mismatch_penalty} | Gap = {gap_penalty} **")
    print("-----------------------------------------------------------")
    print("Alinhamento")
    print(' '.join(aligned_seq1))
    print(' '.join(aligned_seq2))
    print("-----------------------------------------------------------")

if __name__ == "__main__":
    main()
