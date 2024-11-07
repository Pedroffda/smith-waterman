def read_input(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    seq1 = lines[0].strip()
    seq2 = lines[1].strip()
    gap_penalty = int(lines[2].strip())
    mismatch_penalty = int(lines[3].strip())
    match_score = int(lines[4].strip())
    return seq1, seq2, gap_penalty, mismatch_penalty, match_score

def initialize_matrix(rows, cols, gap_penalty):
    H = [[0 for _ in range(cols)] for _ in range(rows)]
    for i in range(1, rows):
        H[i][0] = H[i-1][0] + gap_penalty
    for j in range(1, cols):
        H[0][j] = H[0][j-1] + gap_penalty
    return H

def smith_waterman_global(seq1, seq2, gap_penalty, mismatch_penalty, match_score):
    m, n = len(seq1), len(seq2)
    H = initialize_matrix(m+1, n+1, gap_penalty)

    for i in range(1, m+1):
        for j in range(1, n+1):
            if seq1[i-1] == seq2[j-1]:
                diag = H[i-1][j-1] + match_score
            else:
                diag = H[i-1][j-1] + mismatch_penalty
            delete = H[i-1][j] + gap_penalty
            insert = H[i][j-1] + gap_penalty
            H[i][j] = max(diag, delete, insert)
    max_score = H[0][n]
    max_pos = (0, n)
    for i in range(1, m+1):
        if H[i][n] > max_score:
            max_score = H[i][n]
            max_pos = (i, n)

    return H, max_score, max_pos

def traceback(H, seq1, seq2, start_pos, gap_penalty, mismatch_penalty, match_score):
    aligned_seq1 = ''
    aligned_seq2 = ''
    i, j = start_pos

    while i > 0 and j > 0:
        score_current = H[i][j]
        score_diag = H[i-1][j-1]
        score_left = H[i][j-1]

        if score_current == score_diag + (match_score if seq1[i-1] == seq2[j-1] else mismatch_penalty):
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif score_current == score_left + gap_penalty:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        else:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1

    while i > 0:
        aligned_seq1 = seq1[i-1] + aligned_seq1
        aligned_seq2 = '-' + aligned_seq2
        i -= 1
    while j > 0:
        aligned_seq1 = '-' + aligned_seq1
        aligned_seq2 = seq2[j-1] + aligned_seq2
        j -= 1

    return aligned_seq1, aligned_seq2

def print_matrix(H, seq1, seq2, output_file):
    with open(output_file, 'a') as f:
        f.write("-----------------------------------------------------------** matrix **\n")
        f.write("===========================================================\n")

        seq2 = ' ' + seq2
        header = '      ' + '    '.join(seq2)
        
        for i in range(len(H)-1, -1, -1):
            if i == 0:
                label = '-'
            else:
                label = seq1[i-1]
            row = f"{label:<2} " + ' '.join(f"{cell:4}" for cell in H[i])
            f.write(row + "\n")

        f.write(header + "\n")
        f.write("===========================================================\n")

def main():
    output_file = 'output.txt'
    
    # Limpar o arquivo de saída
    with open(output_file, 'w') as f:
        pass

    # Ler entrada
    seq1, seq2, gap_penalty, mismatch_penalty, match_score = read_input('input.txt')

    # Executar o algoritmo Needleman-Wunsch
    H, max_score, max_pos = smith_waterman_global(seq1, seq2, gap_penalty, mismatch_penalty, match_score)

    # Realizar o traceback
    aligned_seq1, aligned_seq2 = traceback(H, seq1, seq2, max_pos, gap_penalty, mismatch_penalty, match_score)

    # Imprimir resultados
    print_matrix(H, seq1, seq2, output_file)
    with open(output_file, 'a') as f:
        f.write(f"Score = {max_score}\n")
        f.write(f"** Match = {match_score} | mismatch = {mismatch_penalty} | Gap = {gap_penalty} **\n")
        f.write("-----------------------------------------------------------\n")
        f.write("Alinhamento\n")
        f.write(' '.join(aligned_seq1) + "\n")
        f.write(' '.join(aligned_seq2) + "\n")
        f.write("-----------------------------------------------------------\n")

if __name__ == "__main__":
    main()
