def align(seq1, seq2, match, mismatch, gap, mode):
    n, m = len(seq1), len(seq2)

    dp = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[None] * (m + 1) for _ in range(n + 1)]

    # --- Инициализация (только для global) ---
    if mode == "global":
        for i in range(1, n + 1):
            dp[i][0] = i * gap
            trace[i][0] = "up"
        for j in range(1, m + 1):
            dp[0][j] = j * gap
            trace[0][j] = "left"

    # --- Заполнение матрицы ---
    best_i = best_j = 0
    best_score = float("-inf")

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            if seq1[i-1] == seq2[j-1]:
                diag = dp[i-1][j-1] + match
            else:
                diag = dp[i-1][j-1] + mismatch

            up = dp[i-1][j] + gap
            left = dp[i][j-1] + gap

            if mode == "local":
                dp[i][j] = max(0, diag, up, left)
            else:
                dp[i][j] = max(diag, up, left)

            # направление
            if dp[i][j] == diag:
                trace[i][j] = "diag"
            elif dp[i][j] == up:
                trace[i][j] = "up"
            elif dp[i][j] == left:
                trace[i][j] = "left"
            else:
                trace[i][j] = None

            if mode == "local" and dp[i][j] > best_score:
                best_score = dp[i][j]
                best_i, best_j = i, j

    # --- старт traceback ---
    if mode == "global":
        i, j = n, m
        best_score = dp[n][m]
    else:
        i, j = best_i, best_j

    a1, a2 = [], []
    matches = mismatches = gaps = 0

    # --- восстановление ---
    while i > 0 or j > 0:
        if mode == "local" and dp[i][j] == 0:
            break

        move = trace[i][j]

        if move == "diag":
            a1.append(seq1[i-1])
            a2.append(seq2[j-1])
            if seq1[i-1] == seq2[j-1]:
                matches += 1
            else:
                mismatches += 1
            i -= 1
            j -= 1

        elif move == "up":
            a1.append(seq1[i-1])
            a2.append("-")
            gaps += 1
            i -= 1

        elif move == "left":
            a1.append("-")
            a2.append(seq2[j-1])
            gaps += 1
            j -= 1

        else:
            break

    a1.reverse()
    a2.reverse()

    return "".join(a1), "".join(a2), matches, mismatches, gaps, best_score


def read_sequence(prompt):
    """Читает многострочную последовательность до пустой строки"""
    print(prompt)
    lines = []
    while True:
        line = input().strip().upper()
        if not line:  # пустая строка - конец ввода
            break
        lines.append(line)
    return "".join(lines)

# ===== ВВОД =====
seq1 = read_sequence("Введите seq1 (для окончания введите пустую строку):")
seq2 = read_sequence("Введите seq2 (для окончания введите пустую строку):")
print("Введите seq1:")

match = int(input("match: "))
mismatch = int(input("mismatch: "))
gap = int(input("gap: "))

mode = input("global/local (G/L): ").lower()
mode = "global" if mode in ["g", "global"] else "local"

# ===== ЗАПУСК =====
a1, a2, matches, mismatches, gaps, score = align(
    seq1, seq2, match, mismatch, gap, mode
)

# ===== ВЫВОД =====
print("\nВыравнивание:")
print(a1)
print(a2)

print("\nСтатистика:")
print("Совпадений:", matches)
print("Несовпадений:", mismatches)
print("Гэпов:", gaps)

print("\nScore:", score)