def align(seq1, seq2, match, mismatch, gap_open, gap_extend, mode):
    n, m = len(seq1), len(seq2)
    inf = float('inf')

    # --- Инициализация 3-х матриц ---
    # M  - выравнивание заканчивается совпадением/несовпадением
    # Ix - выравнивание заканчивается гэпом в seq2 (символ seq1 напротив '-')
    # Iy - выравнивание заканчивается гэпом в seq1 ('-' напротив символа seq2)
    M = [[0] * (m + 1) for _ in range(n + 1)]
    Ix = [[0] * (m + 1) for _ in range(n + 1)]
    Iy = [[0] * (m + 1) for _ in range(n + 1)]

    # Матрицы для восстановления пути (traceback)
    ptrM = [[None] * (m + 1) for _ in range(n + 1)]
    ptrIx = [[None] * (m + 1) for _ in range(n + 1)]
    ptrIy = [[None] * (m + 1) for _ in range(n + 1)]

    # --- Базовые случаи ---
    if mode == "global":
        M[0][0] = 0
        Ix[0][0] = -inf
        Iy[0][0] = -inf

        for i in range(1, n + 1):
            M[i][0] = -inf
            Ix[i][0] = gap_open + i * gap_extend
            Iy[i][0] = -inf
            ptrIx[i][0] = "Ix"

        for j in range(1, m + 1):
            M[0][j] = -inf
            Ix[0][j] = -inf
            Iy[0][j] = gap_open + j * gap_extend
            ptrIy[0][j] = "Iy"

    # Для local базовая инициализация нулями подходит, оставляем как есть

    best_score_local = 0
    best_i_local, best_j_local = 0, 0
    best_state_local = "M"

    # --- Заполнение матриц ---
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # 1. Считаем Ix (движение вниз, гэп в seq2)
            val_M_Ix = M[i - 1][j] + gap_open + gap_extend if M[i - 1][j] != -inf else -inf
            val_Ix_Ix = Ix[i - 1][j] + gap_extend if Ix[i - 1][j] != -inf else -inf

            Ix[i][j] = max(val_M_Ix, val_Ix_Ix)
            ptrIx[i][j] = "M" if Ix[i][j] == val_M_Ix else "Ix"

            # 2. Считаем Iy (движение вправо, гэп в seq1)
            val_M_Iy = M[i][j - 1] + gap_open + gap_extend if M[i][j - 1] != -inf else -inf
            val_Iy_Iy = Iy[i][j - 1] + gap_extend if Iy[i][j - 1] != -inf else -inf

            Iy[i][j] = max(val_M_Iy, val_Iy_Iy)
            ptrIy[i][j] = "M" if Iy[i][j] == val_M_Iy else "Iy"

            # 3. Считаем M (движение по диагонали)
            score = match if seq1[i - 1] == seq2[j - 1] else mismatch

            val_M = M[i - 1][j - 1] + score if M[i - 1][j - 1] != -inf else -inf
            val_Ix = Ix[i - 1][j - 1] + score if Ix[i - 1][j - 1] != -inf else -inf
            val_Iy = Iy[i - 1][j - 1] + score if Iy[i - 1][j - 1] != -inf else -inf

            max_M = max(val_M, val_Ix, val_Iy)

            if mode == "local":
                max_M = max(0, max_M)

            M[i][j] = max_M

            if mode == "local" and M[i][j] == 0:
                ptrM[i][j] = None
            elif max_M == val_M:
                ptrM[i][j] = "M"
            elif max_M == val_Ix:
                ptrM[i][j] = "Ix"
            else:
                ptrM[i][j] = "Iy"

            # Отслеживаем лучший результат для local
            if mode == "local":
                current_max = max(M[i][j], Ix[i][j], Iy[i][j])
                if current_max > best_score_local:
                    best_score_local = current_max
                    best_i_local, best_j_local = i, j
                    if current_max == M[i][j]:
                        best_state_local = "M"
                    elif current_max == Ix[i][j]:
                        best_state_local = "Ix"
                    else:
                        best_state_local = "Iy"

    # --- Старт Traceback ---
    if mode == "global":
        i, j = n, m
        best_score = max(M[n][m], Ix[n][m], Iy[n][m])
        if best_score == M[n][m]:
            state = "M"
        elif best_score == Ix[n][m]:
            state = "Ix"
        else:
            state = "Iy"
    else:
        i, j = best_i_local, best_j_local
        best_score = best_score_local
        state = best_state_local

    a1, a2 = [], []

    # --- Восстановление пути ---
    while i > 0 or j > 0:
        if mode == "local":
            if state == "M" and M[i][j] == 0: break
            if state == "Ix" and Ix[i][j] <= 0: break
            if state == "Iy" and Iy[i][j] <= 0: break

        if state == "M":
            if i == 0 or j == 0: break
            a1.append(seq1[i - 1])
            a2.append(seq2[j - 1])
            state = ptrM[i][j]
            i -= 1
            j -= 1
        elif state == "Ix":
            if i == 0: break
            a1.append(seq1[i - 1])
            a2.append("-")
            state = ptrIx[i][j]
            i -= 1
        elif state == "Iy":
            if j == 0: break
            a1.append("-")
            a2.append(seq2[j - 1])
            state = ptrIy[i][j]
            j -= 1

    a1.reverse()
    a2.reverse()
    res1, res2 = "".join(a1), "".join(a2)

    # --- Подсчет статистики по готовому выравниванию ---
    matches = mismatches = gap_opens = total_gaps = 0
    in_gap1 = in_gap2 = False

    for k in range(len(res1)):
        c1, c2 = res1[k], res2[k]

        # Считаем совпадения/несовпадения и общее число прочерков
        if c1 != '-' and c2 != '-':
            if c1 == c2:
                matches += 1
            else:
                mismatches += 1
        else:
            total_gaps += 1

        # Считаем открытия гэпов в seq1
        if c1 == '-':
            if not in_gap1:
                gap_opens += 1
                in_gap1 = True
        else:
            in_gap1 = False

        # Считаем открытия гэпов в seq2
        if c2 == '-':
            if not in_gap2:
                gap_opens += 1
                in_gap2 = True
        else:
            in_gap2 = False

    return res1, res2, matches, mismatches, gap_opens, total_gaps, best_score


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
if __name__ == "__main__":
    seq1 = read_sequence("Введите seq1 (для окончания введите пустую строку):")
    seq2 = read_sequence("Введите seq2 (для окончания введите пустую строку):")

    match = int(input("Бонус за совпадение (положительное, например 5): "))
    mismatch = int(input("Штраф за несовпадение (отрицательное, например -4): "))
    gap_open = int(input("Штраф за открытие гэпа (отрицательное, например -9): "))
    gap_extend = int(input("Штраф за продление гэпа (отрицательное, например -1): "))

    mode = input("global/local (G/L): ").lower()
    mode = "global" if mode in ["g", "global"] else "local"

    # ===== ЗАПУСК =====
    a1, a2, matches, mismatches, gap_opens, total_gaps, score = align(
        seq1, seq2, match, mismatch, gap_open, gap_extend, mode
    )

    # ===== ВЫВОД =====
    print("\nОптимальное выравнивание:")
    print(a1)
    print(a2)

    print("\nСтатистика:")
    print(f"Совпадений: {matches}")
    print(f"Несовпадений: {mismatches}")
    print(f"Открытий гэпов: {gap_opens}")
    print(f"Общее число прочерков (длина гэпов): {total_gaps}")

    print(f"\nScore (Целевая функция): {score}")