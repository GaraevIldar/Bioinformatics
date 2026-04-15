import random

# Словарь для транскрипции ДНК в РНК и трансляции в белок
DNA_TO_RNA = {'A': 'A', 'T': 'U', 'G': 'G', 'C': 'C'}
RNA_TO_PROTEIN = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}
START_CODON = 'AUG'
STOP_CODONS = {'UAA', 'UAG', 'UGA'}


def get_complement(dna):
    """Возвращает комплементарную последовательность ДНК"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(complement[base] for base in dna)


def reverse_complement(dna):
    """Возвращает обратно-комплементарную последовательность"""
    return get_complement(dna)[::-1]


def dna_to_rna(dna):
    """Транскрипция ДНК в РНК"""
    return ''.join(DNA_TO_RNA[base] for base in dna)


def translate_rna(rna):
    """Трансляция РНК в белок"""
    protein = []
    for i in range(0, len(rna) - 2, 3):
        codon = rna[i:i + 3]
        if len(codon) == 3:
            aa = RNA_TO_PROTEIN.get(codon, 'X')
            protein.append(aa)
            if aa == '*':  # стоп-кодон
                break
    return ''.join(protein)


def find_orfs_in_sequence(seq, frame, chain_type):
    """Поиск всех ORF в одной последовательности для одной рамки"""
    orfs = []
    rna_seq = dna_to_rna(seq)

    i = frame
    while i < len(rna_seq) - 2:
        codon = rna_seq[i:i + 3]
        if codon == START_CODON:
            # Нашли старт-кодон, ищем стоп-кодон
            for j in range(i, len(rna_seq) - 2, 3):
                stop_codon = rna_seq[j:j + 3]
                if stop_codon in STOP_CODONS:
                    orf_length = j - i + 3
                    # Проверяем минимальную длину белка (10 аминокислот)
                    if orf_length >= 33:  # 10 аминокидайслот * 3 + старт-кодон
                        orfs.append({
                            'start': i,
                            'end': j + 3,
                            'seq': rna_seq[i:j + 3],
                            'frame': frame + 1,
                            'chain': chain_type
                        })
                    break
            i = j + 3
        else:
            i += 3
    return orfs


def find_longest_orf(dna_sequence):
    """Поиск самой длинной ORF во всех 6 рамках"""
    all_orfs = []

    # Прямая цепь
    for frame in range(3):
        orfs = find_orfs_in_sequence(dna_sequence, frame, "прямая")
        all_orfs.extend(orfs)

    # Обратная цепь
    rev_comp = reverse_complement(dna_sequence)
    for frame in range(3):
        orfs = find_orfs_in_sequence(rev_comp, frame, "обратная")
        # Корректируем координаты для обратной цепи
        for orf in orfs:
            orf['start'] = len(dna_sequence) - orf['end']
            orf['end'] = len(dna_sequence) - orf['start']
        all_orfs.extend(orfs)

    # Находим самую длинную ORF
    if not all_orfs:
        return None

    longest_orf = max(all_orfs, key=lambda x: len(x['seq']))
    return longest_orf


def format_output(orf, dna_sequence):
    """Форматирование вывода"""
    if not orf:
        print("ORF не найдена")
        return

    print(f"\nНайдена ORF на {orf['chain']} цепи в рамке {orf['frame']}")
    print(f"Диапазон: {orf['start'] + 1} - {orf['end']}")

    # Разделение на триплеты
    triplets = ' '.join([orf['seq'][i:i + 3] for i in range(0, len(orf['seq']), 3)])
    print(f"ORF (триплеты): {triplets}")

    # Трансляция в белок
    protein = translate_rna(orf['seq'])
    print(f"Белок (однобуквенный код): {protein}")


def generate_dna(length, gc_percent):
    """Генерация случайной ДНК с заданным GC-составом"""
    # Расчет количества нуклеотидов
    gc_count = round(length * gc_percent / 100)
    at_count = length - gc_count

    g_count = gc_count // 2
    c_count = gc_count - g_count
    a_count = at_count // 2
    t_count = at_count - a_count

    # Создание пула нуклеотидов
    nucleotides = ['G'] * g_count + ['C'] * c_count + ['A'] * a_count + ['T'] * t_count

    # Перемешивание
    random.shuffle(nucleotides)

    return ''.join(nucleotides)


def main():
    # Ввод данных с проверкой
    while True:
        try:
            length = int(input("Введите длину ДНК (100-1000): "))
            if 100 <= length <= 1000:
                break
            else:
                print("Ошибка: длина должна быть от 100 до 1000")
        except ValueError:
            print("Ошибка: введите целое число")

    while True:
        try:
            gc = int(input("Введите GC-состав в процентах (20-80): "))
            if 20 <= gc <= 80:
                break
            else:
                print("Ошибка: GC-состав должен быть от 20 до 80")
        except ValueError:
            print("Ошибка: введите целое число")

    # Генерация ДНК
    dna = generate_dna(length, gc)
    print(f"\nСгенерированная ДНК: {dna}")

    # Подсчет реального состава
    g_count = dna.count('G')
    c_count = dna.count('C')
    a_count = dna.count('A')
    t_count = dna.count('T')

    print(f"\nРеальный состав:")
    print(f"G: {g_count} ({g_count / length * 100:.1f}%)")
    print(f"C: {c_count} ({c_count / length * 100:.1f}%)")
    print(f"A: {a_count} ({a_count / length * 100:.1f}%)")
    print(f"T: {t_count} ({t_count / length * 100:.1f}%)")

    # Поиск ORF
    longest_orf = find_longest_orf(dna)
    format_output(longest_orf, dna)


main()
