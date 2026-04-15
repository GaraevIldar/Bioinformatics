import tkinter as tk
from tkinter import filedialog

def read_fasta(file_path):
    lengths = []
    seq = ""

    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq:
                    lengths.append(len(seq))
                    seq = ""
            else:
                seq += line
        if seq:
            lengths.append(len(seq))

    return lengths


def calculate_n50(lengths):
    lengths_sorted = sorted(lengths, reverse=True)
    total_length = sum(lengths_sorted)
    half = total_length / 2

    current = 0
    for l in lengths_sorted:
        current += l
        if current >= half:
            return l


def main():
    root = tk.Tk()
    root.withdraw()

    file_path = filedialog.askopenfilename(
        title="Выберите FASTA файл",
        filetypes=[("FASTA files", "*.fasta *.fa *.fna"), ("All files", "*.*")]
    )

    if not file_path:
        print("Файл не выбран")
        return

    lengths = read_fasta(file_path)

    count = len(lengths)
    min_len = min(lengths)
    max_len = max(lengths)
    avg_len = sum(lengths) / count
    n50 = calculate_n50(lengths)

    print("Количество последовательностей:", count)
    print("Минимальная длина:", min_len)
    print("Максимальная длина:", max_len)
    print("Средняя длина:", avg_len)
    print("N50:", n50)


main()