import tkinter as tk
from tkinter import filedialog
import os


def analyze_fastq(filename):
    try:
        with open(filename, 'r') as file:
            read_count = 0
            total_bases = 0
            total_quality_bases = 0

            line_num = 0
            for line in file:
                line = line.strip()
                line_num += 1

                if line_num % 4 == 2:
                    sequence = line
                    total_bases += len(sequence)
                    read_count += 1

                elif line_num % 4 == 0:
                    for char in line:
                        q_score = ord(char) - 33
                        if q_score >= 30:
                            total_quality_bases += 1

            if read_count > 0:
                avg_length = total_bases / read_count
                if total_bases > 0:
                    quality_percent = (total_quality_bases / total_bases) * 100
                else:
                    quality_percent = 0
            else:
                avg_length = 0
                quality_percent = 0

            print(f"РЕЗУЛЬТАТЫ АНАЛИЗА FASTQ ФАЙЛА")
            print(f"Количество ридов: {read_count}")
            print(f"Общее количество букв: {total_bases}")
            print(f"Средняя длина рида: {avg_length:.2f}")
            print(f"Процент букв с Q>=30: {quality_percent:.2f}%")

    except FileNotFoundError:
        print(f"Ошибка: файл '{filename}' не найден")
    except Exception as e:
        print(f"Ошибка при обработке файла: {e}")


def select_file_and_analyze():

    root = tk.Tk()
    root.withdraw()

    filename = filedialog.askopenfilename(
        title="Выберите FASTQ файл для анализа",
        filetypes=[("FASTQ files", "*.fastq *.fq"), ("All files", "*.*")]
    )

    if filename:
        analyze_fastq(filename)
    else:
        print("Файл не выбран")


if __name__ == "__main__":
    select_file_and_analyze()