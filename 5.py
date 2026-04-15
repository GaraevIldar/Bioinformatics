import subprocess
import os
import sys
from Bio import SeqIO

# ==============================
# 🔧 ПУТИ К BLAST (ТВОИ)
# ==============================
MAKEBLASTDB = r"C:\Program Files\NCBI\blast-2.17.0+\bin\makeblastdb.exe"
BLASTN = r"C:\Program Files\NCBI\blast-2.17.0+\bin\blastn.exe"

# ==============================
# 📁 ВХОДНЫЕ ФАЙЛЫ
# ==============================
file_188 = "188.ffn"
file_190 = "190.ffn"

# ==============================
# 📁 ВЫХОДНЫЕ ФАЙЛЫ
# ==============================
blast_188_vs_190 = "188_vs_190.tsv"
blast_190_vs_188 = "190_vs_188.tsv"
output_file = "orthologs.tsv"

# ==============================
# ✅ ПРОВЕРКИ
# ==============================
def check_files():
    if not os.path.exists(MAKEBLASTDB):
        print("❌ Не найден makeblastdb:")
        print(MAKEBLASTDB)
        sys.exit()

    if not os.path.exists(BLASTN):
        print("❌ Не найден blastn:")
        print(BLASTN)
        sys.exit()

    if not os.path.exists(file_188):
        print(f"❌ Нет файла: {file_188}")
        sys.exit()

    if not os.path.exists(file_190):
        print(f"❌ Нет файла: {file_190}")
        sys.exit()

    print("✅ Все файлы найдены\n")

# ==============================
# 🧱 СОЗДАНИЕ БАЗЫ
# ==============================
def make_db(fasta, db_name):
    subprocess.run([
        MAKEBLASTDB,
        "-in", fasta,
        "-dbtype", "nucl",
        "-out", db_name
    ], check=True)

# ==============================
# 🔍 BLAST
# ==============================
def run_blast(query, db, out):
    subprocess.run([
        BLASTN,
        "-query", query,
        "-db", db,
        "-out", out,
        "-outfmt", "6 qseqid sseqid bitscore",
        "-max_target_seqs", "1"
    ], check=True)

# ==============================
# 📊 ПАРСИНГ
# ==============================
def parse_best_hits(file):
    best_hits = {}
    with open(file) as f:
        for line in f:
            q, s, score = line.strip().split()
            best_hits[q] = s
    return best_hits

# ==============================
# 📖 ОПИСАНИЯ
# ==============================
def get_descriptions(fasta):
    desc = {}
    for record in SeqIO.parse(fasta, "fasta"):
        desc[record.id] = record.description
    return desc

# ==============================
# 🚀 ОСНОВНАЯ ПРОГРАММА
# ==============================
def main():
    check_files()

    print("🧱 Создание BLAST баз...")
    make_db(file_188, "db188")
    make_db(file_190, "db190")

    print("🔍 Запуск BLAST...")
    run_blast(file_188, "db190", blast_188_vs_190)
    run_blast(file_190, "db188", blast_190_vs_188)

    print("📊 Чтение результатов...")
    best_188_to_190 = parse_best_hits(blast_188_vs_190)
    best_190_to_188 = parse_best_hits(blast_190_vs_188)

    print("🧬 Поиск ортологов (RBH)...")
    orthologs = []

    for gene_188, gene_190 in best_188_to_190.items():
        if gene_190 in best_190_to_188:
            if best_190_to_188[gene_190] == gene_188:
                orthologs.append((gene_188, gene_190))

    print("📖 Чтение описаний...")
    desc_188 = get_descriptions(file_188)
    desc_190 = get_descriptions(file_190)

    print("💾 Сохранение результата...")
    with open(output_file, "w", encoding="utf-8") as out:
        out.write("ID_188\tID_190\tdescription_188\tdescription_190\n")
        for g188, g190 in orthologs:
            out.write(
                f"{g188}\t{g190}\t{desc_188.get(g188,'')}\t{desc_190.get(g190,'')}\n"
            )

    # ==============================
    # 📈 СТАТИСТИКА
    # ==============================
    genes_188 = sum(1 for _ in SeqIO.parse(file_188, "fasta"))
    genes_190 = sum(1 for _ in SeqIO.parse(file_190, "fasta"))

    print("\n=== РЕЗУЛЬТАТЫ ===")
    print("🧬 Генов в 188:", genes_188)
    print("🧬 Генов в 190:", genes_190)
    print("🔗 Ортологов 1-к-1:", len(orthologs))
    print(f"📄 Файл сохранён: {output_file}")

# ==============================
# ▶️ ЗАПУСК
# ==============================
if __name__ == "__main__":
    main()