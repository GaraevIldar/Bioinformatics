import gzip

bed_file = "dm6_refseq.bed"

# ---------- 1. Подсчет транскриптов ----------
protein_coding = set()
noncoding = set()

with open(bed_file, "r") as f:
    for line in f:
        cols = line.strip().split("\t")
        name = cols[3]

        if name.startswith("NM_"):
            protein_coding.add(name)
        elif name.startswith("NR_"):
            noncoding.add(name)

print("Protein-coding transcripts:", len(protein_coding))
print("Non-coding transcripts:", len(noncoding))


# ---------- 2. Преобразование BED12 → BED6 (экзоны) ----------
bed6 = []

with open(bed_file, "r") as f:
    for line in f:
        cols = line.strip().split("\t")

        chrom = cols[0]
        chromStart = int(cols[1])
        name = cols[3]
        score = cols[4]
        strand = cols[5]

        blockCount = int(cols[9])
        blockSizes = list(map(int, cols[10].strip(",").split(",")))
        blockStarts = list(map(int, cols[11].strip(",").split(",")))

        # разбиваем на экзоны
        for i in range(blockCount):
            exon_start = chromStart + blockStarts[i]
            exon_end = exon_start + blockSizes[i]

            bed6.append((chrom, exon_start, exon_end, name, score, strand))

# подсчет экзонов по цепям
plus_count = sum(1 for x in bed6 if x[5] == "+")
minus_count = sum(1 for x in bed6 if x[5] == "-")

print("Exons on + strand:", plus_count)
print("Exons on - strand:", minus_count)


# ---------- 3. Разделение ----------
plus_intervals = [x for x in bed6 if x[5] == "+"]
minus_intervals = [x for x in bed6 if x[5] == "-"]


# ---------- пересечение (+) и (-) ----------
from collections import defaultdict

def intersect_fast(a, b):
    result = []

    # группируем по хромосомам
    a_dict = defaultdict(list)
    b_dict = defaultdict(list)

    for chrom, start, end, *_ in a:
        a_dict[chrom].append((start, end))

    for chrom, start, end, *_ in b:
        b_dict[chrom].append((start, end))

    for chrom in a_dict:
        if chrom not in b_dict:
            continue

        a_list = sorted(a_dict[chrom])
        b_list = sorted(b_dict[chrom])

        i = j = 0

        while i < len(a_list) and j < len(b_list):
            a_start, a_end = a_list[i]
            b_start, b_end = b_list[j]

            start = max(a_start, b_start)
            end = min(a_end, b_end)

            if start < end:
                result.append((chrom, start, end))

            # двигаем указатели
            if a_end < b_end:
                i += 1
            else:
                j += 1

    return result


intersected = intersect_fast(plus_intervals, minus_intervals)


# ---------- сортировка ----------
sorted_intervals = sorted(intersected, key=lambda x: (x[0], x[1], x[2]))


# ---------- merge ----------
def merge_intervals(intervals):
    merged = []

    for chrom, start, end in intervals:
        if not merged:
            merged.append([chrom, start, end])
        else:
            last = merged[-1]

            if last[0] == chrom and start <= last[2]:
                last[2] = max(last[2], end)
            else:
                merged.append([chrom, start, end])

    return merged


merged = merge_intervals(sorted_intervals)


# ---------- подсчет нуклеотидов ----------
total_nt = sum(end - start for _, start, end in merged)

print("Nucleotides overlapping on both strands:", total_nt)
