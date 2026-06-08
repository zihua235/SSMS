import bz2
import csv
infile = "/Users/zihuahuang/Documents/humann/metacyc_reactions_level4ec_only.uniref.bz2"

outfile = "metacyc_reactions_level4ec_only_uniref.csv"

with bz2.open(infile, "rt") as fin, open(outfile, "w", newline="") as fout:
    reader = csv.reader(fin, delimiter="\t")
    writer = csv.writer(fout)
    # 写表头
    writer.writerow(["Reaction", "EC", "UniRef_families..."])
    # 写每一行
    for row in reader:
        writer.writerow(row)

print(f"✅ 导出完成: {outfile}")