import sys

def get_bed_file(file):
    reads = {}
    with open(file, 'r') as f:
        for line in f:
            chr, s, e, rid, _, strand = line.strip().split("\t")
            if chr not in reads:
                reads[chr] = []
            reads[chr].append((int(s), int(e), rid, strand))
    return reads

def get_genome(file):
    genome = {}
    with open(file, 'r') as f:
        key = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                key = line[1:].split()[0]
                genome[key] = []
            else:
                genome[key].append(line)
        for key in genome:
            genome[key] = "".join(genome[key])
    return genome

def main():
    if len(sys.argv) < 3:
        print("Usage: python script_name.py input_bed genome_fasta [extract_length] [ratio_cutoff]")
        sys.exit(1)
    
    reads_file = sys.argv[1]
    genome_file = sys.argv[2]
    extract_length = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    ratio_cutoff = float(sys.argv[4]) if len(sys.argv) > 4 else 0.4

    reads = get_bed_file(reads_file)
    genome = get_genome(genome_file)

    for chr in sorted(reads.keys()):
        chr_seq = genome.get(chr, "")
        for s, e, rid, strand in reads[chr]:
            if e < s:
                s, e = e, s
            rid_seq_for = chr_seq[s-extract_length:s].upper()
            rid_seq_back = chr_seq[e+1:e+1+extract_length].upper()[::-1].translate(str.maketrans('ACGT', 'TGCA'))
            
            AA_num = 0
            if strand == '+':
                AA_num = rid_seq_for.count('T')
            elif strand == '-':
                AA_num = rid_seq_back.count('T')
            
            aa_ratio = AA_num / extract_length
            
            if aa_ratio < ratio_cutoff:
                print(f"{chr}\t{s}\t{e}\t{rid}\t{aa_ratio:.2f}\t{strand}")

if __name__ == "__main__":
    main()
