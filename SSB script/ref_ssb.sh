#下载gtf fa repeatmasker


samtools faidx GRCh37.p13.genome.fa

cat gencode.v19.annotation.gtf | awk -F'\t' '$3=="gene" {
    split($9,a,"; ");
    for (i in a) {
        if (a[i] ~ /gene_id/) {id=a[i]} 
        if (a[i] ~ /gene_name/) {name=a[i]} 
        if (a[i] ~ /gene_type/) {type=a[i]}
    }
    sub(/gene_id "/, "", id); 
    sub(/gene_name "/, "", name); 
    sub(/gene_type "/, "", type); 
    gsub(/"/, "", id); 
    gsub(/"/, "", name); 
    gsub(/"/, "", type);
    print $1"\t"($4)"\t"$5"\t"name"\t.""\t"$7"\t"type"\t"id
}' > gencode.v19.annotation.bed


cat gencode.v19.annotation.gtf | awk -F'\t' '$3=="gene" {
    split($9,a,"; ");
    for (i in a) {
        if (a[i] ~ /gene_id/) {id=a[i]} 
        if (a[i] ~ /gene_name/) {name=a[i]} 
        if (a[i] ~ /gene_type/) {type=a[i]}
    }
    sub(/gene_id "/, "", id); 
    sub(/gene_name "/, "", name); 
    sub(/gene_type "/, "", type); 
    gsub(/"/, "", id); 
    gsub(/"/, "", name); 
    gsub(/"/, "", type);
    if ($7 == "+") {
        start = $4-1;  # 对于正链，起始位置是TSS
        end = start + 1;
    } else {
        end = $5;      # 对于负链，结束位置是TSS
        start = end - 1;
    }
    print $1"\t"start"\t"end"\t"name"\t.\t"$7"\t"type"\t"id
}' >  gencode.v19.tss.bed

bedtools slop -i gencode.v19.tss.bed -g GRCh37.p13.genome.fa.fai -b 200 > gencode.v19.tss.200.bed
bedtools slop -i gencode.v19.tss.bed -g GRCh37.p13.genome.fa.fai -b 500 > gencode.v19.tss.500.bed
bedtools slop -i gencode.v19.tss.bed -g GRCh37.p13.genome.fa.fai -b 1000 > gencode.v19.tss.1000.bed
bedtools slop -i gencode.v19.tss.bed -g GRCh37.p13.genome.fa.fai -b 2000 > gencode.v19.tss.2000.bed


awk -F'\t' '$3 == "exon" {
    split($9,a,"; ");
    for (i in a) {
        if (a[i] ~ /gene_id/) {id=a[i]}
        if (a[i] ~ /gene_name/) {name=a[i]}
    }
    sub(/gene_id "/, "", id); 
    sub(/gene_name "/, "", name); 
    gsub(/"/, "", id); 
    gsub(/"/, "", name);
    print $1"\t"$4-1"\t"$5"\t"id"\t"name"\t"$7
}' gencode.v19.annotation.gtf > gencode.v19.exon.bed


bwa index GRCh37.p13.genome.fa -p GRCh37.p13