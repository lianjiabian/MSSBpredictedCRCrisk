in=$1
out=$2

# 使用samtools和awk来过滤并立即转换为BAM格式
samtools view -h -q 20 -F 4 $in.sam | \
awk '($1 ~ /^@/) ($2 == 147) || ($2 == 163) {print}' | \
samtools view -b - > $in.temp.bam

# 使用bedtools将BAM转换为BED
bedtools bamtobed -tag NM -i $in.temp.bam > $in.temp.bed

# 根据链方向调整坐标，并排序
awk '{
    if ($6 == "+") {
        print $1"\t"$2"\t"$2+1"\t"$4"\t"$5"\t"$6
    } else if ($6 == "-") {
        print $1"\t"$3-1"\t"$3"\t"$4"\t"$5"\t"$6
    }
}' $in.temp.bed | sort -k1,1 -k2,2n > $out

# 清除临时文件
rm -f $in.temp.bam $in.temp.bed