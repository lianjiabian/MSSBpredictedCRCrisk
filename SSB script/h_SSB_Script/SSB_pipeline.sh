#!/bin/bash

# 基本目录设置


RAW_DATA_DIR=$1
#BASE_DIR=$(pwd)
BASE_DIR=$2
FILTERED_DIR="${BASE_DIR}/1_Filter"
ALIGN_DIR="${BASE_DIR}/2_Align_Sam"
FILTERED_BED_DIR="${BASE_DIR}/3_Filter_Bed"
BREAKS_BED_DIR="${BASE_DIR}/4_Breaks_Bed"
DEDUP_BED_DIR="${BASE_DIR}/5_Dedup_Bed"
DUPREPEAT_BED_DIR="${BASE_DIR}/6_Duprepeat_Bed"
ANNO_BED_DIR="${BASE_DIR}/7_Anno_Bed"
ANNO_BED_DIR_tss200="${BASE_DIR}/7_Tss200_Anno_Bed"
ANNO_BED_DIR_tss500="${BASE_DIR}/7_Tss500_Anno_Bed"
ANNO_BED_DIR_tss1000="${BASE_DIR}/7_Tss1000_Anno_Bed"
ANNO_BED_DIR_tss2000="${BASE_DIR}/7_Tss2000_Anno_Bed"
ANNO_BED_DIR_exon="${BASE_DIR}/7_Exon_Anno_Bed"
GENCODE_DIR="/data/Index/Gencode/GRCh37.p13"
RMSK_BED_FILE="/data/Index/Gencode/GRCh37.p13/hg19.ucsc.repeatmasker.bed"

get_filename() {
    id=$1
    suffix=$2 # "_1" 或 "_2"
    dir=$3    # 目录路径
    # 尝试找到匹配的文件名
    if [[ -f "${dir}/${id}${suffix}.clean.fq" ]]; then
        echo "${id}${suffix}.clean.fq"
    elif [[ -f "${dir}/${id}${suffix}.fq" ]]; then
        echo "${id}${suffix}.fq"
    elif [[ -f "${dir}/${id}${suffix}.clean.fq.trim" ]]; then
        echo "${id}${suffix}.clean.fq.trim"
    elif [[ -f "${dir}/${id}${suffix}.fq.trim" ]]; then
        echo "${id}${suffix}.fq.trim"
    else
        echo "Error: File not found for ${id}${suffix}" >&2
        exit 1
    fi
}

# 创建目录 ${RAW_DATA_DIR} 
mkdir -p ${FILTERED_DIR} ${ALIGN_DIR} ${FILTERED_BED_DIR} ${BREAKS_BED_DIR} ${DEDUP_BED_DIR} ${DUPREPEAT_BED_DIR} ${ANNO_BED_DIR} ${ANNO_BED_DIR_tss200} ${ANNO_BED_DIR_exon} ${ANNO_BED_DIR_tss500} ${ANNO_BED_DIR_tss1000} ${ANNO_BED_DIR_tss2000}

# 解压缩原始数据
gunzip ${RAW_DATA_DIR}/*.fq.gz

# 生成样本列表，支持两种命名格式
ls -1 ${RAW_DATA_DIR}/*_1.clean.fq ${RAW_DATA_DIR}/*_1.fq 2>/dev/null | xargs -n 1 basename | sed -e 's/_1.clean.fq//g' -e 's/_1.fq//g' > SampleList.txt

# 定义处理样本的函数
process_sample() {
  id=$1
  step=$2

  case $step in
    filter)
      file1=$(get_filename $id "_1" ${RAW_DATA_DIR})
      file2=$(get_filename $id "_2" ${RAW_DATA_DIR})
      python /data/Script/SSB_Script/1_filter_paired_flags.py ${RAW_DATA_DIR}/$file1 ${RAW_DATA_DIR}/$file2 GGGGGGGGGG TTTTTTTTTTTT
      echo "Step1 Filter is Over!"
      ;;
    align)
      file1_trim=$(get_filename $id "_1" ${FILTERED_DIR})
      file2_trim=$(get_filename $id "_2" ${FILTERED_DIR})
      bwa mem -t 40 -M -R '@RG\tID:foo\tSM:bar\tLB:library1' ${GENCODE_DIR}/GRCh37.p13 ${FILTERED_DIR}/$file1_trim ${FILTERED_DIR}/$file2_trim > ${ALIGN_DIR}/${id}.sam
      echo "Step2 align is Over!"
      ;;
    sam_to_bed)
      bash /data/Script/SSB_Script/3_sam_to_bed.sh ${ALIGN_DIR}/${id} ${FILTERED_BED_DIR}/${id}.candidates.bed;
      echo "Step3 sam_to_bed is Over!"
      ;;
    filter_polyA)
      python /data/Script/SSB_Script/5_filtering.internal.polyA.py ${FILTERED_BED_DIR}/${id}.candidates.bed ${GENCODE_DIR}/GRCh37.p13.genome.fa 40 0.4 > ${BREAKS_BED_DIR}/${id}.breaks.bed;
      echo "Step4 filter_polyA is Over!"
      ;;
    dedup)
      sort -k1,1 -k2,2n -k3,3n -k6,6 ${BREAKS_BED_DIR}/${id}.breaks.bed | awk '!seen[$1,$2,$3,$6]++' > ${DEDUP_BED_DIR}/${id}.deduplicated_output.bed;
      echo "Step5 dedup is Over!"
      ;;
    derepeat)
      bedtools intersect -v -a ${DEDUP_BED_DIR}/${id}.deduplicated_output.bed -b ${RMSK_BED_FILE} > ${DUPREPEAT_BED_DIR}/${id}.derepeat.bed;
      echo "Step6 derepeat is Over!"
      ;;
    annotate)
      bedtools intersect -a ${GENCODE_DIR}/gencode.v19.annotation.bed -b ${DUPREPEAT_BED_DIR}/${id}.derepeat.bed -wa -wb > ${ANNO_BED_DIR}/${id}.anno.bed;echo "Step7 annotate is Over!";
      bedtools intersect -a ${GENCODE_DIR}/gencode.v19.tss.200.bed -b ${DUPREPEAT_BED_DIR}/${id}.derepeat.bed -wa -wb > ${ANNO_BED_DIR_tss200}/${id}.anno.bed;
      bedtools intersect -a ${GENCODE_DIR}/gencode.v19.tss.500.bed -b ${DUPREPEAT_BED_DIR}/${id}.derepeat.bed -wa -wb > ${ANNO_BED_DIR_tss500}/${id}.anno.bed;
      bedtools intersect -a ${GENCODE_DIR}/gencode.v19.tss.1000.bed -b ${DUPREPEAT_BED_DIR}/${id}.derepeat.bed -wa -wb > ${ANNO_BED_DIR_tss1000}/${id}.anno.bed;
      bedtools intersect -a ${GENCODE_DIR}/gencode.v19.tss.2000.bed -b ${DUPREPEAT_BED_DIR}/${id}.derepeat.bed -wa -wb > ${ANNO_BED_DIR_tss2000}/${id}.anno.bed;
      bedtools intersect -a ${GENCODE_DIR}/gencode.v19.exon.bed -b ${DUPREPEAT_BED_DIR}/${id}.derepeat.bed -wa -wb > ${ANNO_BED_DIR_exon}/${id}.anno.bed;

      cut -f 4 ${ANNO_BED_DIR}/${id}.anno.bed | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2 -nr > ${ANNO_BED_DIR}/${id}.count.tsv;
      cut -f 4 ${ANNO_BED_DIR_tss200}/${id}.anno.bed | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2 -nr > ${ANNO_BED_DIR_tss200}/${id}.tss.count.tsv;
      cut -f 4 ${ANNO_BED_DIR_tss500}/${id}.anno.bed | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2 -nr > ${ANNO_BED_DIR_tss500}/${id}.tss.count.tsv;
      cut -f 4 ${ANNO_BED_DIR_tss1000}/${id}.anno.bed | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2 -nr > ${ANNO_BED_DIR_tss1000}/${id}.tss.count.tsv;
      cut -f 4 ${ANNO_BED_DIR_tss2000}/${id}.anno.bed | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2 -nr > ${ANNO_BED_DIR_tss2000}/${id}.tss.count.tsv;
      cut -f 5 ${ANNO_BED_DIR_exon}/${id}.anno.bed | sort | uniq -c | awk '{print $2"\t"$1}' | sort -k2 -nr > ${ANNO_BED_DIR_exon}/${id}.exon.count.tsv;
      echo "Step7 annotate is Over!"
      ;;
    *)
      echo "Unknown step: $step"
      ;;
  esac
}

# 处理样本
while read id; do
  process_sample "$id" filter
  mv ${RAW_DATA_DIR}/*trim ${FILTERED_DIR}
  process_sample "$id" align
  process_sample "$id" sam_to_bed
  process_sample "$id" filter_polyA
  process_sample "$id" dedup
  process_sample "$id" derepeat
  process_sample "$id" annotate
done < SampleList.txt

