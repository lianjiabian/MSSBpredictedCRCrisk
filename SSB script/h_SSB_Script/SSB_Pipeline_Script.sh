## SSB Pipeline Script

# 1 Filter Reads by 1_filter_paired_flags.py

cat SampleList.txt | while read id;do time python ../1_filter_paired_flags.py ${id}_1.fq ${id}_2.fq GGGGGGGGGG TTTTTTTTTTTT;done

# 2 Align gencode GRCh37.p13 genome by bwa
cat SampleList.txt | tail -n +2 | while read id;do bwa mem -t 40 -M -R '@RG\tID:foo\tSM:bar\tLB:library1' /data/Index/Gencode/GRCh37.p13/GRCh37.p13 1_Filter/${id}_1.fq.trim 1_Filter/${id}_2.fq.trim > 2_AlignSam/${id}.sam;don