#Dataset download
##ground truth
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab//release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_CHM13v2.0_v5.0q_stvar.vcf.gz
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab//release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_CHM13v2.0_v5.0q_stvar.vcf.gz.tbi
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab//release/AshkenazimTrio/HG002_NA24385_son/v5.0q/HG002_CHM13v2.0_v5.0q_stvar.benchmark.bed

##HiFi
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/20kb/m64011_190830_220126.Q20.fastq
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/working/HPRC_PLUS/HG002/raw_data/PacBio_HiFi/20kb/m64011_190901_095311.Q20.fastq

##ONT
aria2c -x 16 -s 16 -c https://ont-open-data.s3.amazonaws.com/gm24385_2023.12/all_pass.vhg002v1.bam
samtools fastq all_pass.vhg002v1.bam > all_pass.vhg002v1.fastq

##reference genome
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz
gzip -d chm13v2.0.fa.gz
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/annotation/chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed

#mapping
##HIFI
minimap2 -a -t 64 -x map-hifi -Q --eqx --secondary=no -K4G --MD -Y -L chm13v2.0.fa.gz m64011_190830_220126.Q20.fastq m64011_190901_095311.Q20.fastq | samtools view -bh -@ 64 - | samtools sort -@ 64 -o HIFI.primary.sort.bam -
samtools index HIFI.primary.sort.bam
##ONT
minimap2 -a -t 64 -x map-ont -Q --eqx --secondary=no -K4G --MD -Y -L chm13v2.0.fa.gz all_pass.vhg002v1.fastq | samtools view -bh -@ 64 - | samtools sort -@ 64 -o ONT.primary.sort.bam -
samtools index ONT.primary.sort.bam

# Variant selection
# 50-999 bp: 20
bcftools view -H \
    -R HG002_CHM13v2.0_v5.0q_stvar.benchmark.bed \
    -i 'INFO/SVTYPE="INS" && INFO/SVLEN>50 && INFO/SVLEN<1000' \
    HG002_CHM13v2.0_v5.0q_stvar.vcf.gz \
    | shuf --random-source=<(yes 123456789) -n 20 \
    > insertion_GIAB_random50.vcf

# 1-9.9 kb: 20
bcftools view -H \
    -R HG002_CHM13v2.0_v5.0q_stvar.benchmark.bed \
    -i 'INFO/SVTYPE="INS" && INFO/SVLEN>=1000 && INFO/SVLEN<10000' \
    HG002_CHM13v2.0_v5.0q_stvar.vcf.gz \
    | shuf --random-source=<(yes 123456789) -n 20 \
    >> insertion_GIAB_random50.vcf

# >=10 kb: 10
bcftools view -H \
    -R HG002_CHM13v2.0_v5.0q_stvar.benchmark.bed \
    -i 'INFO/SVTYPE="INS" && INFO/SVLEN>=10000' \
    HG002_CHM13v2.0_v5.0q_stvar.vcf.gz \
    | shuf --random-source=<(yes 654987) -n 10 \
    >> insertion_GIAB_random50.vcf
    
awk '{a+=1}{print $1"\t"($2-1)"\t"($2-1)"\tINS"}' insertion_GIAB_random50.vcf > insertion_GIAB_random50.bed
awk '{a+=1}{print $1"\t"($2-1-150)"\t"($2-1+150)"\tINS"NR}' insertion_GIAB_random50.vcf > insertion_GIAB_random50.padding150.bed

#Visualization
time -v python bamsnap-lrs dna --bam HIFI.primary.sort.bam --bam ONT.primary.sort.bam --regions insertion_GIAB_random50.vcf --out-prefix bamsnap/insertion_ONT_HIFI.png --fa chm13v2.0.fa  --show-axis --show-coverage --padding 150  --detail high --overview-detail show > bamsnap_ins.runtime.log 2>&1 

time -v  hawkeye.py sv_browse -g chm13 -i HIFI.primary.sort.bam,ONT.primary.sort.bam  -b insertion_GIAB_random50.bed -r chm13v2.0.fa -o svhawkeye/ -q 0 -I 0 -d 150 -F png -f bed --sv_min_length 0  > svhawkeye_ins.runtime.log 2>&1 

time -v wally-v0.9.2-linux-amd64 region --map-qual 0 --snv-vaf 0 --snv-cov 0 --width 1200 -g chm13v2.0.fa -R insertion_GIAB_random50.padding150.bed HIFI.primary.sort.bam ONT.primary.sort.bam > wally_ins.runtime.log 2>&1

tar  -czvf  wally.INS.tar.gz  INS*png
tar  -czvf bamsnap_lrs.INS.tar.gz bamsnap/*png
tar  -czvf svhawkeye.INS.tar.gz svhawkeye/figure/*png

python annotate_repeat_context.py \
    -v insertion_GIAB_random50.vcf \
    -r chm13v2.0_RepeatMasker_4.1.2p1.2022Apr14.bed \
    --bedtools /share/app/bedtools/2.29.2/bin/bedtools
    
python bamsnap_benchmark_stats_per_bam.py dna \
    --bam HIFI.primary.sort.bam \
    --bam ONT.primary.sort.bam \
    --regions insertion_GIAB_random50.vcf \
    --stats-out ./insertion_ONT_HIFI.stats.tsv \
    --fa chm13v2.0.fa \
    --padding 150 \
    --show-axis \
    --show-coverage \
    --detail high \
    --overview-detail show \
    --bamsnap-src path/to/src/bamsnap_lrs/
