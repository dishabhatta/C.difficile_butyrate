## RNA-seq pipeline

## Adapted from Rooksie's pipelines in the workshop

## DB 
## 10.28.2022

### Files: RNA-seq files (fastq) + CD630 assembly files
### Experiment design: A1,2,3 - early BHI; A7, 8, 9- late BHI; B1,2,3 - early Butyrate; B7,8,9- late Butyrate

####  Process: fastqc, trim-galore, (multiqc both fastqc results), metaphlan, kneaddata, 


## Organising the files and folders:
### We whole genome sequenced the CD630 we used to conduct the growth curves
### CD630 wgs working out of references 
### GCA_000009205.2_ASM920v2_genomic.fna and GCF_000932055.2_ASM93205v2_genomic.fna are the official C. difficile WGS from NCBI
### RNA sequences unzipped from  (zip bomb error)
		jar -xf ./Cd630.But01_RNASeq_2022.02.08.zip 

### 3 sets of sequences: A1, A2, A3 are triplicates of without butyrate early log; A7, A8, A9 are triplicates of without butyrate late log; B1, B2, B3 are triplicates of with butyrate early log and B7, B8, B9 are triplicates of with butyrate late log.

### full workflow for RNAseqs
#### fastqc for quality and trim galore to trim for the raw seqs

```
#!/bin/sh
#PBS -N pipeline_rnaseqs 
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=251gb:interconnect=1g,walltime=100:00:00

module add fastqc/0.11.9
module add trimgalore/0.6.7

cd /scratch1/dbhatta/RNA_seq/seqs/transcripts/

for i in `cat ./sample_id`; do
		fastqc -t 14 ${i}_R1_001.fastq ${i}_R2_001.fastq -o ${i}_bef_fastqc
		trim_galore --paired -j 14 --phred33 --fastqc --length 36 -o ./${i}_trim ${i}_R1_001.fastq ${i}_R2_001.fastq
done
qstat -xf $PBS_JOBID

```
### using multiqc/1.12 on before and after trimming
### compile into excel for stats


### references
#### This one is the wgs our lab did on our strain, so I started by assembling and annotating it
```
#!/bin/sh
#PBS -N pipeline_CD630 
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=251gb:interconnect=1g,walltime=100:00:00 

module add samtools/1.13
module add prokka/1.14.5
module add bowtie2/2.4.5
module add spades/3.15.0

cd /scratch1/dbhatta/RNA_seq/reference/

gunzip *.fastq.gz
for i in `cat ./sample_id`; do
        ~/anaconda3/bin/trim_galore --fastqc --trim-n --paired ${i}_R1_001.fastq ${i}_R2_001.fastq --path_to_cutadapt ~/anaconda3/bin/cutadapt --output_dir ./${i}_trim/
        spades.py -o ./${i}_spades/ -1 ./${i}_trim/${i}_R1_001_val_1.fq -2 ./${i}_trim/${i}_R2_001_val_2.fq -k 55,77,127 --careful -t 8
        ~/anaconda3/bin/bowtie2-build ./${i}_spades/contigs.fasta ./${i}_spades/${i}_index #### I had bowtie2; you can use the one Rooksie has
        ~/anaconda3/bin/bowtie2 -x ./${i}_spades/${i}_index -1 ${i}_R1_001.fastq -2 ${i}_R2_001.fastq -S ./${i}_spades/sam_${i}.sam
        samtools view -b ./${i}_spades/sam_${i}.sam -o ./${i}_spades/bam_${i}.bam
        samtools sort ./${i}_spades/bam_${i}.bam -o ./${i}_spades/sort_${i}.bam
        samtools coverage ./${i}_spades/sort_${i}.bam | awk '{sum += $6}END{ print "Average = ", sum/NR}' > ./${i}_spades/${i}_average.txt
        prokka --outdir ./${i}_spades/prokka_${i} --prefix ${i} ./${i}_spades/contigs.fasta
        awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' ./${i}_spades/prokka_${i}/${i}.ffn > ./${i}_spades/prokka_${i}/${i}_sl.ffn
        grep -A1 "16S ribosomal RNA" ./${i}_spades/prokka_${i}/${i}_sl.ffn > ${i}_16S.txt
done

qstat -xf $PBS_JOBID

```

### originally, kneaddata along with metaphlan was used for removing contamination in these sequences; but that wasn't working when I tried it so I found a different way of removing contimnants by mapping them to the contam seqs and then collecting everything that did not align with the contam ref

### Fixing contamination

### Karken and metaphlan showed that there is a contamination in A7, A8 and A9. A9 had the max contam.

```
### these were trials on A9 whre most contam was observed: went down from 8.58724 to 1.9 in metaphlan
rsem-prepare-reference --bowtie2 ./contam_csel.fasta CONTAM_REF

rsem-calculate-expression --bowtie2 -p 15 --paired-end --strandedness none ./A9_S66_R1_001_val_1.fq A9_S66_R2_001_val_2.fq /scratch1/dbhatta/RNA_seq/reference/rsem_contam/CONTAM_REF A9_trial
samtools view -b -f 4 A9_trial.transcript.bam > unmapped_A9_trial.bam
samtools sort -n -o unmapped_A9_trial_sort.bam unmapped_A9_trial.bam
bedtools bamtofastq -i ./unmapped_A9_trial_sort.bam -fq contamfree_A9_trial_R1.fq -fq2 contamfree_A9_trial_R2.fq
metaphlan contamfree_A9_trial_R1.fq,contamfree_A9_trial_R2.fq -o profile_test_A9.txt --bowtie2out bowtie_A9.bz2 --nproc 14 --input_type fastq --read_min_len 49

metaphlan contamfreeF2_A9_trial_R1.fq,contamfreeF2_A9_trial_R2.fq -o profileF2_test_A9.txt --bowtie2out bowtieF2_A9.bz2 --nproc 14 --input_type fastq --read_min_len 49



## to do this for all strains


#!/bin/sh
#PBS -N contam_removal
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=62gb:interconnect=fdr,walltime=72:00:00

module add rsem/1.3.3
module add perl/5.36.0
module add samtools/1.15.1
module add metaphlan/4.0.2
module add bedtools/2.30.0

cd /scratch1/dbhatta/RNA_seq/reference/rsem_contam/contam_rsem/

for i in `cat sample_id`; do

        mkdir ${i}_contam_rsem
        rsem-calculate-expression --bowtie2 -p 15 --paired-end  --bowtie2-sensitivity-level very_sensitive --strandedness none ./${i}_R1_001_val_1.fq ./${i}_R2_001_val_2.fq /scratch1/dbhatta/RNA_seq/reference/rsem_contam/ ./${i}_contam_r$
        samtools view -b -f 4 ./${i}_contam_rsem/${i}_contam.transcript.bam > ./${i}_contam_rsem/unmapped_${i}.bam
        samtools sort -n -o ./${i}_contam_rsem/unmapped_${i}_sort.bam ./${i}_contam_rsem/unmapped_${i}.bam
        bedtools bamtofastq -i ./${i}_contam_rsem/unmapped_${i}_sort.bam -fq ./${i}_contam_rsem/contamfree_${i}_R1.fq -fq2 ./${i}_contam_rsem/contamfree_${i}_R2.fq
        metaphlan ./${i}_contam_rsem/contamfree_${i}_R1.fq,./${i}_contam_rsem/contamfree_${i}_R2.fq -o ./${i}_contam_rsem/profile_${i}_contam.txt --bowtie2out ./${i}_contam_rsem/bowtie2_aftercontam_${i}.bz2 --nproc 14 --input_type fastq $

done

qstat -xf $PBS_JOBID

```

### back to seqs

### prep the cdiffcile reference for an alignment to cdiff

```
rsem-prepare-reference --bowtie2 ./contigs.fasta ./CD630_reference

#!/bin/sh
#PBS -N rsem_final
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=62gb:interconnect=fdr,walltime=72:00:00

module add rsem/1.3.3
module add perl/5.36.0

cd /scratch1/dbhatta/RNA_seq/reference/rsem_ref/rsem_toref

for i in `cat sample_id`; do
	mkdir ${i}_rsem_ref
	rsem-calculate-expression --bowtie2 -p 15 --paired-end --strandedness none ./contamfree_${i}_R1.fq ./contamfree_${i}_R2.fq /scratch1/dbhatta/RNA_seq/reference/rsem_ref/CD630_reference ./${i}_rsem_ref/${i}_final
	rsem-plot-model ./${i}_rsem_ref/${i}_final ./${i}_rsem_ref/${i}_final_diagnostic.pdf
	
done
qstat -xf $PBS_JOBID

```

### to note down the RNA prelim stats: from samtools, get into ./${i}_rsem_ref directory and run:

	samtools flagstat ./B3_S69_final.transcript.bam 

### for bowtie2 rsem stats: read the .log file 

### collect .gff and .bam files in a single folder run subread/featurecounts

```
module load samtools/1.15.1
module load subread/2.0.3


featureCounts -a CD630_S243.gff -t CDS -g ID -Q 0 -p -P -d 50 -D 600 -B -o All.gene.counts.final.txt -T 15 -s 0 *.transcript.bam
# outputs: All.gene.counts.final.txt & All.gene.counts.final.txt.summary

featureCounts -a GCA_000009205_2_ASM920v2_genomic.gff -t gene -g ID -Q 0 -p -P -d 50 -D 600 -B -o All.gene.counts.final_GCA9205.txt -T 15 -s 0 *.transcript.bam
featureCounts -a GCF_000932055.2_ASM93205v2_genomic.gff -t gene -g ID -Q 0 -p -P -d 50 -D 600 -B -o All.gene.counts.final.txt -T 15 -s 0 *.transcript.bam


```

### Rest is in RDesktop

### I repeated the entire process from rsem-prepare-reference if I used the NCBI references for a better featureCounts profiles: GCA_000009205.2_ASM920v2_genomic.fna and GCF_000932055.2_ASM93205v2_genomic.fna; their stats have been separately saved

### our CD630 .gff doesnot work with -t gene

### Finally GCA_000009205.2_ASM920v2_genomic.fna was used in the final figures




## Blast2go, Interproscan and EGGNOG-mapper
### Over-representation of KEGG and gene set enrichment analysis of KEGG/KO values for supplementary

### the faa was annotated with EGGNOG-MAPPER for ko and KEGG
### additional annotation can be done using Blast2GO and blastp and interproscan
### Blast2go was downloaded on the local machine

```
Interproscan

#!/bin/sh
#PBS -N interproscan
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=62gb:interconnect=1g,walltime=100:00:00

module load interproscan/5.59-91.0

cd /scratch1/dbhatta/RNA_seq/reference/blast2go_files/

interproscan.sh -i CD630_S243.faa -t p -d interpro_results -dp -pa -appl Pfam,PANTHER,CDD,Hamap,PRINTS,PIRSR,TIGRFAM --goterms --iprlookup

qstat -xf $PBS_JOBID

BLAST

#!/bin/sh
#PBS -N blast
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=62gb:interconnect=1g,walltime=100:00:00

source activate blast

cd /scratch1/dbhatta/RNA_seq/reference/blast2go_files/

blastp -num_threads 12 -max_target_seqs 10 -evalue 0.001 -outfmt 14 -db /zfs/gcl/share/dblibs/ncbi-blast/db/nr -query CD630_S243.faa -out CD_blast2go.xml

qstat -xf $PBS_JOBID


EGGNOG-MAPPER


#!/bin/sh
#PBS -N eggnogmapper
#PBS -j oe 
#PBS -m abe
#PBS -l select=1:ncpus=16:mem=251gb:interconnect=1g,walltime=100:00:00

source activate eggnogmapper

cd /scratch1/dbhatta/RNA_seq/reference/eggnogmap_cd630/gca9205

emapper.py -i GCA_000009205.2_ASM920v2_protein.faa --itype proteins --data_dir /scratch1/dbhatta/eggnog_mapper_db --output_dir /scratch1/dbhatta/RNA_seq/reference/eggnogmap_cd630/gca9205 -o cd630_gca9205

qstat -xf $PBS_JOBID
```

























