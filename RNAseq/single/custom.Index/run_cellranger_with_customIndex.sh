####
## build index for custom 10x mm10 index
##     1. add Spp1-EGFR
##     2. add Spp1-tdTomato
##
## common code to run cellranger 
##     --libraries   
##     (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#libraries-csv)
##          fastqs,sample,library_type
##          path-to-GEX-fastqs,fastq-head-name,Gene Expression
##          path-to-FB-fastqs,fastq-head-name,Antibody Capture
##
##     --feature-ref
##     (https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#feature-ref)
##          id,name,read,pattern,sequence,feature_type          
##          totalSeqB301,sample_name,R2,5PNNNNNNNNNN(BC)NNNNNNNNN,ACCCACCAGTAAGAC,Antibody Capture
##          ...
##
## extract small bam files for IGV
##
####


#### 1. 
## new fasta
# vi Spp1-EGFP.fa
cp ../refdata-gex-mm10-2020-A/fasta/genome.fa ./genome.Spp1-EGFP.fa
cat Spp1-EGFP.fa >> genome.Spp1-EGFP.fa

## new gtf
# vi Spp1-EGFP.gtf
cp ../refdata-gex-mm10-2020-A/genes/genes.gtf ./genes.Spp1-EGFP.gtf
cat Spp1-EGFP.gtf >> genes.Spp1-EGFP.gtf

## code to submit for building index
# build_2020-A-Spp1-EGFP.slurm

#!/bin/bash

#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH -J sc10x_LYN
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 12
#SBATCH --mem-per-cpu=6GB
#SBATCH -o %j_slurm.out
#SBATCH -e %j_slurm.err

/pipelines/cellranger-7.0.0/bin/cellranger mkref \
--genome=refdata-gex-mm10-2020-A-Spp1-EGFP \
--fasta=genome.Spp1-EGFP.fa \
--genes=genes.Spp1-EGFP.gtf \
--nthreads=12

## common code to run cellranger 
# 
/pipelines/cellranger-7.0.0/bin/cellranger count \
--id=sample_output \
--transcriptome=path-to-new_index/refdata-gex-mm10-2020-A-Spp1-EGFP \
--libraries=libraries_sample.csv \
--feature-ref=feature.Seqx10.csv \
--include-introns false \
--localcores=12



#### 2.
## new fasta
# vi Spp1-tdt.fa
cp ../refdata-gex-mm10-2020-A/fasta/genome.fa ./genome.Spp1-tdt.fa
cat Spp1-tdt.fa >> genome.Spp1-tdt.fa

## new gtf
# vi Spp1-tdt-fixed.gtf
cp ../refdata-gex-mm10-2020-A/genes/genes.gtf ./genes.Spp1-tdt-fixed.gtf
cat Spp1-tdt-fixed.gtf >> genes.Spp1-tdt-fixed.gtf

## code to submit for building index
# build_2020-A-Spp1-tdt-fixed.slurm

#!/bin/bash

#SBATCH -p amd-ep2,amd-ep2-short
#SBATCH -J sc10x_LYN
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task 12
#SBATCH --mem-per-cpu=6GB
#SBATCH -o %j_slurm.out
#SBATCH -e %j_slurm.err

/pipelines/cellranger-7.0.0/bin/cellranger mkref \
--genome=refdata-gex-mm10-2020-A-Spp1-tdt-fixed \
--fasta=genome.Spp1-tdt.fa \
--genes=genes.Spp1-tdt-fixed.gtf \
--nthreads=12

## common code to run cellranger 

/pipelines/cellranger-7.0.0/bin/cellranger count \
--id=sample_output \
--transcriptome=path-to-new_index/refdata-gex-mm10-2020-A-Spp1-tdt-fixed \
--libraries=libraries_sample.csv \
--feature-ref=feature.Seqx10.csv \
--include-introns false \
--localcores=12


#### extract small bam files for IGV
##
cd ../sample_output/outs
mkdir check_exogene && cd check_exogene

## Spp1-EGFP
# raw Spp1 and exogenous Spp-EGFP could be both extracted by 'grep Spp1'
#    just change ref-genome to check them
#    'head -n 78': to keep bam header lines, may should check the number first
samtools view -@12 ../possorted_genome_bam.bam |grep Spp1 |\
cat <(samtools view -@12 -h ../possorted_genome_bam.bam |head -n 78) - |\
samtools view -@12 -bS |samtools sort -@12 - -o Spp1.bam
samtools index Spp1.bam

## Spp1-Tdt
samtools view -@12 ../possorted_genome_bam.bam |grep -Ew "Spp1|CreERT2|tdTomato" |\
cat <(samtools view -@12 -h ../possorted_genome_bam.bam |head -n 78) - |\
samtools view -@12 -bS |samtools sort -@12 - -o Spp1tdt.fixed.bam
samtools index Spp1tdt.fixed.bam


#### end