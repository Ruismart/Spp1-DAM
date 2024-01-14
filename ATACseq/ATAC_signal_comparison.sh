####
##
## a. whole genome correlation
## b. plot signal around differential peaks
## c. fasta for meme-chip
##
####


## a. 
multiBigwigSummary bins \
-p 4 \
-bs 5000 \
-b \
../../output/bigwig_file/Hapx_b2/Hapx_b2.CTRA_s2.nodup.rpkm.bw \
../../output/bigwig_file/Hapx_b2/Hapx_b2.CTRA_s3.nodup.rpkm.bw \
../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTN_s2.nodup.rpkm.bw \
../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTN_s3.nodup.rpkm.bw \
../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTP_s2.nodup.rpkm.bw \
../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTP_s3.nodup.rpkm.bw \
-l \
h2_CTRA_s2 \
h2_CTRA_s3 \
h2_TDTN_s2 \
h2_TDTN_s3 \
h2_TDTP_s2 \
h2_TDTP_s3 \
-o ./20230925_ATAC_LYN.filtered.correlation.npz 

#
plotCorrelation \
-in 20230925_ATAC_LYN.filtered.correlation.npz \
-o 20230925_ATAC_LYN.filtered.correlation.png \
--corMethod pearson \
-p scatterplot \
--removeOutliers


## b.
# save coordinates of differential results as bed files
#     used for GREAT analysis (http://great.stanford.edu)
#     column 4 could be added as chrN_start_end 
#     then GREAT associated genes get easier to target back to peaks
# get midpoint.bed
cat ../sorted.Cup789.bed |\
awk 'OFS="\t"{printf "%s\t%d\t%d\n" ,$1,$2+($3-$2)/2,$2+1+($3-$2)/2}' \
> sorted.Cup789.midpoint.bed

cat ../sorted.Pup44.bed |\
awk 'OFS="\t"{printf "%s\t%d\t%d\n" ,$1,$2+($3-$2)/2,$2+1+($3-$2)/2}' \
> sorted.Pup44.midpoint.bed

# deeptools
computeMatrix reference-point -S \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.CTRA_s2.nodup.rpkm.bw \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.CTRA_s3.nodup.rpkm.bw \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTN_s2.nodup.rpkm.bw \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTN_s3.nodup.rpkm.bw \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTP_s2.nodup.rpkm.bw \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTP_s3.nodup.rpkm.bw \
-R ./sorted.Cup789.midpoint.bed \
--skipZeros \
-o ./sorted.Cup789.mat.gz \
-p 6 \
-a 3000 -b 3000 \
--referencePoint center \
--samplesLabel \
CTRA_s2 \
CTRA_s3 \
TDTN_s2 \
TDTN_s3 \
TDTP_s2 \
TDTP_s3

computeMatrix reference-point -S \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.CTRA_s2.nodup.rpkm.bw \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.CTRA_s3.nodup.rpkm.bw \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTN_s2.nodup.rpkm.bw \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTN_s3.nodup.rpkm.bw \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTP_s2.nodup.rpkm.bw \
../../../output/bigwig_file/Hapx_b2/Hapx_b2.TDTP_s3.nodup.rpkm.bw \
-R ./sorted.Pup44.midpoint.bed \
--skipZeros \
-o ./sorted.Pup44.mat.gz \
-p 6 \
-a 3000 -b 3000 \
--referencePoint center \
--samplesLabel \
CTRA_s2 \
CTRA_s3 \
TDTN_s2 \
TDTN_s3 \
TDTP_s2 \
TDTP_s3

plotHeatmap \
-m ./sorted.Cup789.mat.gz \
-out ./Cup789.pdf \
--sortUsing sum \
--startLabel "Peak Start" \
--endLabel "Peak End" \
--xAxisLabel "" \
--regionsLabel "Peaks" \
--colorList 'white, #EE2124' \
--heatmapWidth 5 \
--heatmapHeight 20 \
--dpi 200 \
--sortUsingSamples 1 2

plotHeatmap \
-m ./sorted.Pup44.mat.gz \
-out ./Pup44.pdf \
--sortUsing sum \
--startLabel "Peak Start" \
--endLabel "Peak End" \
--xAxisLabel "" \
--regionsLabel "Peaks" \
--colorList 'white, #EE2124' \
--heatmapWidth 5 \
--heatmapHeight 20 \
--dpi 200 \
--sortUsingSamples 5 6

## c.
fa_file='/Shared/genomics/mouse/GRCm38_vM25/GRCm38.p6.genome.fa'

bedtools getfasta -fi ${fa_file} -bed \
../sorted.Cup789.bed \
-fo sorted.Cup789.fa

bedtools getfasta -fi ${fa_file} -bed \
../sorted.Pup44.bed \
-fo sorted.Pup44.fa

# submit to https://meme-suite.org/meme/tools/meme-chip
#     Maximum width: 30
#     database select:  
#         Human and Mouse(HOCOMOCO v11 FULL)
## end

