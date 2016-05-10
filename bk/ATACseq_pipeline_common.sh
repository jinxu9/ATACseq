for file in `ls  ../Raw_data/merge_file/*_R1.fq.gz.trimmed.gz`
do
file2=`echo $file |sed 's/R1/R2/g'` 
echo $file
echo $file2
name=`basename $file`
output=`echo $name| awk -F_ '{print $1}'`
echo $output
bowtie2   -p 15  -x /home/jinxu/DB/mmu10/mm10_allchr_bowtie2Index/mm10 -1 $file -2 $file2 -S  $output.sam 1>$output.mapping.log 2>$output.mapping.err 
done


for file in `ls *.sam`
do
echo $file
file2=`echo $file|sed 's/\.sam//g'`
echo $file2
{
awk '$3!="chrM"' $file | samtools view -S -b -F 4  -   -o $file2.bam  
samtools sort $file2.bam $file2.sort
}&
done 

wait

# merge.bam has beedn sorted and rm chrM 
for file in `ls *.bam`
do
{
/home/jinxu/bin/genomeCoverageBed -bg -split -ibam $file -g /home/jinxu/DB/mmu10/mm10_allchr_bowtie2Index/mm10.all.chrSize > $file.bedGraph 
/home/jinxu/bin/norm_bedGraph.pl $file.bedGraph $file.bedGraph.norm
/home/jinxu/bin/bedGraphToBigWig $file.bedGraph.norm   /home/jinxu/DB/mmu10/mm10_allchr_bowtie2Index/mm10.all.chrSize   $file.bedGraph.norm.bw 
}& 
done

wait

for file in `ls *.bam`
do
{
bam2bed  $file # using good paired-end reads only.
bedfile=`echo $file|sed 's/bam/bed/g'`
outfile=`echo $file |sed 's/\.merge\.bam//g'`
path=`pwd`
macs2  callpeak -t  $file.bed.rmchrM -f BED  -g mm --outdir $path  -q 0.01 -n $outfile --nomodel  --shift 0  
} &
done


