outfile="RaiKO-ATACseq"

peaks=`ls *.narrowPeaks`

mergeZINBAPeaks  $outfile.peakbed $peaks
zinbapeakbed2txt.pl /home/jinxu/DB/mmu10/annotation/refGene.txt $outfile.peakbed  $outfile.mergedpeaks_anno.txt 

for file in `ls *.merge.sort.bam`
do
{

rnaexp_raw_fast.pl   $outfile.mergedpeaks_anno.txt   /home/jinxu/DB/mmu10/Combin_129_CASTEiJ/129CAST_genome_index/129S1_CASTEiJ_mm10.size  $file.bedGraph $file.exp
}&
done

wait

files=`ls *.exp`
paste $files >$outfile.merge_raw_exp
select_ataccolumn.pl  $outfile.merge_raw_exp  $outfile.merge_Peak_exp.txt



