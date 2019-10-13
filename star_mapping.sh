scriptdir=$(dirname $(realpath -s $0))
echo $scriptdir
while getopts 'o:s:d:' c; do
	case $c in
		o) optional=$OPTARG ;;
		s) sampfile=$OPTARG ;;
		d) workdir=$OPTARG ;;
	esac
done 

cd $workdir 

mkdir STAR_output
mkdir gene_counts

tail -n+2 $sampfile | while read line; do
	sample=$(echo $line | cut -d ' ' -f 1) 
	
	command="STAR --runThreadN 24 --genomeDir /ebio/ecnv_projects/common_resourses/data/reference_genome/GRCz11/ --readFilesIn $sample.1.fastq.gz $sample.2.fastq.gz --readFilesCommand gunzip -c --sjdbInsertSave All --quantMode GeneCounts --outReadsUnmapped Fastq --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --outFileNamePrefix STAR_output/$sample --genomeLoad LoadAndKeep --limitBAMsortRAM 20000000000 $optional"
	
	echo $command
	
	$command

	tail -n +5 STAR_output/${sample}ReadsPerGene.out.tab > gene_counts/${sample}.count

done

Rscript $scriptdir/prepare_matrix.R -e "setwd('$workdir')"
