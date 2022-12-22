#!/bin/sh

hicFiles=$(ls /proj/phanstiel_lab/Data/processed/YY1P/hic/novaseq/hg38/*dietJuicerCore/YY1P*/*inter_30.hic)
sip_jar="/proj/phanstiel_lab/Software/SIP/SIP_HiCv1.6.2.jar"
chromSizes="/proj/phanstiel_lab/References/chromSizes/hg38_chromSizes.txt"
juicer_tools_jar="/proj/phanstiel_lab/Software/juicer/scripts/juicer_tools_1.14.08.jar"

for f in $hicFiles
do
 name=$(basename $f .hic)
 echo $name
 jid=`sbatch <<- SIP | egrep -o -e "\b[0-9]+$"
	#!/bin/sh
	#SBATCH -J $name
	#SBATCH -p general
	#SBATCH -n 1
	#SBATCH -N 1
	#SBATCH --mem=8G
	#SBATCH -t 4320
	#SBATCH -o %x_%j.out
	#SBATCH -e %x_%j.err

	java -jar $sip_jar \
	hic $f \
	$chromSizes \
	$name \
	$juicer_tools_jar \
	-g 2.0 -t 2000 -fdr 0.05

SIP`
	echo Submitted jobid: $jid
done
