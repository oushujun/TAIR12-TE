#!/bin/bash
#SBATCH --job-name=rum_merge
#SBATCH --output=output/%x_output.txt
#SBATCH --error=errors/%x_error.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=100G
#SBATCH --time=08:00:00

# ------------------------------------
# Modules
module load anaconda3/2019.03


# ------------------------------------
# Main
source activate pannagram

P_SIM=80
P_COV=95
N_CNT=4

elements=("Mutator" "LINE" "hAT_TIR" "DNA_transposon" "LTR")

PATH_PAN="/users/anna.igolkina/work/tools/pan/pannagram/"

PATH_OUTPUT="/users/anna.igolkina/work/help/tair12/"
PATH_GFF_OUT="${PATH_OUTPUT}simsearch_${P_SIM}_${P_COV}_${N_CNT}_gff/"
mkdir -p "${PATH_GFF_OUT}"
file_gff_parent="${PATH_GFF_OUT}gff_merged.gff"
if [ -f  ${file_gff_parent} ]; then
    rm ${file_gff_parent}
fi


FILE_GFF="/users/anna.igolkina/work/help/tair12/data/GCA_028009825.2_Col-CC_genomic.fna.mod.EDTA.TEanno.clean.rename.ORFinfo.soloLTR.ATHILA.ONSEN.Evade.TRASH.update_AGI.rm_overlap.20240925_AI.gff3"

for element in "${elements[@]}"
do
	folder_output="${PATH_OUTPUT}simsearch_${P_SIM}_${P_COV}_${N_CNT}_${element}/"

	echo "Output folder: ${folder_output}"
	${PATH_PAN}inst/merge_hits.sh -i  ${FILE_GFF} \
	    -g /users/anna.igolkina/work/help/tair12/genome/GCA_028009825.2_Col-CC_genomic.fasta \
	    -o  ${folder_output} \
	    -s ${P_SIM} -c ${P_COV} -r 20 -n ${N_CNT} --keepblast -p ${element}

	echo ${file_gff_element}
	file_gff_element="${folder_output}gff_out/gff_merged.gff"
	if [ -f  ${file_gff_element} ]; then
	    cat ${file_gff_element} >> "${file_gff_parent}"
	fi

done

if [ -f  ${file_gff_parent} ]; then

	Rscript ${PATH_PAN}inst/merge/merge_05_comb_gff.R  \
			--path.out ${PATH_GFF_OUT} \
			--file.gff ${FILE_GFF} \
			--file.gff.parent ${file_gff_parent}

fi


conda deactivate






