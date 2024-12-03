#!/bin/bash
PATH_SRC="src_merge/"
# ----------------------------------------------------------------------------
#            ERROR HANDLING BLOCK
# ----------------------------------------------------------------------------

# Exit immediately if any command returns a non-zero status
set -e

# Keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG

# Define a trap for the EXIT signal
trap 'catch $?' EXIT

# Function to handle the exit signal
catch() {
    # Check if the exit code is non-zero
    if [ $1 -ne 0 ]; then
        echo "\"${last_command}\" command failed with exit code $1."
    fi
}

# ----------------------------------------------------------------------------
#             FUNCTIONS
# ----------------------------------------------------------------------------

# Function to display a short help message
print_usage() {
    cat << EOF

Usage:
${0##*/} -i F_GFF -g F_GENOME -o PATH_OUT \
		[-s SIMILARITY] [-d DISTANCE] [-p PATTERNS] [-n COPY_NUM] [-k] [-h]

Options:
    -h, --help         Display this help message and exit.
    
    -i, --file_gff F_GFF         Path to the input gff file with TEs.
    -g, --file_genome F_GENOME   Path to the genome file that will be used for the analysis.
    -o, --path_out PATH_OUT      Path to the directory where all results and intermediate files will be saved.
    
    -s, --sim  SIMILARITY        Similarity cutoff value used in the analysis. Default: 90.
    -c, --coverage  COVERAGE     Similarity cutoff value used in the analysis. Default: equal to sim.
    -d, --distance  DISTANCE     Distance between two hits. Default: 1000.
    -p, --patterns PATTERNS      Patterns of repeats to analyse. Default: LTR.
    -n, --copy_number COPY_NUM   Minimum number of copies per genome. Default: 4.
    -r, --max_rounds MAX_ROUND   Maximum number of rounds of merging. Default: 20.
    -k, --keepblast              Keep the intermediate BLAST files
    -v, --visualise              Create figures for final merges


EOF
}

# ----------------------------------------------------------------------------
#            PARAMETERS: parsing
# ----------------------------------------------------------------------------

unrecognized_options=()

sim_cutoff=90
# coverage=85
distance=1000
patterns="LTR"
copy_number=4
max_rounds=20
keepblast=""
visualise="FALSE"

PATH_PAN=""  # path to the pannagram folder

while [ $# -gt 0 ]
do
    # echo $1
    case $1 in
        -h | --help) print_usage; exit ;;
        
        # Required
        -i | --file_gff)     file_gff=$2;    shift 2 ;;  
		-g | --file_genome)  file_genome=$2; shift 2 ;;  
        -o | --path_out)     path_out=$2;    shift 2 ;;  
		
		# Optional
		-s | --sim) 		 sim_cutoff=$2;  shift 2 ;;
        -c | --coverage)     coverage=$2;    shift 2 ;;
		-d | --distance) 	 distance=$2;    shift 2 ;;  
        -p | --patterns)     patterns=$2;    shift 2 ;;
		-n | --copy_number)  copy_number=$2; shift 2 ;;
        -r | --max_rounds)   max_rounds=$2;  shift 2 ;;
        --pan )              PATH_PAN=$2;    shift 2 ;;

        -k | --keepblast)    keepblast=' -keepblast ';    shift 1 ;;
        -v | --visualise)    visualise="TRUE";            shift 1 ;;

        *) unrecognized_options+=("$1"); shift 1 ;;
    esac
done

# Output of Unrecognized Parameters
if [[ ${#unrecognized_options[@]} -gt 0 ]]; then
    print_usage
    echo "Unrecognized options:"
    for option in "${unrecognized_options[@]}"; do
        echo "$option"
    done
    exit 1
fi


# Check if coverage parameter is provided. If not - set qeual to sim
if [ -z "$coverage" ]; then
    coverage=${sim_cutoff}
fi

echo "Coverage: ${coverage}"
echo "Similarity: ${sim_cutoff}"

# ----------------------------------------------------------------------------
#            PARAMETERS: checking
# ----------------------------------------------------------------------------


# Check for required parameters
if [ -z "$file_gff" ] || [ -z "$file_genome" ] || [ -z "$path_out" ]; then
    echo "Error: Missing required parameters (-i, -g, -o)."
    print_usage
    exit 1
fi

# Check the output directory
if [ -d "$path_out" ]; then   # Exist ?
    if [ "$(ls -A "$path_out")" ]; then  # Empty ?
        echo "Warning: The directory '$path_out' exists and is not empty."
    fi
else
    mkdir -p "$path_out"  

    # if the directory was not created successfully
    if [ ! -d "$path_out" ]; then
        echo "Error: Failed to create the directory '$path_out'."
        exit 1
    fi
fi

# Ensure $path_out ends with /
if [ "${path_out: -1}" != "/" ]; then
    path_out="$path_out/"
fi


# path_out_gff="${path_out}simsearch_${sim_cutoff}_${coverage}_${copy_number}_gff/"
# mkdir -p "${path_out_gff}"
file_gff_parent="${path_out}gff_merged_only.gff"
if [ -f  ${file_gff_parent} ]; then
    rm ${file_gff_parent}
fi

# ----------------------------------------------------------------------------
#            MAIN
# ----------------------------------------------------------------------------

# ----------------------------------------
# Extract patterns to analyse
IFS=',' read -ra elements <<< "$patterns"

echo "* The following types will be analyzed:"
for element in "${elements[@]}"; do
    echo "  $element"
done

for element in "${elements[@]}"
do
    path_out_elem="${path_out}simsearch_${sim_cutoff}_${coverage}_${copy_number}_${element}/"
    mkdir -p "${path_out_elem}"

    # echo "Output folder: ${path_out_elem}"

    # ----------------------------------------
    # Read the gff file, and get sequences

    echo "* Stage: Extract candidates for merging for type ${element}."

    file_merged_seqs="${path_out_elem}merged_seqs_1.fasta"

    Rscript ${PATH_SRC}merge_01_extract_hits.R \
            --file.gff=${file_gff} \
            --file.genome=${file_genome} \
            --file.seqs=${file_merged_seqs} \
            --patterns=${element} \
            --len.gap=${distance}

    # ----------------------------------------
    # Simrearch and Merge

    echo "* Stage: Search for copies in the genome over several rounds."

    file_merged_seqs_fixed="${path_out_elem}merged_seqs_fixed.txt"

    for ((i=1; i<=max_rounds; i++)) ; do
        echo "  == Round ${i} =="
    	file_merged_seqs="${path_out_elem}merged_seqs_${i}.fasta"

    	if [ ! -f "${file_merged_seqs}" ]; then
            break
        fi
        
        path_simsearch="${path_out_elem}simseqrch_seqs_${i}/"

    	# Run simsearch

    	simsearch \
        -in_seq ${file_merged_seqs}    \
        -on_genome ${file_genome} \
        -out "${path_out_elem}simseqrch_seqs_${i}/" \
        -sim ${sim_cutoff} \
        -cov ${coverage} \
        ${keepblast}

    	# Get Collapsed sequences - neighbours only

        file_merged_seqs_next="${path_out_elem}merged_seqs_$((i+1)).fasta"

        file_cnt=$(find ${path_simsearch} -type f -name "*${sim_cutoff}_${coverage}.cnt")

        if [ -z "$file_cnt" ]; then
            echo "Error: No file ending with ${sim_cutoff}_${coverage}.cnt found." >&2
            exit 1
        elif [ $(echo "$file_cnt" | wc -l) -ne 1 ]; then
            echo "Error: More than one file matching the pattern *${sim_cutoff}_${coverage}.cnt found." >&2
            exit 1
        fi

        echo "* Analyse counts.."

        Rscript ${PATH_SRC}merge_02_new_hits.R \
            --file.cnt ${file_cnt} \
            --file.genome ${file_genome} \
            --file.seqs ${file_merged_seqs_next} \
            --file.fix ${file_merged_seqs_fixed} \
            --copy.number=${copy_number}

    done

    # echo ${file_merged_seqs_fixed}
    if [ ! -s ${file_merged_seqs_fixed} ];  then
        continue
    fi

    # ----------------------------------------
    # Get the merges and search for them one more time
    echo  "* Stage: Combine merges."

    file_fix_seqs="${path_out_elem}seqs_fix.fasta"

    Rscript ${PATH_SRC}merge_03_get_hits.R \
        --path.out ${path_out_elem} \
        --file.fix ${file_merged_seqs_fixed} \
        --file.fix.seqs=${file_fix_seqs}

    if grep -q "^>" ${file_fix_seqs}; then
        simsearch \
            -in_seq ${file_fix_seqs}    \
            -on_genome ${file_genome} \
            -out "${path_out_elem}simseqrch_seqs_fix/" \
            -sim ${sim_cutoff} \
            -cov ${coverage} \
            ${keepblast}
    else
        echo "Nothing to merge"
    fi

    # ----------------------------------------
    # Visualise and get information
    echo  "* Stage: Get pre-gff file, only with parents."

    Rscript ${PATH_SRC}merge_04_visualisation.R \
            --path.out  ${path_out_elem} \
            --file.genome ${file_genome} \
            --file.gff=${file_gff} \
            --plot ${visualise}

    # ----------------------------------------
    # Merge all outputs together

    # echo ${file_gff_element}
    file_gff_element="${path_out_elem}gff_out/gff_merged.gff"
    if [ -f  ${file_gff_element} ]; then
        cat ${file_gff_element} >> "${file_gff_parent}"
    fi

done

# ----------------------------------------
# Final GFF
echo  "* Stage: Final gff-file."

if [ -f  ${file_gff_parent} ]; then

    Rscript ${PATH_SRC}merge_05_comb_gff.R  \
            --path.out ${path_out} \
            --file.gff ${file_gff} \
            --file.gff.parent ${file_gff_parent}
fi








