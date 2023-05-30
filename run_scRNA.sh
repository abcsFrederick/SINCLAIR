species=$1
datatype=$2
resume=$3
args=$4

module load nextflow/23.04.1

echo "--Workflow started: $datatype"

if [[ $resume == "Y" ]]; then
    echo "----Resuming pipeline"
    nextflow run main.nf -resume \
    -entry $datatype \
    -profile biowulf \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --outdir /data/sevillas2/scRNA_test \
    --species $species \
    $args

elif [[ $resume == "N" ]]; then
    nextflow run main.nf \
    -entry $datatype \
    -profile biowulf \
    --input assets/input_manifest.csv \
    --contrast assets/contrast_manifest.csv \
    --outdir /data/sevillas2/scRNA_test \
    --species $species \
    $args
fi

# sh run_scRNA.sh hg19 GEX Y
# sh run_scRNA.sh hg19 GEX N
# sh run_scRNA.sh hg19 GEX Y -stub-run
# sh run_scRNA.sh hg19 GEX N -stub-run
