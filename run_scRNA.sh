species=$1
datatype=$2
resume=$3

module load nextflow

echo "--Workflow started: $datatype"

if [[ $resume == "Y" ]]; then
    echo "----Resuming pipeline"
    nextflow run main.nf -resume \
    -entry $datatype \
    -profile biowulf \
    --input assets/input_manifest.csv \
    --outdir /data/sevillas2/scRNA_test \
    --species $species \

elif [[ $resume == "N" ]]; then
    nextflow run main.nf \
    -entry $datatype \
    -profile biowulf \
    --input assets/input_manifest.csv \
    --outdir /data/sevillas2/scRNA_test \
    --species $species
fi

# sh run_scRNA.sh hg19 GEX Y