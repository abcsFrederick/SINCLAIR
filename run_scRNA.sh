species=$1
datatype=$2
resume=$3

module load nextflow


if [[ "$resume" == "Y" ]]; then
    echo "--Workflow started: $datatype"

    nextflow run main.nf -resume \
    -entry $datatype \
    -profile docker \
    --input assets/input_manifest.csv \
    --outdir /data/sevillas2/scRNA_test \
    --species $species \

elif [[ "$resume" == "N" ]]; then
    nextflow run main.nf \
    -entry $datatype \
    -profile docker \
    --input assets/input_manifest.csv \
    --outdir /data/sevillas2/scRNA_test \
    --species $species
fi

# sh run_scRNA.sh hg19 GEX Y