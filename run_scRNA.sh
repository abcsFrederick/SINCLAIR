species=$1
datatype=$2
outDir=$3
resume=$4
run_cellranger=$5
stubrun=$6

module load nextflow/23.04.1

echo "--Workflow started: $datatype"

if [[ $run_cellranger == "Y" ]]; then
    input="assets/input_manifest.csv"
    args="--run_cellranger Y"
else
    input="assets/input_manifest_cellranger.csv"
    args="--run_cellranger N"
fi

if [[ $stubrun == "Y" ]]; then
    input="assets/input_manifest.csv"
    args="$args -stub-run"
fi

if [[ $resume == "Y" ]]; then
    echo "----Resuming pipeline"
    nextflow run main.nf -resume \
    -entry $datatype \
    -profile biowulf \
    --input $input \
    --contrast assets/contrast_manifest.csv \
    --outdir $outDir \
    --species $species \
    $args

elif [[ $resume == "N" ]]; then
    nextflow run main.nf \
    -entry $datatype \
    -profile biowulf \
    --input $input \
    --contrast assets/contrast_manifest.csv \
    --outdir $outDir \
    --species $species \
    $args
fi

# with and without cellranger
# rm -r /data/sevillas2/scRNA_test/*; sh run_scRNA.sh hg19 GEX /data/sevillas2/scRNA_test N Y N; rm -r /data/sevillas2/scRNA_test/*; sh run_scRNA.sh hg19 GEX /data/sevillas2/scRNA_test N Y N

# with and without resume
# rm -r /data/sevillas2/scRNA_test/*; sh run_scRNA.sh hg19 GEX /data/sevillas2/scRNA_test Y N N; rm -r /data/sevillas2/scRNA_test/*; sh run_scRNA.sh hg19 GEX /data/sevillas2/scRNA_test Y N N

# with resume, with cellranger
#  sh run_scRNA.sh hg19 GEX /data/sevillas2/scRNA_test Y Y N