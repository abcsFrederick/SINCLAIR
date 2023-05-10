#!/bin/bash

flag=$1

# nextflow run main.nf -entry GEX --input assets/input_manifest.csv --outdir /data/sevillas2/scRNA_test

###################################################################################################################
# GEX
###################################################################################################################
if [[ $flag == "gex_prep" ]]; then
    test_dir="/data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs"

    # untar file
    tar -xvf \
    /data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/data_storage/WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs.tar \
    -C $test_dir
    
    # subset reads
    for f in $test_dir/*.gz; do
        filename=`echo $f | awk -F"/" '{print $NF}'`
        echo "--$filename"
        zcat $f | head -n 500000 > $test_dir/sub_$filename
    done

    # move into subdirs
    out_1="$test_dir/WB_Lysis_1"
    if [[ ! -d $out_1 ]]; then mkdir -p $out_1; fi
    out_2="$test_dir/WB_Lysis_2"
    if [[ ! -d $out_2 ]]; then mkdir -p $out_2; fi

    mv $test_dir/sub*L001* $out_1
    mv $test_dir/sub*L002* $out_2

    # cleanup
    rm $test_dir/*gz

    # review
    tree $test_dir
fi

if [[ $flag == "gex_count" ]]; then

    output="/data/sevillas2/scRNA/gex"
    genome_dir="/fdb/cellranger/refdata-cellranger-3.0.0/hg19"
    SLURM_CPUS_PER_TASK=16
    test_dir="/data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir"
    sample_manifest="$test_dir/gex_samplesheet.csv"

    #This command is run for each sample
    #FASTQ directory should contain paired-end read \
    #files with some unknown reads in the format R1, R2, or I1. There are normally 3 or 4 sets per sample

    #The genome reference directory should contain 4 directories \
    #(fasta, genes, pickle, star) and a reference.json. 
    #Refer to /fdb/cellranger/refdata-cellranger-3.0.0/mm10 for an example

    # cell ranger doesn't allow for output redirection; must move to output dir location
    #https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ct
    
    # remove header
    tail -n +2 $sample_manifest > $output/noheader_manifest.csv
    echo "\n" >> $output/noheader_manifest.csv

    sh_dir="$output/sh"
    if [[ ! -d $sh_dir ]]; then mkdir $sh_dir; fi

    submit_batch(){
		sbatch --cpus-per-task=32 --verbose \
		--output=$sh_dir/%j.out \
		--mem=150g --gres=lscratch:450 --time 12:00:00 \
		--error=$sh_dir/%j.err $sh_file
    }

    while IFS="," read -r sampleID inputDir
    do
        if [[ ! -d $output/$sampleID ]]; then mkdir -p $output/$sampleID; fi
        echo "--$sampleID"
        sh_file="$sh_dir/sh_${sampleID}_gex.sh"

        echo "#!/bin/sh
        
        module load cellranger/7.1.0
        cd $output/$sampleID

        cellranger count \
            --id $sampleID \
            --fastqs $inputDir \
            --transcriptome=$genome_dir \
            --localcores=$SLURM_CPUS_PER_TASK \
            --localmem=32 \
            --jobmode=slurm \
            --maxjobs=20" > $sh_file

        # submit sh
        # cat $sh_file
        submit_batch $sh_file
    done < $output/noheader_manifest.csv
fi