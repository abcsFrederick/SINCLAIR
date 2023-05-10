#!/usr/bin/env python

# GEX
# checks that all samples in INPUTDIR follow the format
# [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
# This script is based on the example at: 
# https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

# Testing 
# bin/check_samplesheet.py assets/input_manifest.csv /data/sevillas2/tmp/scRNA/project

import os
import sys
import errno
import argparse
import gzip
import re
from os import listdir
from os.path import isfile, join

def parse_args(args=None):
    Description = "Samplesheet file - check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)

def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,INPUTDIR
    WB_Lysis_1,/data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/WB_Lysis_1
    WB_Lysis_2,/data/CCBR_Pipeliner/Pipelines/TechDev_scRNASeq_Dev2023/test_dir/WB_Lysis_Granulocytes_3p_Introns_8kCells_fastqs/WB_Lysis_2

    For an example see:
    TO DO
    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:
        ###################################################################################################
        # ALL: Check manifest
        ###################################################################################################
        ## Check header
        MIN_COLS = 2
        HEADER = ["sample", "data_type","read_category","input_dir"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            # split lines
            SAMPLE, DATATYPE, READTYPE, INPUTDIR = lspl[: len(HEADER)]
            
            ## Check sample name entries
            SAMPLE = SAMPLE.replace(" ", "_")
            if not SAMPLE:
                print_error("Sample entry has not been specified!", "Line", line)

            ###################################################################################################
            # GEX: check manifest
            ###################################################################################################
            if DATATYPE == "gex":
                ## Check input dir exists
                if os.path.exists(INPUTDIR):
                    # check all files in input
                    onlyfiles = [f for f in listdir(INPUTDIR) if isfile(join(INPUTDIR, f))]
                    sample_list=list()
                        
                    # for every tile
                    for fileid in onlyfiles:

                        # ensure that the file follows the structure
                        ## [Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq.gz
                        # create sample list - will check all samples are the same name
                        sampleID=fileid.split("_S1")[0]
                        sample_list.append(sampleID)

                        # check S1_L00 is within file
                        sID="S1_L00" in fileid
                        if sID == False:
                            print_error("Input file must include S1_L00 in name:", fileid)

                        # spit after S1_L00
                        postID=fileid.split("_L00")[1]
                        
                        # check lane, readype and extension
                        laneID=postID.split("_")[0]
                        if not laneID.isnumeric():
                            print_error("Input file laneID must be numeric: %s" %fileid)
                            
                        accepted_read_type=["R1","R2","I1","I2"]
                        read_type=postID.split("_")[1]
                        if not read_type in accepted_read_type:
                            print_error("Input file read_type must be R1/R2/I1/I2 but %s was given" %read_type)
                            
                        ext=fileid.split("001.")[1]
                        if ext != "fastq.gz":
                            print_error("Input file extension must be .fastq.gz: %s" %ext)
                    
                    # ensure all samples names are the same
                    if(len(set(sample_list)) != 1):
                        print_error("All input files within given dir must have the same sample ID %s" %sample_list)

                ## Auto-detect paired-end/single-end
                sample_info = []  ## [DATATYPE, READTYPE, INPUTDIR]
                if READTYPE == "paired":  ## Paired-end short reads
                    sample_info = ["gex","0", INPUTDIR]
                elif READTYPE == "single":  ## Single-end short reads
                    sample_info = ["gex","1", INPUTDIR]
                else:
                    print_error("Invalid combination of columns provided!", "Line", line)

            ###################################################################################################
            # ALL: Create SAMPLE DICT
            ###################################################################################################
            ## Create sample mapping dictionary = { sample: [ single_end, INPUTDIR ] }
            if SAMPLE not in sample_mapping_dict:
                sample_mapping_dict[SAMPLE] = {}
                sample_mapping_dict[SAMPLE][DATATYPE] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[SAMPLE][DATATYPE]:
                    print_error("Samplesheet contains duplicate sample names!", "Line", line)
                else:
                    sample_mapping_dict[SAMPLE][DATATYPE].append(sample_info)

    ###################################################################################################
    # ALL: Create OUTPUT
    ###################################################################################################
    # Write validated samplesheet with appropriate columns for 
    ## cellranger input
    ## hash creation
    if len(sample_mapping_dict) > 0:
        for sample in sorted(sample_mapping_dict.keys()):
            for dt in sorted(sample_mapping_dict[sample].keys()):
                # GEX samplesheet
                if dt=="gex":
                    fname=file_out + "_gex_samplesheet.csv"
                    if not os.path.isfile(fname):
                        with open(file_out, "w") as fout:
                            fout.write(",".join(["sample","gex_input_dir"]) + "\n")
                            fout.close()
                    with open(fname, 'a+') as fout:
                        print(sample_mapping_dict[sample][dt])
                        for dt, rt, idir in sample_mapping_dict[sample][dt]:
                            fout.write(sample + "," + idir + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))

def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())