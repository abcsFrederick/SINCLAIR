/*
===================================================================
Assign Local Modules
===================================================================
*/
include {SAMPLESHEET_CHECK      } from '../../modules/local/samplesheet_check'

workflow INPUT_CHECK {
    take:
        samplesheet
    
    main:
        SAMPLESHEET_CHECK (samplesheet)
        .csv
        .splitCsv(heater:true,sep:',')
        .map {create_fastq_channels(it)}
        .set{reads}
    
    emit:
        reads
        valid_samplesheet = SAMPLESHEET_CHECK.out.csv
}

def create_fastq_channel(LinkedHashMap row){
    def meta=[:]
        meta.id             = row.sample
        meta.single_end     = row.single_end.toBoolean()

    def array = []
        if (!file(row.fastq_1).exists()){
            exit1, "ERROR: missing FASTQ file R1\n${row.fastq_1}"
        }
        if (meta.single_end){
            array = [meta,[file(row.fastq_1)]]
        } else {
            if(!file(row.fastq_2).exists()){
                exit1, "ERROR: missing FASTQ file R2\n${row.fastq_2}"
            }
            array=[meta,[file(row.fastq_1),file(row.fastq_2)]]
        }
    return array
}