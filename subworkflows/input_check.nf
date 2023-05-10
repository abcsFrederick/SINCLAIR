/*
===================================================================
Assign Local Modules
===================================================================
*/
include {SAMPLESHEET_CHECK      } from '../modules/local/samplesheet_check'

/*
=======================================================================================================
RUN MAIN WORKFLOW
=======================================================================================================
*/
workflow INPUT_CHECK {
    take:
        samplesheet
    
    main:
        SAMPLESHEET_CHECK (samplesheet)
        // .csv
        // .splitCsv(heater:true,sep:',')
        // .view { row -> "${row[0]} - ${row[1]} - ${row[2]}" }
        //.map {create_mapped_channel(it)}
        //.set { reads }
    
    emit:
        // reads               // channel: [ val(meta), [ reads ] ]
        samplesheet = SAMPLESHEET_CHECK.out.csv
}

// def create_mapped_channel(LinkedHashMap row){
//     def meta=[:]
//         meta.id             = row.sample
//         meta.single_end     = row.single_end.toBoolean()

//     def array = []
//         if (!file(row.fastq_1).exists()){
//             exit1, "ERROR: missing FASTQ file R1\n${row.fastq_1}"
//         }
//         if (meta.single_end){
//             array = [meta,[file(row.fastq_1)]]
//         } else {
//             if(!file(row.fastq_2).exists()){
//                 exit1, "ERROR: missing FASTQ file R2\n${row.fastq_2}"
//             }
//             array=[meta,[file(row.fastq_1),file(row.fastq_2)]]
//         }
//     return array
// }