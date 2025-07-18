#!/usr/bin/env nextflow

process queryUCSC {
    tag "${species}"
    label "utils"

    memory '1GB'
    cpus 1

    input:
    tuple val(species), path(list_UCSC_script)

    output:
    tuple val(species), stdout, emit: result
    tuple val(species), val("done"), emit: done //dummy value, consumed by filterUCSC in main.nf so that it doesn't start before check complete
    

    script:
    """
    python ${list_UCSC_script} -s ${species}
    """
}

workflow{
    if (!params.species){
        error "Please provide a species to search for UCSC tracks"
    }

    Channel
        .value(params.species)
        .set { species_ch }

    Channel
        .fromPath('./list_UCSC_bed.py')
        .set { script_ch }
    query_ch=species_ch.combine(script_ch)
    queryUCSC(query_ch)
    queryUCSC.out.map { _species, stdout -> stdout } .view()
    
}