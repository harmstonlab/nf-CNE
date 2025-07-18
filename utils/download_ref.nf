#!/usr/bin/env nextflow

process downloadReference {
    tag "${species}"
    label "utils"

    storeDir "${params.outdir}/reference/${species}"
    memory '1GB'
    cpus 1

    input:
    val species

    output:
    tuple val(species), path("${species}.fa"), emit: fa
    tuple val(species), path("${species}.2bit"), emit: bit
    tuple val(species), path("${species}.fa.fai"), emit: fai

    script:
    
    """
    wget https://hgdownload.soe.ucsc.edu/gbdb/${species}/${species}.2bit
    twoBitToFa ${species}.2bit ${species}.fa
    samtools faidx ${species}.fa
    """
}

workflow{
    if (!params.species){
        error "Please provide a species to download the reference genome for"
    }

    Channel
        .value(params.species)
        .set { species_ch }
    
    downloadReference(species_ch)
}