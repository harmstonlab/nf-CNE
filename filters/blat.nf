#!/usr/bin/env nextflow

process filterBLAT {
    tag "${species}"
    label "filters"
    storeDir "${params.outdir}/BLAT_filtered"
    memory params.memory

    input:
    tuple val(species), path(reference), path(bed), val(blat_cuttoff), val(blat_identity), val(filter_script)

    output:
    tuple val(species), path("${bed.name.replace('.bed', '_BLAT_filtered.bed')}")
    

    script:
    """
    # pre touch file to create it, as if bed is empty it won't create the input file
    touch input.txt result.psl
    gawk -F' ' '{ print "${reference}:"\$1":"\$2"-"\$3 > "input.txt" }'  ${bed}
    # "or true" is a bit of a hack... if bed is empty, pblat fails to read it and gives a 255 exit code, crashing pipeline
    blat ${reference} input.txt result.psl || true

    python ${filter_script} \
        --bed ${bed} \
        --cutoff ${blat_cuttoff} \
        --cut_identity ${blat_identity} \
        --psl result.psl
    """
}

def resolveReference2Bit(reference) {
    if (reference && reference.endsWith('.fai')){
        reference = file(reference).parent
    }
    // Check if the input is a directory
    if (reference && file(reference).isDirectory()) {
        // Look for .fa files in the directory
        def faFiles = file(reference).listFiles().findAll { it.name.endsWith('.2bit') }
        if (faFiles.size() != 1) {
            error "Expected exactly one .2bit file in directory '${reference}', but found ${faFiles.size()}"
        }
        return faFiles[0].toAbsolutePath()
    }
    // If not a directory, assume it is a specific .2bit file
    else if (reference && reference.endsWith('.2bit')) {
        return file(reference).toAbsolutePath()
    } else {
        error "Invalid input for reference: '${reference}'. Provide either a .2bit file or a directory containing one."
    }
}

workflow {

    if (!(params.threads instanceof Integer) || params.threads <= 0) {
        error "Invalid value for 'threads': ${params.threads}. It must be an integer greater than 0."
    }

    if (!(params.blat_cuttoff instanceof Integer) || params.blat_cuttoff<= 0) {
        error "Invalid value for 'blat_cuttoff': ${params.blat_cuttoff}. It must be an integer greater than 0."
    }

    def reference_Path = resolveReference2Bit(params.reference)

    Channel
        .fromPath(reference_Path)
        .set { reference_ch }

    Channel
        .fromPath(params.bed)
        .set { bed_ch }

    Channel
        .value(params.blat_cuttoff)
        .set { blat_cuttoff_ch }

    Channel
        .value(params.blat_identity)
        .set { blat_identity_ch }

    Channel
        .fromPath('./BLAT_filter.py')
        .set { blat_filter_script_ch }
        
    Channel
        .value("${reference_Path.baseName}")
        .set {species_ch}
    
    filterBLAT(species_ch.combine(reference_ch).combine(bed_ch).combine(blat_cuttoff_ch).combine(blat_identity_ch).combine(blat_filter_script_ch))
}
