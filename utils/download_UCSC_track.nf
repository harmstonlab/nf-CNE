#!/usr/bin/env nextflow

process downloadUCSC {
    tag "${species}-${track}"
    label "utils"
    storeDir "${params.outdir}/UCSC_bed"

    memory '8GB'
    cpus 1

    input:
    //triggier is a dummy value, emmited by queryUCSC so that it doesn't start before check complete as could lead to crashes if invalid track provided
    tuple val(species), val(trigger), val(track), path(script)

    output:
    tuple val(species), path("${species}_${track}.bed"), emit: bed

    script:
    """
    python ${script} --species ${species} --track ${track}
    """
}

process combineUCSCbeds{
    tag "${species}"
    label "utils"
    cpus 1
    storeDir "${params.outdir}/UCSC_bed"

    input:
    tuple val(species), path(beds)

    output:
    tuple val(species), file(finalName), emit: merged_bed

    script:

    def bedBaseNames = beds.collect { it.baseName }
    def prefix       = bedBaseNames[0] 
    def prefixCommon = prefix.substring(0, prefix.indexOf('_')+1) 
    
    def restList = bedBaseNames.drop(1).collect { it - prefixCommon }   
    finalName = prefix + (restList ? '_' + restList.join('_') : '') + '.bed'
    """
    cat ${beds} \
      | sort -k1,1 -k2,2n \
      | awk -F '\t' 'BEGIN {OFS = FS} { print \$1, \$2, \$3 }' \
      | bedtools merge  > "${finalName}"
    """
}

workflow{
    if (!params.species){
        error "Please provide a species to use when downloading the UCSC track"
    }

    if (!params.track){
        error "Please provide a UCSC track to download"
    }


    Channel
        .value(params.species)
        .set { species_ch }


    tracks = params.track.split(',')
            .findAll { it}
            .collect { it} 
    Channel
            .fromList(tracks)
            .set { track_ch }

    Channel
        .fromPath("./get_UCSC_bed.py")
        .set { ucsc_filter_script_ch }

    Channel
        .value(true)
        .set { trigger_ch }

    files=downloadUCSC(species_ch.combine(trigger_ch).combine(track_ch).combine(ucsc_filter_script_ch))
    
    combineUCSCbeds(files.groupTuple(by:0))
    
}