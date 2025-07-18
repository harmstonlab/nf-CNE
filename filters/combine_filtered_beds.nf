#!/usr/bin/env nextflow

process combineFilteredBed {
    tag "${count} ${count == 1 ? 'file' : 'files'}"
    label "filters"
    memory '2GB'
    cpus 1
    publishDir "${params.outdir}/COMBINED_filtered", mode: 'copy'

    input:
    tuple val(species), path(unfiltered), path(filtered), val(count)

    output:
    path "${unfiltered.name.replace('.bed', "_combined_filtered.bed")}", emit: combined_filter_bed
    

    script:
    """
    num_files=\$(ls *filtered.bed | wc -l)
    empty_files=\$(wc -l ./*.bed | grep "^\\s*0 ./" | grep -v "total" | wc -l)
    # If there are empty files, create an empty output file and skip the intersection
    if [[ \$empty_files -gt 0 ]]; then
        touch ${unfiltered.name.replace('.bed', "_combined_filtered.bed")}
    else
        # Perform the intersection only if there are no empty files
        bedtools intersect -a ${unfiltered} -b *filtered.bed -c -f 1 -F 1 | grep "\${num_files}\\\$" > ${unfiltered.name.replace('.bed', "_combined_filtered.bed")}
    fi    

    """
}



process combineFinalBed {
    tag "${count} ${count == 1 ? 'file' : 'files'}"
    label "filters"
    memory '2GB'
    cpus 1
    publishDir "${params.outdir}/FINAL_FILTERED", mode: 'copy'

    input:

    val count
    path filtered
    path script

    output:
    path "*FINAL_filtered.bed"

    script:
    """
    python3 ${script}

    """
}