#!/usr/bin/env nextflow

process gunzipFile {
    tag "${input_file}"
    label "utils"
    memory '1GB'
    cpus 1

    input:
    tuple val(entry), path(input_file)

    output:
    tuple val(entry), path("${input_file.name.replace('.gz', "")}"), emit: out_file

    script:
    """
    
    gunzip -f ${input_file}

    """
}
