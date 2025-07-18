#!/usr/bin/env nextflow

process identifyConservedRegions {
    tag "${input_1.name} - ${input_2.name}"
    storeDir "${params.outdir}/unfiltered"
    memory '4GB'
    cpus 2

    input:
    tuple path(reference_1), path(reference_2), path(input_1), path(input_2), path(filter_1_path), path(filter_2_path), val(identity), val(columns), path(identify_script)

    output:
    tuple val(species_1), path("${output_1_name}"), emit: bed_1
    tuple val(species_2),  path("${output_2_name}"), emit: bed_2

    script:
    species_1="${reference_1.name.replace('.fa.fai','')}"
    species_2="${reference_2.name.replace('.fa.fai','')}"

    input_1_base = input_1.name.replace('.net.axt', '')
    input_2_base = input_2.name.replace('.net.axt', '')

    def dummy_name = 'empty_placeholder.bed'
    filter_1_suffix = filter_1_path.getName() != dummy_name ? filter_1_path.name : ''
    filter_2_suffix = filter_2_path.getName() != dummy_name ? filter_2_path.name : ''

    output_1_suffix = filter_1_suffix? "_FILTERED_${filter_1_suffix}" : '.bed'
    output_2_suffix = filter_2_suffix? "_FILTERED_${filter_2_suffix}" : '.bed'

    // Construct output file names
    output_1_name = "${input_1_base}_${identity}I_${columns}col${output_1_suffix}"
    output_2_name = "${input_2_base}_${identity}I_${columns}col${output_2_suffix}"

    filter_1_cli_arg = filter_1_path.getName() != dummy_name ? "--filter_1 ${filter_1_path}" : ''
    filter_2_cli_arg = filter_2_path.getName() != dummy_name ? "--filter_2 ${filter_2_path}" : ''


    """

    python ${identify_script} \
        --identity ${identity} \
        --columns ${columns} \
        --reference_1 ${reference_1} \
        --reference_2 ${reference_2} \
        --input_1 ${input_1} \
        --input_2 ${input_2} \
        ${filter_1_cli_arg} \
        ${filter_2_cli_arg} 
    """
}



def resolveReferenceFai(reference) {
    // Check if the input is a directory
    if (reference && reference.endsWith('.fa')){
        reference = file(reference).parent
    }
    if (reference && file(reference).isDirectory()) {
        // Look for .fa files in the directory
        def faiFiles = file(reference).listFiles().findAll { it.name.endsWith('.fai') }
        if (faiFiles.size() != 1) {
            error "Expected exactly one .fai file in directory '${reference}', but found ${faiFiles.size()}"
        }
        return faiFiles[0].toAbsolutePath()
    }
    // If not a directory, assume it is a specific .fa file
    else if (reference && reference.endsWith('.fai')) {
        return file(reference).toAbsolutePath()
    } else {
        error "Invalid input for reference: '${reference}'. Provide either a .fai file or a directory containing one."
    }
}

workflow {
    if (params.identity <= 0 || params.columns <= 0) {
        error "Both identity and columns must be greater than 0."
    }

    if (params.identity > params.columns) {
        error "Identity must be less than or equal to columns"
    }

    if ((params.input_1 && !params.input_2) || (params.input_2 && !params.input_1) || (!params.input_1 || !params.input_2)){
        error "Both input_1 and input_2 must be provided"
    }


    def reference_1_fai_Path = resolveReferenceFai(params.reference_1)
    def reference_2_fai_Path = resolveReferenceFai(params.reference_2)

    Channel
        .fromPath(reference_1_fai_Path)
        .set { reference_1_fai_ch }

    Channel
        .fromPath(reference_2_fai_Path)
        .set { reference_2_fai_ch }

    Channel
        .fromPath(params.input_1)
        .set { input_1_ch }

    Channel
        .fromPath(params.input_2)
        .set { input_2_ch }


    (params.filter_1 ? Channel.fromPath(params.filter_1) : Channel.empty())
      .map { file -> file.toString() }
      .ifEmpty([])
      .set { filter_1_ch }

    (params.filter_2 ? Channel.fromPath(params.filter_2) : Channel.empty())
      .map { file -> file.toString() }
      .ifEmpty([])
      .set { filter_2_ch }

    Channel
        .value(params.identity)
        .set { identity_ch }

    Channel
        .value(params.columns)
        .set { columns_ch }

    Channel
        .fromPath('./identify_pairwise.py')
        .set { script_ch }


    identifyConservedRegions(
        (reference_1_fai_ch.combine(reference_2_fai_ch).combine(input_1_ch).combine(input_2_ch).combine(filter_1_ch).combine(filter_2_ch).combine(identity_ch).combine(columns_ch).combine(script_ch))
    )

}