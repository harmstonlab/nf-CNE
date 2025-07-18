#!/usr/bin/env nextflow
process resolveReference {
    tag "${reference}"
    label "utils"

    memory '1GB'
    cpus 1
    
    input:
    path reference

    
    output:
    tuple val("${reference.name.takeWhile { it != '.' }}"), path("*.fa"), optional: true, emit: fa
    tuple val("${reference.name.takeWhile { it != '.' }}"), path("*.fa.fai"), optional: true, emit: fai
    tuple val("${reference.name.takeWhile { it != '.' }}"), path("*.2bit"), optional: true, emit: bit
    path "ERROR.txt", optional: true, emit: error_file

    
    script:
    """
    set -eo pipefail
    shopt -s nullglob

    ref_input="\$(realpath "$reference")"

    # Set up resolved paths
    resolved_fa=""
    resolved_fai=""
    resolved_2bit=""

    # Case 1: input is a .fa file
    if [[ "\$ref_input" == *.fa ]]; then
        dir=\$(dirname "\$ref_input")
        base=\$(basename "\$ref_input" .fa)

        resolved_fa="\$ref_input"
        resolved_fai="\$dir/\$base.fa.fai"
        resolved_2bit="\$dir/\$base.2bit"

    # Case 2: input is a directory
    elif [[ -d "\$ref_input" ]]; then
        fa_files=( "\$ref_input"/*.fa )
        fai_files=( "\$ref_input"/*.fa.fai )
        bit_files=( "\$ref_input"/*.2bit )

        if [ \${#fa_files[@]} -ne 1 ]; then
            echo "Expected exactly one .fa file in directory '\$ref_input', found \${#fa_files[@]}" > ERROR.txt
            exit 0
        fi
        if [ \${#fai_files[@]} -ne 1 ]; then
            echo "Expected exactly one .fa.fai file in directory '\$ref_input', found \${#fai_files[@]}" > ERROR.txt
            exit 0
        fi
        if [ \${#bit_files[@]} -ne 1 ]; then
            echo "Expected exactly one .2bit file in directory '\$ref_input', found \${#bit_files[@]}" > ERROR.txt
            exit 0
        fi

        resolved_fa="\${fa_files[0]}"
        resolved_fai="\${fai_files[0]}"
        resolved_2bit="\${bit_files[0]}"

    else
        echo "Input must be a .fa file or a directory: got '\$ref_input'" > ERROR.txt
        exit 0
    fi

    # Final sanity check
    for f in "\$resolved_fa" "\$resolved_fai" "\$resolved_2bit"; do
        if [ ! -f "\$f" ]; then
            echo "Missing required file: \$f" > ERROR.txt
            exit 0
        fi
    done

    rm -f ERROR.txt

    ln -s -f "\$resolved_fa" .
    ln -s -f "\$resolved_fai" .
    ln -s -f "\$resolved_2bit" .
    """
}



workflow{
    if (!params.reference){
        error "Please reference genome path"
    }

    Channel
        .value(params.reference)
        .set { reference_ch }
    
    result=resolveReference(reference_ch)
    result.fa.view()
    result.fai.view()
}