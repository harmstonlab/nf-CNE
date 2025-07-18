#!/usr/bin/env nextflow

include { filterBLAT               } from './filters/blat.nf'
include { combineFilteredBed ; combineFinalBed } from './filters/combine_filtered_beds.nf'
include { identifyConservedRegions } from './pairwise_identification/pairwise_identification.nf'
include { downloadReference        } from './utils/download_ref.nf'
include { gunzipFile               } from './utils/gunzip_file.nf'
include { queryUCSC                } from './utils/list_ucsc_bed.nf'
include { downloadUCSC ; combineUCSCbeds } from './utils/download_UCSC_track.nf'
include { resolveReference         } from './utils/resolve_ref.nf'

workflow {
    if (params.help) {
        helpMessage()
        exit(0)
    }
    // array of allowed filters/tools, any new filters should be added IN UPPER CASE
    allowed_tool_filters = ["BLAT"]

    valid_input_headers = ["species", "reference", "chain", "UCSC_filters", "Tool_Filters", "BLAT_hits_threshold", "BLAT_identity_threshold"]

    input_tsv = file(params.input_tsv ?: 'input.tsv')
    def lines = input_tsv.text.readLines()
    
    def headers = lines[0].split('\t')

    def invalid_headers = headers.findAll { !(it in valid_input_headers) }

    if (invalid_headers) {
        error("${params.input_tsv} cotanins invalid headers: ${invalid_headers.join(', ')}.\nValid headers are: ${valid_input_headers.join(', ')}")
    }

    def input_info = lines[1..-1].collect { line ->
        def fields = line.split('\t', -1).collect { it.trim() == '' ? null : it.trim() }
        def rowMap = [:]
        headers.eachWithIndex { key, i ->
            rowMap[key] = fields[i]
        }
        return rowMap
    }
    input_info.eachWithIndex { entry, idx ->
        entry.order = idx
    }

    if (params.identity <= 0 || params.columns <= 0) {
        error("Both identity and columns must be greater than 0.")
    }

    if (params.identity > params.columns) {
        error("Identity must be less than or equal to columns")
    }

    if (input_info.size() < 2) {
        error("At least two rows must be provided in ${params.input_tsv}. Provide file only contained ${input_info.size()}")
    }

    if (input_info.size() > 2) {
        error("${params.input_tsv} contains ${input_info.size()} species. MSA not yet supported")
    }



    input_info.each { entry ->
        def tools_raw = entry['Tool_Filters'] ?: ''
        def tools_list = tools_raw.split(',').collect { it.trim() }

        def invalid_tools = tools_list.findAll { !(it in allowed_tool_filters) }

        if (invalid_tools && invalid_tools != [] && tools_raw != "") {
            error("Invalid tools found for ${entry.species}: ${invalid_tools.join(', ')}\nValid tool filters are: ${allowed_tool_filters.join(', ')}")
        }
    }

    input_info.each { entry ->
        // parse the comma-list of tool names
        def tools = (entry['Tool_Filters'] ?: '')
            .split(',')
            .collect { it.trim() }
            .findAll { it }

        // only if BLAT was requested
        if (tools.contains('BLAT')) {
            // hits threshold: must be an integer > 0
            def hitsStr = entry['BLAT_hits_threshold']
            def int hits
            try {
                hits = hitsStr as Integer
            }
            catch (Exception _e) {
                error("Invalid BLAT_hits_threshold for species '${entry.species}': '${hitsStr}' is not an integer.")
            }
            if (hits <= 0) {
                error("Invalid BLAT_hits_threshold for species '${entry.species}': ${hits} (must be > 0).")
            }

            // identity threshold: must be a number in (0,1)
            def idStr = entry['BLAT_identity_threshold']
            def double ident
            try {
                ident = idStr as Float
            }
            catch (Exception _e) {
                error("Invalid BLAT_identity_threshold for species '${entry.species}': '${idStr}' is not a number.")
            }
            if (ident <= 0.0 || ident > 1.0) {
                error("Invalid BLAT_identity_threshold for species '${entry.species}': ${ident} (must be >0 and ≤1).")
            }
        }
    }


    refs_to_resolve_ch = Channel.empty()
    refs_to_get_ch = Channel.empty()
    refs_fa_ch = Channel.empty()
    refs_fai_ch = Channel.empty()
    refs_2bit_ch = Channel.empty()
    filtered_beds_ch = Channel.empty()
    input_info.each { entry ->
        def ref = entry.reference

        if (!ref) {
            // Case 1: Empty or null reference — download it
            println("Reference missing for ${entry.species}, will download")
            refs_to_get_ch = refs_to_get_ch.concat(Channel.value(entry.species))
        }
        else if (ref.contains('/')) {
            // Case 2: Looks like a path — check if it exists
            def f = file(ref)
            if (!f.exists()) {
                error("${ref} does not exist and appears to be a path.\nPlease supply a valid path to the reference file for ${entry.species} or leave blank and nf-CNE will download it.")
            }
            else {
                println("Using local reference file for ${entry.species}: ${ref}")
                refs_to_resolve_ch = refs_to_resolve_ch.concat(Channel.value(f))
            }
        }
        else {
            // Case 3: Not a path — treat as species name to download
            println("Reference given as species name for ${entry.species}: ${ref}")
            refs_to_get_ch = refs_to_get_ch.concat(Channel.value(ref))
        }
    }


    downloaded_refs = downloadReference(refs_to_get_ch)
    resolved_refs = resolveReference(refs_to_resolve_ch)
    resolved_refs.error_file
        .filter { it.exists() }
        .map { file -> error("Reference resolution failed: ${file.text.trim()}") }

    refs_fa_ch = downloaded_refs.fa.concat(resolved_refs.fa)
    refs_fai_ch = downloaded_refs.fai.concat(resolved_refs.fai)
    refs_2bit_ch = downloaded_refs.bit.concat(resolved_refs.bit)


    refs_fai_ch
        .map { species, fai_path -> [(species): fai_path] }
        .collect()
        .map { mapList ->
            mapList.inject([:]) { acc, m -> acc + m }
        }
        .set { fai_map_ch }

    // Check and download any UCSC tracks
    Channel.fromList(input_info)
        .map { entry -> tuple(entry.species, file('./utils/list_UCSC_bed.py')) }
        .distinct()
        .set { query_ucsc_ch }

    UCSC_result = queryUCSC(query_ucsc_ch)
    UCSC_result.result
        .map { species, stdout ->
            def tracks = stdout
                .readLines()
                .collect { it.trim() }
                .findAll { it }
            tuple(species, tracks)
        }
        .subscribe { species, validTracks ->
            def entry = input_info.find { it.species == species }
            // pull the raw filters (might be null)
            def raw = entry.UCSC_filters
            if (!raw) {
                // nothing to check, skip
                return null
            }

            // clean up the valid list
            def cleanValid = validTracks.collect { it.trim() }

            // strip quotes, split, trim
            def requested = raw
                .replaceAll('"', '')
                .replaceAll("'", "")
                .split(',')
                .collect { it.trim() }
                .findAll { it }

            // now only if they actually asked for something do we check
            if (requested) {
                def missing = requested.findAll { !cleanValid.contains(it) }
                if (missing) {
                    error("• ${cleanValid.join('\n• ')}\nERROR: The following UCSC filters were not found for species '${species}':\n• ${missing.join('\n  • ')}\nValid filters are above")
                }
            }
        }



    // Turn List<Map> into a channel of (species,track) tuples
    // one tuple per requested UCSC track
    download_requests_ch = Channel.fromList(input_info)
        .flatMap { entry ->
            def raw = entry.UCSC_filters
            if (!raw) {
                // no UCSC filters → skip this entry
                return []
            }
            raw
                .replaceAll(/["']/, '')
                .split(',')
                .collect { it.trim() }
                .findAll { it }
                .collect { track ->
                    tuple(entry.species, track)
                }
        }

    // Cross with the Python script path
    ready_to_download_ch = UCSC_result.done
        .cross(download_requests_ch)
        .map { doneTuple, reqTuple ->
            // unpack
            def (spDone, doneVal) = doneTuple
            def (spReq, track) = reqTuple

            // only keep when species match
            if (spDone == spReq) {
                tuple(spReq, doneVal, track)
            }
            else {
                null
            }
        }
        .filter { it != null }

    download_inputs_ch = ready_to_download_ch.map { species, doneVal, track ->
        tuple(
            species,
            doneVal,
            track,
            file('./utils/get_UCSC_bed.py'),
        )
    }

    UCSC_all_files = downloadUCSC(download_inputs_ch)

    // Group per-species, sort by filename, then combine into one bed file
    UCSC_all_files
        .groupTuple(by: 0)
        .map { species, files ->
            def sorted = files.sort { it.name }
            tuple(species, sorted)
        }
        .set { sorted_ucsc_files_ch }
    final_bed_ch = combineUCSCbeds(sorted_ucsc_files_ch)

    Channel.from(input_info)
        .map { entry ->

            tuple(entry.species.toString(), entry)
        }
        .set { info_ch }

    final_bed_ch
        .map { species, bedFile -> tuple(species.toString(), bedFile) }
        .set { bed_ch }

    info_ch
        .join(bed_ch)
        .map { _species, entry, bedFile ->
            entry + [UCSC_filters: bedFile.toString()]
        }
        .collect()
        .set { input_info_UCSC_list_ch }

    input_info_UCSC = input_info_UCSC_list_ch
        .flatMap { it }
        .map { entry -> tuple(entry.species, entry) }


    if (input_info.size() == 2) {
        // Check that we have two files that exist. If both don't exist pairwise alignment will run (TODO)
        def (entry1, entry2) = input_info

        def f1 = entry1.chain ? file(entry1.chain) : null
        def f2 = entry2.chain ? file(entry2.chain) : null

        def exists1 = f1?.exists()
        def exists2 = f2?.exists()

        if (!exists1 && !exists2) {
            error("Generating pairwise alignments from scratch is not yet supported — please provide chain files.")
        }
        if ((exists1 && !exists2) || (!exists1 && exists2)) {
            error("One chain file is missing. Found: ${entry1.species}: ${f1?.name}, ${entry2.species}: ${f2?.name}. Check input TSV.")
        }


        def gzipped = input_info.findAll { it.chain.endsWith('.gz') }
        def plain = input_info.findAll { !it.chain.endsWith('.gz') }

        Channel.fromList(gzipped)
            .map { entry -> tuple(entry, file(entry.chain)) }
            .set { gzipped_ch }
        gunzipFile(gzipped_ch)
            .map { entry, unzipped_path ->
                entry.chain = unzipped_path.toString()
                return entry
            }
            .concat(Channel.fromList(plain))
            .collect()
            .set { gunzipped_sample_info_ch }

        chains_ch = gunzipped_sample_info_ch
            .flatMap { it }
            .map { entry -> tuple(entry.species, entry) }

        ucsc_ch = input_info_UCSC.map { species, entry -> tuple(species, entry) }


        final_ch = chains_ch
            .join(ucsc_ch)
            .map { _species, gzEntry, ucscEntry ->
                // merge but preserve order
                def merged = gzEntry + [
                    UCSC_filters: ucscEntry.UCSC_filters,
                    order: gzEntry.order,
                ]
                return merged
            }
        pair_maps_ch = final_ch
            .collect()
            .map { list ->
                list.sort { it.order }
            }

        def dummy_filter_path = file("${workflow.workDir}/empty_placeholder.bed")
        dummy_filter_path.text = ""


        identify_inputs_ch = pair_maps_ch
            .combine(fai_map_ch)
            .map { records ->
                def (m1, m2) = records

                def filter1 = m1.UCSC_filters ?: dummy_filter_path
                def filter2 = m2.UCSC_filters ?: dummy_filter_path
                def fai_1 = fai_map_ch[m1.species]
                def fai_2 = fai_map_ch[m2.species]
                tuple(
                    file(fai_1.value),
                    file(fai_2.value),
                    file(m1.chain),
                    file(m2.chain),
                    file(filter1),
                    file(filter2),
                    params.identity,
                    params.columns,
                    file('./pairwise_identification/identify_pairwise.py'),
                )
            }


        (bed_1, bed_2) = identifyConservedRegions(identify_inputs_ch)

        bed_files_ch = bed_1.concat(bed_2)
    }

    input_info_ch = Channel.fromList(input_info)

    blat_params_ch = input_info_ch
        .filter { entry ->
            entry['Tool_Filters']
                .replaceAll('"', '')
                .replaceAll("'", "")
                .split(',')
                .collect { it.trim() }
                .contains('BLAT')
        }
        .map { entry ->
            def cutoff = (entry['BLAT_hits_threshold'] as Integer)
            def identity = (entry['BLAT_identity_threshold'] as Float)

            tuple(entry.species, cutoff, identity)
        }

    ref_bed_ch = refs_2bit_ch
        .join(bed_files_ch)
        .map { species, ref2bit, bed ->
            tuple(species, ref2bit, bed)
        }

    to_filter_ch = ref_bed_ch
        .join(blat_params_ch)
        .map { species, ref2bit, bed, cutoff, identity ->
            tuple(species, ref2bit, bed, cutoff, identity)
        }


    filter_inputs_ch = to_filter_ch.map { species, ref2bit, bed, cutoff, identity ->
        tuple(
            species,
            ref2bit,
            bed,
            cutoff,
            identity,
            file('./filters/BLAT_filter.py'),
        )
    }

    BLAT_filtered_beds_ch = filterBLAT(filter_inputs_ch)
    filtered_files_ch = filtered_beds_ch.concat(BLAT_filtered_beds_ch)



    filtered_group_size_ch = filtered_files_ch
        .groupTuple(by: 0)
        .map { species, filtered_files ->
            tuple(species, filtered_files.size())
        }

    filter_input_ch = bed_files_ch.join(filtered_files_ch.groupTuple(by: 0).join(filtered_group_size_ch))
    filtered_beds = combineFilteredBed(filter_input_ch)

    final_bed_ch = filtered_beds.collect()

    Channel.fromPath('./filters/filter_missing_CNEs.py')
        .set { combine_final_filters_script_ch }

    combineFinalBed(input_info.size(), final_bed_ch, combine_final_filters_script_ch)
}

def helpMessage() {
    println(
        """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --columns 50 --identity 50 --input_tsv input.tsv

        Arguments:
        --columns                       The size of the sliding window to be used when searching for conserved elements (Default: 50)
        --identity                      The minimum number of matching bases within the columns window (Default: 50)
        --input_tsv                     The path to a TSV file containing species specific options (see README for details - Default: input.tsv)
        --outdir                        The location of where to store the outputs (Default: ./results)
        --memory                        The amount of memory to allocate to the BLAT process (Default: 8GB)
        --help                          This usage statement
        """
    )
}