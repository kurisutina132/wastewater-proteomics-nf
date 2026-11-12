#!/usr/bin/env nextflow
/*
 * wastewater-proteomics-nf
 * Lightweight pipeline for shotgun proteomics of wastewater:
 * PRIDE -> mzML -> Sage (search + FDR + LFQ) -> peptide taxonomy ->
 * protein quant -> biomarker / immunoglobulin panel -> stats -> HTML report
 */
nextflow.enable.dsl = 2

def helpMessage() {
    log.info """
    wastewater-proteomics-nf

    Input (choose one):
      --pride PXD059455        download a subset of RAW files from PRIDE
      --samplesheet file.csv   columns: sample,file,group  (file = .raw/.mzML, local path or URL)
      -profile test            synthetic data, no download

    Main options:
      --max_files      number of RAW files to take from PRIDE   [${params.max_files}]
      --group_regex    regex (1 capture group) on file names -> group
      --db_set         minimal | full                             [${params.db_set}]
      --fasta          use your own FASTA instead of building one
      --fdr            PSM/peptide FDR                            [${params.fdr}]
      --outdir         results folder                             [${params.outdir}]
    """.stripIndent()
}

/* ------------------------------------------------------------------ */
process MAKE_TEST_DATA {
    label 'python'
    output:
    path "test_data/*.mzML",              emit: mzml
    path "test_data/test_samplesheet.csv", emit: samplesheet
    path "test_data/test_db.fasta",        emit: fasta
    path "test_data/test_taxmap.tsv",      emit: taxmap
    script:
    """
    make_test_data.py --outdir test_data
    """
}

process FETCH_PRIDE {
    label 'python'
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'
    output:
    path "samplesheet.csv", emit: samplesheet
    path "pride_files.tsv"
    script:
    def inc = params.include_regex ? "--include_regex '${params.include_regex}'" : ''
    def grp = params.group_regex   ? "--group_regex '${params.group_regex}'"     : ''
    """
    fetch_pride.py --accession ${params.pride} --max_files ${params.max_files} \\
        --max_file_mb ${params.max_file_mb} ${inc} ${grp}
    """
}

process RAW_TO_MZML {
    tag "$sample"
    label 'convert'
    input:
    tuple val(sample), val(group), path(raw)
    output:
    tuple val(sample), val(group), path("${sample}.mzML")
    script:
    """
    ${params.trfp_bin} -i ${raw} -o . -f 2
    mv ${raw.baseName}.mzML ${sample}.mzML 2>/dev/null || true
    """
}

process BUILD_DB {
    label 'python'
    publishDir "${params.outdir}/database", mode: 'copy'
    input:
    path proteomes
    output:
    path "db.fasta",        emit: fasta
    path "taxmap.tsv",      emit: taxmap
    path "db_manifest.tsv", emit: manifest
    script:
    """
    build_db.py --proteomes ${proteomes} --set ${params.db_set} \\
        --contaminants ${params.contaminants}
    """
}

process SAGE_SEARCH {
    label 'sage'
    publishDir "${params.outdir}/sage", mode: 'copy'
    input:
    path fasta
    path mzmls
    output:
    path "sage_out/results.sage.tsv", emit: psms
    path "sage_out/lfq.tsv",          emit: lfq
    path "sage_out/results.json"
    path "sage_config.json"
    script:
    """
    make_sage_config.py --fasta ${fasta} \\
        --missed_cleavages ${params.missed_cleavages} \\
        --precursor_ppm ${params.precursor_ppm} --fragment_ppm ${params.fragment_ppm} \\
        --carbamidomethyl ${params.carbamidomethyl} --oxidation ${params.oxidation} \\
        --min_matched_peaks ${params.min_matched_peaks}
    ${params.sage_bin} sage_config.json -o sage_out ${mzmls}
    """
}

process ANALYZE {
    label 'python'
    publishDir "${params.outdir}", mode: 'copy'
    input:
    path psms
    path lfq
    path fasta
    path taxmap
    path samplesheet
    path panel
    output:
    path "results/*"
    path "results/report.html", emit: report
    script:
    def tm = taxmap.name != 'NO_TAXMAP' ? "--taxmap ${taxmap}" : ''
    """
    analyze.py --psms ${psms} --lfq ${lfq} --fasta ${fasta} ${tm} \\
        --samplesheet ${samplesheet} --panel ${panel} \\
        --fdr ${params.fdr} --min_peptides ${params.min_peptides} \\
        --title "${params.report_title}" --outdir results
    """
}

/* ------------------------------------------------------------------ */
workflow {
    if (params.help) { helpMessage(); return }

    // ---------- inputs -> (sample, group, file) ----------------------
    if (params.test_data) {
        MAKE_TEST_DATA()
        by_name = MAKE_TEST_DATA.out.mzml.flatten().map { f -> [f.name, f] }
        ch_runs = MAKE_TEST_DATA.out.samplesheet
            .splitCsv(header: true)
            .map { r -> [r.file, r.sample, r.group] }
            .join(by_name)
            .map { fname, s, g, f -> [s, g, f] }
        ch_fasta  = MAKE_TEST_DATA.out.fasta
        ch_taxmap = MAKE_TEST_DATA.out.taxmap
    } else {
        if (params.pride) {
            FETCH_PRIDE()
            ch_sheet = FETCH_PRIDE.out.samplesheet
        } else if (params.samplesheet) {
            ch_sheet = Channel.fromPath(params.samplesheet, checkIfExists: true)
        } else {
            helpMessage(); error "Provide --pride, --samplesheet or -profile test"
        }
        ch_runs = ch_sheet.splitCsv(header: true)
            .map { r -> [r.sample, r.group ?: 'all', file(r.file)] }

        if (params.fasta) {
            ch_fasta  = Channel.fromPath(params.fasta, checkIfExists: true)
            ch_taxmap = params.taxmap ? Channel.fromPath(params.taxmap, checkIfExists: true)
                                      : Channel.value(file("${projectDir}/assets/NO_TAXMAP"))
        } else {
            BUILD_DB(Channel.fromPath(params.proteomes, checkIfExists: true))
            ch_fasta  = BUILD_DB.out.fasta
            ch_taxmap = BUILD_DB.out.taxmap
        }
    }

    // ---------- RAW -> mzML ------------------------------------------
    branched = ch_runs.branch { s, g, f ->
        raw:  f.name.toLowerCase().endsWith('.raw')
        mzml: true
    }
    ch_mzml = RAW_TO_MZML(branched.raw).mix(branched.mzml)

    // ---------- samplesheet for the analysis step ---------------------
    ch_analysis_sheet = ch_mzml
        .map { s, g, f -> "${s},${g},${f.name}" }
        .collectFile(name: 'analysis_samplesheet.csv', newLine: true, sort: true,
                     seed: 'sample,group,mzml', storeDir: "${params.outdir}/pipeline_info")

    // ---------- search + analysis -------------------------------------
    SAGE_SEARCH(ch_fasta, ch_mzml.map { s, g, f -> f }.collect())
    ANALYZE(SAGE_SEARCH.out.psms, SAGE_SEARCH.out.lfq, ch_fasta, ch_taxmap,
            ch_analysis_sheet, file(params.panel))
}

workflow.onComplete {
    log.info(workflow.success ? "\nDone. Report: ${params.outdir}/results/report.html\n"
                              : "\nPipeline failed: ${workflow.errorMessage}\n")
}
