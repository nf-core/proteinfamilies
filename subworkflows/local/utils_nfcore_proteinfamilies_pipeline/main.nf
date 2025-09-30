//
// Subworkflow with functionality specific to the nf-core/proteinfamilies pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        "",
        false,
        false,
        false,
        "",
        "",
        ""
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //

    ch_samplesheet = Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def quality_check_text = [
        "Amino acid sequence statistics were generated with SeqFu (Telatin et al. 2021).",
        params.skip_preprocessing ? "" : "Input sequences were preprocessed with SeqKit (gap trimming, length filtering, validation, duplicate removal) (Shen et al. 2024)."
    ].join(' ').trim()

    def clustering_text = "Amino acid sequence clustering was performed with MMseqs2 (Mirdita et al. 2021)."

    def alignment_text = [
        "Multiple Sequence Alignment (MSA) was performed with ",
        params.alignment_tool == 'famsa' ? "FAMSA (Deorowicz et al. 2026)." : "",
        params.alignment_tool == 'mafft' ? "mafft (Katoh et al. 2013)." : ""
    ].join(' ').trim()

    def clipping_text = "MSAs were gap-trimmed with ClipKIT (Steenwyk et al. 2020)."

    def model_text = "Family Hidden Markov Models (HMMs) were built with hmmer (Potter et al. 2018)."

    def postprocessing_text = "Run statistics were reported using MultiQC (Ewels et al. 2016)."

    def citation_text = [
        quality_check_text,
        clustering_text,
        alignment_text,
        params.trim_msa ? clipping_text : "",
        model_text,
        postprocessing_text
    ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def quality_check_text = [
        '<li>Telatin, A., Fariselli, P., & Birolo, G. (2021). SeqFu: a suite of utilities for the robust and reproducible manipulation of sequence files. Bioengineering, 8(5), 59. doi: <a href="https://doi.org/10.3390/bioengineering8050059">10.3390/bioengineering8050059</a></li>',
        params.skip_preprocessing ? '' : '<li>Shen, W., Sipos, B., & Zhao, L. (2024). SeqKit2: A Swiss army knife for sequence and alignment processing. Imeta, 3(3), e191. doi: <a href="https://doi.org/10.1002/imt2.191">10.1002/imt2.191</a></li>'
    ].join(' ').trim()

    def clustering_text = '<li>Mirdita, M., Steinegger, M., Breitwieser, F., Söding, J., & Levy Karin, E. (2021). Fast and sensitive taxonomic assignment to metagenomic contigs. Bioinformatics, 37(18), 3029-3031. doi: <a href="https://doi.org/10.1093/bioinformatics/btab184">10.1093/bioinformatics/btab184</a></li>'

    def alignment_text = [
        params.alignment_tool == 'famsa' ? '<li>Deorowicz, S., Debudaj-Grabysz, A., & Gudyś, A. (2016). FAMSA: Fast and accurate multiple sequence alignment of huge protein families. Scientific reports, 6(1), 33964. doi: <a href="https://doi.org/10.1038/srep33964">10.1038/srep33964</a></li>' : "",
        params.alignment_tool == 'mafft' ? '<li>Katoh, K., & Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. Molecular biology and evolution, 30(4), 772-780. doi: <a href="https://doi.org/10.1093/molbev/mst010">10.1093/molbev/mst010</a></li>' : "",
    ].join(' ').trim()

    def clipping_text = '<li>Steenwyk, J. L., Buida III, T. J., Li, Y., Shen, X. X., & Rokas, A. (2020). ClipKIT: a multiple sequence alignment trimming software for accurate phylogenomic inference. PLoS biology, 18(12), e3001007. doi: <a href="https://doi.org/10.1371/journal.pbio.3001007">10.1371/journal.pbio.3001007</a></li>'

    def model_text = '<li>Eddy, S. R. (2011). Accelerated profile HMM searches. PLoS computational biology, 7(10), e1002195. doi: <a href="https://doi.org/10.1371/journal.pcbi.1002195">10.1371/journal.pcbi.1002195</a></li>'

    def postprocessing_text = '<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19), 3047–3048. doi: <a href="https://doi.org/10.1093/bioinformatics/btw354">10.1093/bioinformatics/btw354</a></li>'

    def reference_text = [
        quality_check_text,
        clustering_text,
        alignment_text,
        params.trim_msa ? clipping_text : "",
        model_text,
        postprocessing_text
    ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

//
// Validate that two folders (HMMs and MSAs for update) contain the same number of files and matching base filenames
//
def validateMatchingFolders(channel1, channel2) {
    // Fetch the contents of the channels
    channel1
        .join(channel2)
        .map { meta, folder1, folder2 ->
            def files1 = folder1.listFiles()
            def files2 = folder2.listFiles()

            // Check if the number of files matches
            if (files1.size() != files2.size()) {
                error("[nf-core/proteinfamilies] ERROR: Folder mismatch: ${folder1} has ${files1.size()} files, but ${folder2} has ${files2.size()} files.")
            }

            // Extract base filenames (without extensions) and sort
            def baseNames1 = files1.collect { it.getSimpleName() }.sort()
            def baseNames2 = files2.collect { it.getSimpleName() }.sort()

            // Check if base filenames match one to one
            if (baseNames1 != baseNames2) {
                error("[nf-core/proteinfamilies] ERROR: Filename mismatch: Expected matching files in ${folder1} and ${folder2}. Base filenames do not match.")
            }
    }
}
