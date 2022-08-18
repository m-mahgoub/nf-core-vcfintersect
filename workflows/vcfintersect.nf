/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowVcfintersect.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param && params.multiple_vcf_directory) { Channel.fromPath("${params.input}/*.vcf", checkIfExists: true) } }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
// include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { BCFTOOLS_MERGE                   } from '../modules/local/bcftools_merge.nf'
include { BCFTOOLS_VIEW                  } from '../modules/local/bcftools_view.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                          } from '../modules/nf-core/modules/multiqc/main'
include { TABIX_BGZIPTABIX                 } from '../modules/nf-core/modules/tabix/bgziptabix/main.nf'
include { TABIX_TABIX                      } from '../modules/nf-core/modules/tabix/tabix/main.nf'
include { CUSTOM_DUMPSOFTWAREVERSIONS      } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary


workflow VCFINTERSECT {

    ch_versions = Channel.empty()

    // Define VCF file meta including meta info
    Channel.fromPath("${params.input}/*.vcf")
    .map { 
        it -> 
        def meta = [:]
        meta.id = it.getBaseName()
        [ meta, it ] }
    .set { ch_vcf_from_directory }

    // Channel.fromPath("${params.input}/*.vcf")
    // .map { it -> [ it.getBaseName(), it ] }
    // .set { ch_vcf_from_directory }

    //
    // MODULE: Index VCF files (and compress if uncompressed)
    //
    if (params.compressed_vcf) {
        TABIX_TABIX (
            ch_vcf_from_directory
        )
        ch_compressed_indexed_vcf = TABIX_TABIX.out
        ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())
    }
    else if (!params.compressed_vcf) {
        TABIX_BGZIPTABIX (
            ch_vcf_from_directory
        )
        ch_compressed_indexed_vcf = TABIX_BGZIPTABIX.out
        ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())   
    }

    //
    // MODULE: Merge VCF files (if multiple files provided)
    //
    if (params.multiple_vcf_directory) {
        // Collect all .gz files in one channel
        ch_compressed_indexed_vcf.gz_tbi
        .collect( { it -> it[1] } )
        .set { ch_gz_files }

        // Collect all .tbi files in one channel
        ch_compressed_indexed_vcf.gz_tbi
        .collect( { it -> it[2] } )
        .set { ch_tbi_files }

        BCFTOOLS_MERGE (
            ch_gz_files, 
            ch_tbi_files
        )
        ch_merged_vcf = BCFTOOLS_MERGE.out.merged_variants
        ch_versions = ch_versions.mix(BCFTOOLS_MERGE.out.versions)
    }

    //
    // MODULE: Filter merged VCF file
    BCFTOOLS_VIEW (
        ch_merged_vcf
    )
    ch_filtered_vcf = BCFTOOLS_VIEW.out.vcf
    ch_versions = ch_versions.mix(BCFTOOLS_VIEW.out.versions) 

    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowVcfintersect.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)


    // Emit for testing purpose
    // emit: Channel.empty()
    // emit: ch_compressed_indexed_vcf.gz_tbi
    //  emit: ch_gz_files
    // emit : ch_merged_vcf.merged_variants
    emit: ch_filtered_vcf
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
