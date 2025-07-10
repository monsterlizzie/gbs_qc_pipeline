nextflow.enable.dsl=2
// import modules for process
include { FILE_VALIDATION; PREPROCESS; READ_QC } from "$projectDir/modules/preprocess"
include { ASSEMBLY_UNICYCLER; ASSEMBLY_SHOVILL; ASSEMBLY_ASSESS; ASSEMBLY_QC } from "$projectDir/modules/assembly"
include { GET_REF_GENOME_BWA_DB; MAPPING; SAM_TO_SORTED_BAM; SNP_CALL; HET_SNP_COUNT; MAPPING_QC } from "$projectDir/modules/mapping"
include { GET_KRAKEN2_DB; TAXONOMY; BRACKEN;TAXONOMY_QC } from "$projectDir/modules/taxonomy"
include { OVERALL_QC } from "$projectDir/modules/overall_qc"
include { GENERATE_SAMPLE_REPORT; GENERATE_OVERALL_REPORT } from "$projectDir/modules/output"

// Add this utility process to ensure 'databases/' exists
process INIT_DB_DIR {
    label 'bash_container'
    label 'farm_local'
    publishDir "${params.db}"
    

    output:
    val "${params.db}", emit: db_dir
    path "do_not_modify", emit: dummy

    
    script:
    """
    mkdir do_not_modify
    """
}

// Main pipeline workflow
workflow {

    main:

    INIT_DB_DIR()
    db_dir_ch = INIT_DB_DIR.out.db_dir
    
    // Get path and prefix of Reference Genome BWA Database, generate from assembly if necessary
    GET_REF_GENOME_BWA_DB(params.ref_genome, db_dir_ch)

    // Get path to Kraken2 Database, download if necessary
    GET_KRAKEN2_DB(params.kraken2_db_remote, db_dir_ch)

    // Bracken too, GET_KRAKEN2_DB(params.kraken2_db_remote, params.db)
 
// Get read pairs into Channel raw_read_pairs_ch
    raw_read_pairs_ch = Channel.fromFilePairs("$params.reads/*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}", checkIfExists: true)

    // Basic input files validation
    // Output into Channel FILE_VALIDATION.out.result
    FILE_VALIDATION(raw_read_pairs_ch)

    // From Channel raw_read_pairs_ch, only output valid reads of samples based on Channel FILE_VALIDATION.out.result
    VALID_READS_ch = FILE_VALIDATION.out.result.join(raw_read_pairs_ch, failOnDuplicate: true)
                        .filter { it[1] == 'PASS' }
                        .map { it[0, 2..-1] }

    // Preprocess valid read pairs
    // Output into Channels PREPROCESS.out.processed_reads & PREPROCESS.out.json
    PREPROCESS(VALID_READS_ch)

    // From Channel PREPROCESS.out.json, provide Read QC status
    // Output into Channels READ_QC.out.bases, READ_QC.out.result, READ_QC.out.report
    READ_QC(PREPROCESS.out.json, params.length_low, params.depth)

    // From Channel PREPROCESS.out.processed_reads, only output reads of samples passed Read QC based on Channel READ_QC.out.result
    READ_QC_PASSED_READS_ch = READ_QC.out.result.join(PREPROCESS.out.processed_reads, failOnDuplicate: true)
                        .filter { it[1] == 'PASS' }
                        .map { it[0, 2..-1] }

    // From Channel READ_QC_PASSED_READS_ch, assemble the preprocess read pairs
    // Output into Channel ASSEMBLY_ch, and hardlink (default) the assemblies to $params.output directory
    switch (params.assembler) {
        case 'shovill':
            ASSEMBLY_ch = ASSEMBLY_SHOVILL(READ_QC_PASSED_READS_ch, params.min_contig_length, params.assembler_thread)
            break

        case 'unicycler':
            ASSEMBLY_ch = ASSEMBLY_UNICYCLER(READ_QC_PASSED_READS_ch, params.min_contig_length, params.assembler_thread)
            break
    }

    // From Channel ASSEMBLY_ch, assess assembly quality
    // Output into Channel ASSEMBLY_ASSESS.out.report
    ASSEMBLY_ASSESS(ASSEMBLY_ch)

    // From Channel ASSEMBLY_ASSESS.out.report and Channel READ_QC.out.bases, provide Assembly QC status
    // Output into Channels ASSEMBLY_QC.out.result & ASSEMBLY_QC.out.report
    ASSEMBLY_QC(
        ASSEMBLY_ASSESS.out.report
        .join(READ_QC.out.bases, failOnDuplicate: true),
        params.contigs,
        params.length_low,
        params.length_high,
        params.depth
    )

    // From Channel READ_QC_PASSED_READS_ch map reads to reference
    // Output into Channel MAPPING.out.sam
    MAPPING(GET_REF_GENOME_BWA_DB.out.path, GET_REF_GENOME_BWA_DB.out.prefix, READ_QC_PASSED_READS_ch)

    // From Channel MAPPING.out.sam, Convert SAM into sorted BAM and calculate reference coverage
    // Output into Channels SAM_TO_SORTED_BAM.out.sorted_bam and SAM_TO_SORTED_BAM.out.ref_coverage
    SAM_TO_SORTED_BAM(MAPPING.out.sam, params.lite)

    // From Channel SAM_TO_SORTED_BAM.out.sorted_bam calculates non-cluster Het-SNP site count
    // Output into Channel HET_SNP_COUNT.out.result
    SNP_CALL(params.ref_genome, SAM_TO_SORTED_BAM.out.sorted_bam, params.lite)
    HET_SNP_COUNT(SNP_CALL.out.vcf)

    // Merge Channels SAM_TO_SORTED_BAM.out.ref_coverage & HET_SNP_COUNT.out.result to provide Mapping QC Status
    // Output into Channels MAPPING_QC.out.result & MAPPING_QC.out.report
    MAPPING_QC(
        SAM_TO_SORTED_BAM.out.ref_coverage
        .join(HET_SNP_COUNT.out.result, failOnDuplicate: true),
        params.ref_coverage,
        params.het_snp_site
    )

    // From Channel READ_QC_PASSED_READS_ch assess Streptococcus agalactiae percentage in reads
    // Output into Channel TAXONOMY.out.report
    TAXONOMY(GET_KRAKEN2_DB.out.path, params.kraken2_memory_mapping, READ_QC_PASSED_READS_ch)
//Bracken too, estimate the abundance of Streptococcus agalactiae
    BRACKEN(
        GET_KRAKEN2_DB.out.path,
        TAXONOMY.out.report,
        params.read_len,
        params.classification_level,
        params.threshold
    )
    // From Channel TAXONOMY.out.report, provide taxonomy QC status
    // Output into Channels TAXONOMY_QC.out.result & TAXONOMY_QC.out.report
    TAXONOMY_QC(BRACKEN.out.bracken_report, params.sagalactiae_percentage, params.top_non_agalactiae_species_percentage)

    
    // Merge Channels FILE_VALIDATION.out.result & READ_QC.out.result & ASSEMBLY_QC.out.result & MAPPING_QC.out.result & TAXONOMY_QC.out.result to provide Overall QC Status
    // Output into Channel OVERALL_QC.out.result & OVERALL_QC.out.report
    OVERALL_QC(
        raw_read_pairs_ch.map{ it[0] }
        .join(FILE_VALIDATION.out.result, failOnDuplicate: true, remainder: true)
        .join(READ_QC.out.result, failOnDuplicate: true, remainder: true)
        .join(ASSEMBLY_QC.out.result, failOnDuplicate: true, remainder: true)
        .join(MAPPING_QC.out.result, failOnDuplicate: true, remainder: true)
        .join(TAXONOMY_QC.out.result, failOnDuplicate: true, remainder: true)
    )

    // From Channel READ_QC_PASSED_READS_ch, only output reads of samples passed overall QC based on Channel OVERALL_QC.out.result
    OVERALL_QC_PASSED_READS_ch = OVERALL_QC.out.result.join(READ_QC_PASSED_READS_ch, failOnDuplicate: true)
                        .filter { it[1] == 'PASS' }
                        .map { it[0, 2..-1] }

    // From Channel ASSEMBLY_ch, only output assemblies of samples passed overall QC based on Channel OVERALL_QC.out.result
    OVERALL_QC_PASSED_ASSEMBLIES_ch = OVERALL_QC.out.result.join(ASSEMBLY_ch, failOnDuplicate: true)
                            .filter { it[1] == 'PASS' }
                            .map { it[0, 2..-1] }

   
    // Generate sample reports by merging outputs from all result-generating modules
    GENERATE_SAMPLE_REPORT(
        raw_read_pairs_ch.map{ it[0] }
        .join(READ_QC.out.report, failOnDuplicate: true, remainder: true)
        .join(ASSEMBLY_QC.out.report, failOnDuplicate: true, remainder: true)
        .join(MAPPING_QC.out.report, failOnDuplicate: true, remainder: true)
        .join(TAXONOMY_QC.out.report, failOnDuplicate: true, remainder: true)
        .join(OVERALL_QC.out.report, failOnDuplicate: true, failOnMismatch: true)
        .map{[ it[0], it[1..-1].minus(null)]}
       
    )
    //GENERATE_OVERALL_REPORT(GENERATE_SAMPLE_REPORT.out.report.collect()) for later combining qc and tyoer reports 

}
