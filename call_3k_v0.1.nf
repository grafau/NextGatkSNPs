#!/usr/bin/env nextflow

// ######################################
// #    Rafal Gutaker, NYU, 2018        #
// #    Pipeline that performs:         #
// #    - splitting gvcf per chromosome #
// #    - Genotypes chromosomes		#
// #                                    #
// ######################################


params.help=""
params.exe="local"

if (params.help) {
        log.info " "
        log.info "This is Nextflow pipeline for genotyping SNPS in rice 3k dataset. Most parameters are stored in file call_3k.conf"
        log.info " "
        log.info "Usage: nextflow run call_3k.nf -c call_3k.conf --ref --list --(other options)"
        log.info "Options:"
        log.info "--help\t[BOOLEAN]\tShow this help message"
        log.info "--ref\t[STRING]\tPath to the indexed referece fasta file [OBLIGATORY]"
        log.info "--list\t[STRING]\tPath to the file with run IDs and fastq files to be processed [OBLIGATORY]"
        log.info "--exe\t[STRING]\tExecutor mode, -local- or -slurm- [DEFUALT: local]"
        log.info " "
        exit 1
}

// Initalize Input
DAT = file("${params.list}")
REF = file("${params.ref}")

// Create Input Channel
SampleData = Channel.fromPath("${DAT}").splitCsv(header: ['CHR','BED','FILE'], skip: 0, by:1, sep:",")


// ########### Split into chromosomes  #############
process split {

	errorStrategy 'retry'
	maxRetries 2
	maxErrors 13	

        executor = "${params.exe}"
	if ("${params.exe}" == "slurm")
	{
        	clusterOptions = "--qos=cpuplus --cpus-per-task=${params.sv_cpu} --time=${params.sv_rt} --mem=${params.sv_vmem}"
	}

        input:
         val(GVCF) from SampleData

        output:
         set val(CHR), file({ "${CHUNK}" }) into splat

        script:

	 CHR = "${GVCF.CHR}"
	 BED = file("${GVCF.BED}")
         CHUNK = "${GVCF.FILE}_${GVCF.CHR}.g.vcf"

         IN = file("${GVCF.FILE}")

         """
         ${params.sv_bin} ${params.sv_param} -V ${IN} -R ${REF} -L ${BED} -o ${CHUNK}
	 """
}


// ####### MERGE CHROMOSOME CHANNELS #########


splat
	.map {chr, file -> tuple(chr, file) }
	.groupTuple()
	.set { merg_chrom }



// ########### Genotype GVCF ###############
process genotype {

	errorStrategy 'retry'  
        maxRetries 2
        maxErrors 2

        executor = "${params.exe}"
        if ("${params.exe}" == "slurm")
        {
        	clusterOptions = "--qos=cpuplus --cpus-per-task=${params.gv_cpu} --time=${params.gv_rt} --mem=${params.gv_vmem}"
	}

        input:
         set val(CHR), file(CHUNK) from merg_chrom

        output:
         file({ "${CALLD}" }) into called

        script:

         ALL_IN = CHUNK.collect { "-V $it" }.join(' ')
         CALLD = "${CHR}.vcf"

         """
         ${params.gv_bin} ${params.gv_param} ${ALL_IN} -R ${REF} -o ${CALLD}
         """
}


// Return filtered gVCFs
called.subscribe { println "[OUT] $it" }



