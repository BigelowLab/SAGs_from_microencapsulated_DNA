#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.indir ="./input/"
process sayHello {
    input:
    val greeting

    output:
    stdout

    script:
    """
    echo '${greeting} world!'
    """
}

workflow {
    CH_fastq = channel.fromPath("${params.indir}/*.{fastq,fastq.gz,fq,fq.gz}", checkIfExists: true)//channel.of('Bonjour', 'Ciao', 'Hello', 'Hola')
    sayHello(CH_fastq).view()
}

