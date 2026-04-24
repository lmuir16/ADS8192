#!/usr/bin/env nextflow

// Pasted from:
// https://github.com/nextflow-io/training/blob/master/nextflow-run/1-hello.nf


// Student: these should be adjusted based on how your Rapp command-line interface accepts inputs
params.counts = 'counts.txt'  // pass your file with: nextflow run main.nf --counts ../r-package/tests/testdata/counts.tsv

include { your_package } from "./modules/your_package.nf"

workflow {
    main:
        // Student: you may need to adjust the following line based on the params (input argument) you have above
        ch_counts = channel.fromPath(params.counts, checkIfExists: true)
        your_package(ch_counts)

        // Student: if "your_package" needs more than one inputs, e.g., counts + metadata instead of just counts, you would do
        // ch_counts  = channel.fromPath(params.counts, checkIfExists: true)
        // ch_meta    = channel.fromPath(params.metadata, checkIfExists: true)
        // your_package(ch_counts, ch_meta)


    publish:
        results = your_package.out
}

output {
    results {
        path './results'
        mode 'copy'
    }
}