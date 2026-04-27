#!/usr/bin/env nextflow

// Student: these should be adjusted based on how your Rapp command-line interface accepts inputs
params.counts = 'data/test/counts.tsv'  // pass your file

include { run_network } from "./modules/coexpression_network.nf"
include { plot_results } from "./modules/coexpression_network.nf"

workflow {
    main:

        ch_counts = channel.fromPath(params.counts, checkIfExists: true)
            .filter { f -> f.name.endsWith('.tsv') || f.name.endsWith('.csv') } // verify that the file is a .tsv or .csv

        run_network(ch_counts)
        plot_results(ch_counts) 

    publish:
        results = run_network.out
}

output {
    results {
        path './results'
        mode 'copy'
    }
}