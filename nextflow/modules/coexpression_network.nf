// Run gene network analysis (corinet processes) on counts data
process run_network {
    container "net-lmuir:0.0.1"

    input:
    path counts

    output:
    path 'results/'

    script:
    """
        corinet network --counts ${counts} --output results/ --n-top 50
    """
}

// Plot results (heatmap and network plot)
process plot_results {
    container "net-lmuir:0.0.1"

    input:
    path counts

    output:
    path 'results/'

    script:
    """
        corinet heatmap --counts ${counts} --output results/ --n-top 50
        corinet plot-network --counts ${counts} --output results/ --n-top 50 --n-plot 20
    """
}
