process your_package {
    container "hw2-liane-m:0.0.2"  // Student: change based on the tag you gave to your image

    input:
    // Student: adjust the input(s) based on the input of your Rapp CLI
    path counts

    output:
    path 'results/'

    // Student: adjust the following block based on how your Rapp CLI is designed
    script:
    """
        corinet network --counts ${counts} --output results/
        corinet heatmap --counts ${counts} --output results/
        corinet plot_network --counts ${counts} --output results/
    """
}
