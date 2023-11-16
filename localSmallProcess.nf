params.inputs="$baseDir/fasta_files/*.fasta"

process smallProcess {
    
    executor = "local"
    memory = '500M'    

    input:
      path 'inputFile'

    output:
      stdout emit: headerInfo

    """
    head -n 2 inputFile
    """

}

workflow {
    Channel.fromPath(params.inputs, checkIfExists: true) \
      | smallProcess 
    smallProcess.out.headerInfo.view()
}
