#!/usr/bin/env nextflow

/*
 *
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */
 

/*
 * Defines the pipeline inputs parameters (giving a default value for each for them) 
 * Each of the following parameters can be specified as command line options
 */
params.query = "$baseDir/data/OmpA.fasta"
params.db = "$baseDir/blast-db/ompa/ompa"
params.out = "results"
params.chunkSize = 100 

workflow {

    /*
     * Create a channel emitting the given query fasta file(s).
     * Split the file into chunks containing as many sequences as defined by the parameter 'chunkSize'.
     * Finally, assign the resulting channel to the variable 'ch_fasta'
     */
    Channel
        .fromPath(params.query)
        .splitFasta(by: params.chunkSize, file:true)
        .set { ch_fasta }

    db_name = file(params.db).name
    db_dir = file(params.db).parent

    /*
     * Execute a BLAST job for each chunk emitted by the 'ch_fasta' channel
     * and emit the resulting BLAST matches.
     */
    blast_data = blast(ch_fasta, db_name, db_dir)
    ch_hits = top_hits(blast_data)

    /*
     * Each time a file emitted by the 'blast' process, an extract job is executed,
     * producing a file containing the matching sequences.
     */
    ch_sequences = extract(ch_hits, db_name, db_dir)

    /*
     * Collect all the sequences files into a single file
     * and print the resulting file contents when complete.
     */
   
    ch_sequences
        .collectFile(name: params.out)
        .view { file -> "matching sequences:\n ${file.text}" }
    
    gpu_test()
    
}


process blast {
    input:
    path ch_fasta
    val db_name_in
    path db_dir_in

    output:
    path 'blast_result'

    publishDir params.out, mode:'copy', overwrite: true

    """
    blastp -db $db_dir_in/$db_name_in -query $ch_fasta -outfmt 6 > blast_result
    """
}

process top_hits {
    input:
    path blast_result

    output:
    path top_hits

    publishDir params.out, mode:'copy', overwrite: true

    """
    cat $blast_result | head -n 10 | cut -f 2 > top_hits
    """

}


process extract {
    input:
    path 'top_hits'
    val db_name_in
    path db_dir_in

    output:
    path 'sequences'

    publishDir params.out, mode:'copy', overwrite: true

    """
    blastdbcmd -db $db_dir_in/$db_name_in -entry_batch top_hits | head -n 10 > sequences
    """
}

process gpu_test {

    label 'gpu_process'

    """
    ls /projects
    nvidia-smi
    """

}


