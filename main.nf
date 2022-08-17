#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process split_region {
    input:
        val region
    output:
        tuple val(region), path(bedfile)
    script:
        bedfile = region.replaceAll('\t', '_') + '.atomic.bed'
        """
        echo "${region}" > temp.bed
        bedops --chop 100 temp.bed > "${bedfile}"
        """
}

process mpileup {
    input:
        tuple val(region), path(bedfile)
    output:
        path count_list
    script:
        count_list = region.replaceAll('\t', '_') + '.counts.bed'
        """
        python3 ${projectDir}/count_tags.py "${params.bam_file}" < "${bedfile}" > tags_counts.txt

        echo ${region} > region.txt
        tr '\n' ',' < tags_counts.txt > tags_joined.txt
        paste region.txt tags_joined.txt > ${count_list} 
        """
}

process sort {
    publishDir "output", mode: 'symlink'
    input:
        path(unsorted_collect)
    output:
        path result_bed
    script:
        result_bed = 'counts_by_splits.bed'
        """
        sort-bed ${unsorted_collect} > ${result.bed}
        """
}

workflow {
    input = Channel.fromPath(params.nanosv_regions)
    split_region(input) | mpileup
    merge(mpileup.out.collectFile(name: 'counts_by_splits.tsv', newLine: true))
}