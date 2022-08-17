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

    tag region

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
        sort-bed ${unsorted_collect} > ${result_bed}
        """
}

workflow {
    input = Channel.fromPath(params.nanosv_regions).splitText()
    .map(it -> it.trim()).take(10)
    split_region(input) | mpileup
    sort(mpileup.out.collectFile(name: 'counts_by_splits.tsv', newLine: true))
}