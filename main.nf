#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

process split_region {
    tag "${region}"

    input:
        val region
    output:
        tuple val(region), path(bedfile)
    script:
        bedfile = region + '.atomic.bed'
        region_bed = region.replaceAll('_', '\t')
        """
        echo "${region_bed}" > temp.bed
        bedops --chop "${params.slice_size}" temp.bed | bedops -n 1 - "${params.blacklist}" > "${bedfile}"
        """
}

process mpileup {
    tag "${region}"

    input:
        tuple val(region), path(bedfile)
    output:
        path count_list
    script:
        count_list = region + '.counts.bed'
        region_bed = region.replaceAll('_', '\t')
        """
        python3 ${projectDir}/count_tags.py "${params.bam_file}" < "${bedfile}" > tags_counts.txt

        paste "${bedfile}" tags_counts.txt > ${count_list} 
        """
}

process sort {
    publishDir params.outdir, mode: 'symlink'
    input:
        path(unsorted_collect)
    output:
        path result_bed
    script:
        result_bed = 'counts_by_splits.bed'
        """
        sort -k 1,1 -k2,2n ${unsorted_collect} > ${result_bed}
        """
}

workflow {
    input = Channel.fromPath(params.nanosv_regions).splitText()
    .map(it -> it.trim().replaceAll('\t', '_'))
    regions = split_region(input)
    num_regions = regions.count()
    mp = regions.take(num_regions - 1) | mpileup
    sort(mp.collectFile(name: 'counts_by_splits.tsv'))
}