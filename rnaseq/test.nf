params.reads

process FASTQC {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',

    input:
    set val(name), file(reads) from trimmed_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

process STAR {  
    tag "$name"
    publishDir "${params.outdir}/star", mode: 'copy',

    input:
    set val(name), file(reads) from trimmed_reads_star

    output:
    file "*.{bam,bai}" into star_results

    script:
    """
    STAR --runThreadN ${task.cpus} --genomeDir ${params.star_index} --readFilesIn $reads --outFileNamePrefix ${name}_ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes All --outSAMattrRGline ID:${name} SM:${name} PL:ILLUMINA --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --sjdbScore 1 --readFilesCommand zcat
    """
}

proces kallisto {
    tag "$name"
    publishDir "${params.outdir}/kallisto", mode: 'copy',

    input:
    set val(name), file(reads) from trimmed_reads_kallisto

    output:
    file "*.{abundance.tsv,abundance.h5}" into kallisto_results

    script:
    """
    kallisto quant -i ${params.kallisto_index} -o ${name}_ -t ${task.cpus} $reads
    """
}

process RSeQC {
    tag "$name"
    publishDir "${params.outdir}/rseqc", mode: 'copy',

    input:
    set val(name), file(reads) from trimmed_reads_rseqc

    output:
    file "*.{pdf,txt}" into rseqc_results

    script:
    """
    bam_stat.py -i ${name}_Aligned.sortedByCoord.out.bam
    geneBody_coverage.py -i ${name}_Aligned.sortedByCoord.out.bam -r ${params.rseqc_gtf} -o ${name}_geneBodyCoverage
    inner_distance.py -i ${name}_Aligned.sortedByCoord.out.bam -r ${params.rseqc_gtf} -o ${name}_innerDistance
    junction_annotation.py -i ${name}_Aligned.sortedByCoord.out.bam -r ${params.rseqc_gtf} -o ${name}_junctionAnnotation
    junction_saturation.py -i ${name}_Aligned.sortedByCoord.out.bam -r ${params.rseqc_gtf} -o ${name}_junctionSaturation
    read_distribution.py -i ${name}_Aligned.sortedByCoord.out.bam -r ${params.rseqc_gtf} -o ${name}_readDistribution
    read_duplication.py -i ${name}_Aligned.sortedByCoord.out.bam -o ${name}_readDuplication
    read_GC.py -i ${name}_Aligned.sortedByCoord.out.bam -o ${name}_readGC
    read_NVC.py -i ${name}_Aligned.sortedByCoord.out.bam -o ${name}_readNVC
    read_quality.py -i ${name}_Aligned.sortedByCoord.out.bam -o ${name}_readQuality
    """

}