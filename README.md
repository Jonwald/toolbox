# toolbox
Various small scripts for handling seqeunce data files and perfomring common bioinformatics analysis jobs.

to do:
 - add all R scripts for RNAseq / reportable range

## fasta
 
 - fasta_length.py: print name and length of each sequence in a fasta file
 - remove_fasta_wrap.sh: remove line wrapping from a fasta file
 - fake_reads.py: Simulate PE reads in fastq format for a given reference sequence, can control number of fragments, readlength and random error probability
 - sim_reads.py: same as fake reads.py with more options, can control fragment size range, and provide a list of specific SNPs to model

## microbiology
 - get_taxpaths.py: print full lineage of taxon IDs for one or more NCBI taxon IDs
 - names_from_taxids/py: print full lineage of taxon names for one or more NCBI taxon IDs
 - parse_prodigal.py: Parse Prodigal gff output to produce simple TSV output of contig id to gene ID

## misc
 - bwameth_mergeqc.py: merge various qc outputs from bwa-meth pipleine into a single tsv file
 - merge_picard_dup_stats.py: merge qc metrics from picard markdupicates stats files into a single tsv file
 - filter_tsv.py: filter specific rows from a tsv files matching values from a list, useful for for quickly filtering gene expresssion matrices
 - merge_tsv.py: merge specific columns from multiple tsv files on a common key column
 - meth_per_chrom.py: summarise methylation percentage per-cromosome from methyldackel output
 - uniq_counts_codeff.py: summarise the number of unique start-sites for each detected fusion in STAR-fusion output

## RNA-seq
 - rna_saturation.py: work in progress, quickly plot RNA saturation for given gene(s) from a bam file

## TSO500
 - tso_parse.py: parse multiple combinedvariantoutput and metricsoutput files from TSO500 panel, compile into single TSV file for import into LIMs

## VCF
 - cyt_reformat.py: Apply some general germline and somatic filtering variant logic to annotated tumor only vcf files (not to be used as a substitute for matched normals / panel of normals)
 - gc_over_window.py: report average GC% of bam alignments over specific coordinate ranges defined in a bed file
 - reportable_range.py: some functions for calculation reportable range of an assay from BAM alignmnets, (todo: flesh out into re-usable script)
