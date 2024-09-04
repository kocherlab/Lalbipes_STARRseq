configfile: 'LALB_v3.1.yml'

workdir: config['workdir']

rule all:
	input:
		expand(os.path.join(config['paths']['rnaseq_aligned_bam_dir'], "{sample}.Aligned.out.bam"), sample=config['samples']),
		expand(os.path.join(config['paths']['rnaseq_sorted_bam_dir'], "{sample}.sortedByCoord.out.bam"), sample=config['samples']),
		os.path.join(config['paths']['rnaseq_bam_dir'], f"{config['species']}_{config['assembly_version']}.bam"),
		os.path.join(config['paths']['assembly_dir'], f"{config['species']}_{config['assembly_version']}.fa.masked"),
		os.path.join(config['paths']['annotations_dir'], 'BRAKER3', 'braker.gff3'),
		os.path.join(config['paths']['downloads_dir'], 'compleasm', f"{config['busco_database']}.done"),
		os.path.join(config['paths']['assembly_dir'], 'compleasm', 'summary.txt'),
		os.path.join(config['paths']['annotations_dir'], 'BUSCO', 'summary.txt')

include: "snakemake_files/rna_seq_hisat2"
include: "snakemake_files/rna_seq_sort"
include: "snakemake_files/rna_seq_merge"
include: "snakemake_files/repeatmasking"
include: "snakemake_files/annotate_braker3"
include: "snakemake_files/download_compleasm"
include: "snakemake_files/assembly_compleasm"
include: "snakemake_files/annotations_busco"
