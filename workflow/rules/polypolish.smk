# RUN POLYPOLISH THAT WAY
# FOR POLYPOLISH, WE NEED TO ALIGN ALL FASTQ SEPARATED WITH THE -a OPTION!!!! SO; ALIGN TO ALL POSITIONS!!!!

# see: https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish

bwa index draft.fasta
bwa mem -t 16 -a draft.fasta reads_a_1.fastq.gz > alignments_a_1.sam
bwa mem -t 16 -a draft.fasta reads_a_2.fastq.gz > alignments_a_2.sam
bwa mem -t 16 -a draft.fasta reads_b_1.fastq.gz > alignments_b_1.sam
bwa mem -t 16 -a draft.fasta reads_b_2.fastq.gz > alignments_b_2.sam
polypolish filter --in1 alignments_a_1.sam --in2 alignments_a_2.sam --out1 filtered_a_1.sam --out2 filtered_a_2.sam
polypolish filter --in1 alignments_b_1.sam --in2 alignments_b_2.sam --out1 filtered_b_1.sam --out2 filtered_b_2.sam
polypolish polish draft.fasta filtered_*.sam > polished.fasta


rule polypolish__run:
    """Short-read polishing with Polypolish"""
    input:
        asm=MEDAKA / "{assembly_id}.medaka.fa.gz",
        fa=BWA / "{assembly_id}" / "aln.bam",
    output:
        fa=POLYPOLISH / "{assembly_id}.polypolish.fa.gz",
    log:
        POLYPOLISH / "{assembly_id}.polypolish.log",
    container:
        docker["polypolish"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        partition=config["resources"]["partition"]["small"],
    params:
        tmp=POLYPOLISH / "{assembly_id}.pp",
    shell:
        r"""
        polypolish {input.fa} <(zcat {input.asm}) | gzip -c > {output.fa} 2>> {log}
        """

rule polypolish:
    """Collect all Medaka results"""
    input:
        [POLYPOLISH / f"{assembly_id}.polypolish.fa.gz" for assembly_id in ASSEMBLIES],
