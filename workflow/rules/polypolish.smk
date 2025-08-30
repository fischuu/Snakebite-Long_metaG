rule polypolish__run:
    """Short-read polishing with Polypolish"""
    input:
        asm=MEDAKA / "{assembly_id}.medaka.fa.gz",
        r1=get_forwards_from_assembly_id,
        r2=get_reverses_from_assembly_id,
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
        mkdir -p {params.tmp}
        # Stringentes Mapping für Polypolish:
        bwa index <(zcat {input.asm}) 2> {log} || true
        bwa mem -t {threads} <(zcat {input.asm}) {input.r1} {input.r2} \
          | samtools sort -@ {threads} -o {params.tmp}/aln.bam
        samtools index {params.tmp}/aln.bam
        polypolish {params.tmp}/aln.bam <(zcat {input.asm}) \
          | gzip -c > {output.fa} 2>> {log}
        """

rule polypolish:
    """Collect all Medaka results"""
    input:
        [POLYPOLISH / f"{assembly_id}.polypolish.fa.gz" for assembly_id in ASSEMBLIES],
