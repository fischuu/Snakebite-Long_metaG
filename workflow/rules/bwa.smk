rule bwa__run:
    """Align short reads to medaka assembly"""
    input:
        asm=MEDAKA / "{assembly_id}.medaka.fa.gz",
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        bam=BWA / "{assembly_id}" / "aln.bam",
    log:
        BWA / "{assembly_id}.bwa.log",
    container:
        docker["polypolish"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        partition=config["resources"]["partition"]["small"],
    params:
        tmp=lambda w: BWA / w.assembly_id,
        forwards=aggregate_forwards_for_bwa,
        reverses=aggregate_reverses_for_bwa,
    shell:
        r"""
        mkdir -p {params.tmp}
        
        gunzip -c {input.asm} > {params.tmp}/medaka.fa

        bwa index {params.tmp}/medaka.fa 2> {log}

        bwa mem -t {threads} {params.tmp}/medaka.fa {input.forwards} {input.reverses} > {params.tmp}/aln.bam 2>> {log}
        
        samtools sort {params.tmp}/aln.bam -o {params.tmp}/aln_sorted.bam 2>> {log}

        samtools index {params.tmp}/aln_sorted.bam 2>> {log}
        
        rm {params.tmp}/aln.bam
        rm {params.tmp}/medaka.fa
        """
