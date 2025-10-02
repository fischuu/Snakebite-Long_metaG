rule bwa_index__run:
    """Create the bwa index for aligning the reads"""
    input:
        asm=MEDAKA / "{assembly_id}.medaka.fa.gz",
    output:
        mock=directory(BWA_INDEX / "{assembly_id}"),
    log:
        BWA_INDEX / "{assembly_id}.log",
    container:
        docker["polypolish"]
    threads: esc("cpus", "bwa_index__run")
    resources:
        runtime=esc("runtime", "bwa_index__run"),
        mem_mb=esc("mem_mb", "bwa_index__run"),
        cpus_per_task=esc("cpus", "bwa_index__run"),
        slurm_partition=esc("partition", "bwa_index__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'bwa_index__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("bwa_index__run"))
    params:
        tmp=lambda w: BWA_INDEX / w.assembly_id,
    shell:
        r"""
        
        set -x 
        mkdir -p {params.tmp}
        
        gunzip -c {input.asm} > {params.tmp}/medaka.fa 2> {log}

        bwa index {params.tmp}/medaka.fa 2>> {log} 1>&2
        
        echo "Exit code bwa: $?" 2>> {log} 1>&2
        """
        
rule bwa_index__build:
    """Index all medaka assemblies"""
    input:
        [BWA_INDEX / f"{assembly_id}" for assembly_id in ASSEMBLIES],

rule bwa_pe__run:
    """Align paired short reads to medaka assembly"""
    input:
        asm=MEDAKA / "{assembly_id}.medaka.fa.gz",
        forward_=FASTP / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=FASTP / "{sample_id}.{library_id}_2.fq.gz",
        mock=BWA_INDEX / "{assembly_id}",
    output:
        bam=BWA_PAIRED / "{assembly_id}" / "{sample_id}.{library_id}.bam",
        unsorted=temp(BWA_PAIRED / "{assembly_id}" / "{sample_id}.{library_id}_unsorted.bam"),
    log:
        BWA_PAIRED / "{assembly_id}" / "{sample_id}.{library_id}.bwa.log",
    container:
        docker["polypolish"]
    threads: esc("cpus", "bwa_pe__run")
    resources:
        runtime=esc("runtime", "bwa_pe__run"),
        mem_mb=esc("mem_mb", "bwa_pe__run"),
        cpus_per_task=esc("cpus", "bwa_pe__run"),
        slurm_partition=esc("partition", "bwa_pe__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'bwa_pe__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("bwa_pe__run"))    
    params:
        index=lambda w: BWA_INDEX / w.assembly_id,
        tmp=lambda w: BWA_PAIRED / w.assembly_id,
    shell:
        r"""
        mkdir -p {params.tmp}
        
        bwa mem -t {threads} {params.index}/medaka.fa {input.forward_} {input.reverse_} > {output.unsorted} 2>> {log}
        
        samtools sort {output.unsorted} -o {output.bam} 2>> {log}

        samtools index {output.bam} 2>> {log}
        
        """

rule bwa_pe:
    """Map all samples to all the assemblies that they belong to"""
    input:
        [
            BWA_PAIRED / f"{assembly_id}" / f"{sample_id}.{library_id}.bam"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],


rule bwa_se__run:
    """Align paired short reads one by one to medaka assembly"""
    input:
        asm=MEDAKA / "{assembly_id}.medaka.fa.gz",
        forward_=FASTP / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=FASTP / "{sample_id}.{library_id}_2.fq.gz",
        mock=BWA_INDEX / "{assembly_id}",
    output:
        sam_1=BWA_SINGLE / "{assembly_id}" / "{sample_id}.{library_id}_1.sam",
        sam_2=BWA_SINGLE / "{assembly_id}" / "{sample_id}.{library_id}_2.sam",
    log:
        BWA_SINGLE / "{assembly_id}" / "{sample_id}.{library_id}.bwa.log",
    container:
        docker["polypolish"]
    threads: esc("cpus", "bwa_se__run")
    resources:
        runtime=esc("runtime", "bwa_se__run"),
        mem_mb=esc("mem_mb", "bwa_se__run"),
        cpus_per_task=esc("cpus", "bwa_se__run"),
        slurm_partition=esc("partition", "bwa_se__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'bwa_se__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("bwa_se__run"))
    params:
        index=lambda w: BWA_INDEX / w.assembly_id,
        tmp=lambda w: BWA_SINGLE / w.assembly_id,
    shell:
        r"""
        mkdir -p {params.tmp}
        
        bwa mem -t {threads} -a {params.index}/medaka.fa {input.forward_} > {output.sam_1} 2>> {log}
        bwa mem -t {threads} -a {params.index}/medaka.fa {input.reverse_} > {output.sam_2} 2>> {log}
        """

rule bwa_se:
    """Map all samples one by one to all the assemblies that they belong to"""
    input:
        [
            BWA_SINGLE / f"{assembly_id}" / f"{sample_id}.{library_id}_1.sam"
            for assembly_id, sample_id, library_id in ASSEMBLY_SAMPLE_LIBRARY
        ],
