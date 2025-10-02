rule seqtk__run:
    """Filter the flye output for length"""
    input:
        long=FLYE_LONG / "{assembly_id}.fa.gz",
    output:
        fasta=SEQTK / "{assembly_id}.length_filtered.fa.gz",
    log:
        log=SEQTK / "{assembly_id}.log",
    container:
        docker["seqtk"]
    threads: esc("cpus", "seqtk__run")
    resources:
        runtime=esc("runtime", "seqtk__run"),
        mem_mb=esc("mem_mb", "seqtk__run"),
        cpus_per_task=esc("cpus", "seqtk__run"),
        slurm_partition=esc("partition", "seqtk__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'seqtk__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("seqtk__run"))
    params:
        out_dir=lambda w: FLYE_LONG / w.assembly_id,
        additional_options=params["assemble"]["seqtk"]["kept_length"],
    shell:
        """
        seqtk seq -L {params.additional_options} -l 80 {input} > {params.out_dir}/filtered.fasta 2> {log}

        # Compress output contigs
        gzip -c {params.out_dir}/filtered.fasta > {output}
        """


rule seqtk:
    """Remove too short sequence from the flye assembly"""
    input:
        [SEQTK / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
