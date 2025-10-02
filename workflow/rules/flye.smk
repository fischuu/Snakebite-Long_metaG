rule flye__run:
    """Run Flye over one sample using long reads only"""
    input:
        long=get_longreads_from_assembly_id,
    output:
        fasta=FLYE_LONG / "{assembly_id}.fa.gz",
        gfa=FLYE_LONG / "{assembly_id}" / "assembly_graph.gfa",
    log:
        log=FLYE_LONG / "{assembly_id}.log",
    container:
        docker["flye"]
    threads: esc("cpus", "flye__run")
    resources:
        runtime=esc("runtime", "flye__run"),
        mem_mb=esc("mem_mb", "flye__run"),
        cpus_per_task=esc("cpus", "flye__run"),
        slurm_partition=esc("partition", "flye__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'flye__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("flye__run"))
    params:
        out_dir=lambda w: FLYE_LONG / w.assembly_id,
        iterations=params["assemble"]["flye"]["polishing_iterations"],
        additional_options=params["assemble"]["flye"]["additional_options"],
    shell:
        """
        flye --nano-raw {input.long} \
             --threads {threads} \
             --out-dir {params.out_dir} \
             --meta \
             --iterations {params.iterations} \
             {params.additional_options} \
              2> {log} 1>&2

        # Compress output contigs
        gzip -c {params.out_dir}/assembly.fasta > {output.fasta}
        """


rule flye:
    """Collect long-read Flye assemblies"""
    input:
        [FLYE_LONG / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
