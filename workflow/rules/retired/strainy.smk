rule strainy__run:
    """Run Strainy on Flye assembly to resolve strain variation"""
    input:
        assembly=FLYE_LONG / "{assembly_id}.fa.gz",
        long=get_longreads_from_assembly_id,
        gfa=FLYE_LONG / "{assembly_id}" / "assembly_graph.gfa",
    output:
        fasta=STRAINY / "{assembly_id}.fa.gz",
    log:
        log=STRAINY / "{assembly_id}.log",
    container:
        docker["flye"]
    threads: esc("cpus", "strainy__run")
    resources:
        runtime=esc("runtime", "strainy__run"),
        mem_mb=esc("mem_mb", "strainy__run"),
        cpus_per_task=esc("cpus", "strainy__run"),
        slurm_partition=esc("partition", "strainy__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'strainy__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("strainy__run"))
    params:
        out_dir=lambda w: STRAINY / w.assembly_id,
        additional_options=params["assemble"]["strainy"]["additional_options"],
    shell:
        """
        mkdir -p {params.out_dir}

        # Run Strainy
        strainy \
            --gfa_ref {input.gfa} \
            --fastq {input.long} \
            --threads {threads} \
            --output {params.out_dir} \
            -m nano \
            {params.additional_options} \
        2> {log} 1>&2

        # Compress main output assembly
        gzip -c {params.out_dir}/strain_contigs.gfa > {output.fasta}
        """
        

rule strainy:
    """Collect all Strainy results"""
    input:
        [STRAINY / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
