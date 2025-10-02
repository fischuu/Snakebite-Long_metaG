rule medaka_run:
    """Medaka polishing after Racon"""
    input:
        fa=RACON / "{assembly_id}.racon.fa",
        long=get_longreads_from_assembly_id,
    output:
        fa=MEDAKA / "{assembly_id}.medaka.fa.gz",
    log:
        MEDAKA / "{assembly_id}.medaka.log",
    container:
        docker["medaka"]
    threads: esc("cpus", "medaka_run")
    resources:
        runtime=esc("runtime", "medaka_run"),
        mem_mb=esc("mem_mb", "medaka_run"),
        cpus_per_task=esc("cpus", "medaka_run"),
        slurm_partition=esc("partition", "medaka_run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'medaka_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("medaka_run"))
    params:
        out=lambda w: MEDAKA / w.assembly_id,
        model=params["assemble"]["medaka"]["model"],  # z.B. r1041_e82_400bps_sup_v4.2.0
    shell:
        r"""
        mkdir -p {params.out}

        medaka_consensus \
          -t {threads} \
          -m {params.model} \
          -i {input.long} \
          -d {input.fa} \
          -o {params.out} 2> {log}
        
        gzip -c {params.out}/consensus.fasta > {output.fa}
        """

rule medaka:
    """Collect all Medaka results"""
    input:
        [MEDAKA / f"{assembly_id}.medaka.fa.gz" for assembly_id in ASSEMBLIES],
