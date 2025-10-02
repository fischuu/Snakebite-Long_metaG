rule racon_run:
    """1st Racon polishing round on Flye (or Strainy) assembly"""
    input:
        asm=SEQTK / "{assembly_id}.length_filtered.fa.gz",
        long=get_longreads_from_assembly_id,
    output:
        fa=RACON / "{assembly_id}.racon.fa",
    log:
        RACON / "{assembly_id}.racon.log",
    container:
        docker["medaka"]
    threads: esc("cpus", "racon_run")
    resources:
        runtime=esc("runtime", "racon_run"),
        mem_mb=esc("mem_mb", "racon_run"),
        cpus_per_task=esc("cpus", "racon_run"),
        slurm_partition=esc("partition", "racon_run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'racon_run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("racon_run"))
    params:
        tmp=lambda w: RACON / w.assembly_id,
        rounds=params["assemble"]["racon"]["polishing_rounds"],
    shell:
        r"""
        mkdir -p {params.tmp}
        input_fa={input.asm}
        
        # Loop through the steps
        for i in $(seq 1 {params.rounds}); do
            minimap2 -t {threads} -x map-ont <(zcat $input_fa) {input.long} > {params.tmp}/aln_$i.paf 2>> {log}
            racon -t {threads} {input.long} {params.tmp}/aln_$i.paf $input_fa > {params.tmp}/racon_$i.fa 2>> {log}
            bgzip {params.tmp}/racon_$i.fa
            input_fa={params.tmp}/racon_$i.fa.gz
        done

        cp {params.tmp}/racon_{params.rounds}.fa.gz {output.fa}
        """


rule racon:
    """Collect all Racon results"""
    input:
        expand(RACON / "{assembly_id}.racon.fa" , assembly_id=ASSEMBLIES),
        
