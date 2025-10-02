rule compare_assemblies__quast:
    """Compare all assemblies for one sample using metaQUAST"""
    input:
        flye=FLYE_LONG / "{assembly_id}.fa.gz",
        polypolish=POLYPOLISH / "{assembly_id}.polypolish.fa.gz",
        operams=OPERAMS_HYBRID / "{assembly_id}.fa.gz",
        spades=SPADES_HYBRID / "{assembly_id}.fa.gz",
        pilon = PILON / "{assembly_id}.pilon.fa.gz",
        medaka = MEDAKA / "{assembly_id}.medaka.fa.gz",
        racon = RACON / "{assembly_id}.racon.fa",
    output:
        report=COMPARE / "{assembly_id}/report.html",
    log:
        COMPARE / "{assembly_id}/quast.log",
    container:
        docker["metaquast"]
    threads: esc("cpus", "compare_assemblies__quast")
    resources:
        runtime=esc("runtime", "compare_assemblies__quast"),
        mem_mb=esc("mem_mb", "compare_assemblies__quast"),
        cpus_per_task=esc("cpus", "compare_assemblies__quast"),
        slurm_partition=esc("partition", "compare_assemblies__quast"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'compare_assemblies__quast')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("compare_assemblies__quast"))
    params:
        out_dir=lambda w: COMPARE / w.assembly_id,
        racon_rounds=lambda w: " ".join(
            str(RACON / w.assembly_id / f"racon_{i}.fa.gz")
            for i in range(1, params["assemble"]["racon"]["polishing_rounds"] + 1)
        ),
    shell:
        r"""
        mkdir -p {params.out_dir}
        metaquast.py -t {threads} \
            {input.flye} \
            {input.polypolish} \
            {input.operams} \
            {input.spades} \
            {input.pilon} \
            {input.medaka} \
            {params.racon_rounds} \
            -o {params.out_dir} > {log} 2>&1
        """

rule metaquast:
    """Collect all Strainy results"""
    input:
        [COMPARE / f"{assembly_id}/report.html" for assembly_id in ASSEMBLIES],
