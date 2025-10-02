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
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        partition=config["resources"]["partition"]["small"],
    params:
        out_dir=lambda w: COMPARE / w.assembly_id,
        racon_rounds=lambda w: ",".join(
            str(RACON / w.assembly_id / f"racon_{i}.fa.gz")
            for i in range(1, params["assemble"]["racon"]["polishing_rounds"] + 1)
        ),
    shell:
        r"""
        mkdir -p {params.out_dir}
        metaquast.py -t {threads} \
            -r {input.flye} \
            -r {input.polypolish} \
            -r {input.operams} \
            -r {input.spades} \
            -r {input.pilon} \
            -r {input.medaka} \
            -r {params.racon_rounds} \
            -o {params.out_dir} > {log} 2>&1

        # metaQUAST schreibt report.html in out_dir
        cp {params.out_dir}/report.html {output.report}
        """

rule metaquast:
    """Collect all Strainy results"""
    input:
        [COMPARE / f"{assembly_id}/report.html" for assembly_id in ASSEMBLIES],
