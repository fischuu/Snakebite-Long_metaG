rule compare_assemblies__quast:
    """Compare all assemblies for one sample using metaQUAST"""
    input:
        flye=FLYE_LONG / "{assembly_id}.fa.gz",
        strainy=STRAINY / "{assembly_id}.fa.gz",
        polypolish=POLYPOLISH / "{assembly_id}.polypolish.fa.gz",
        operams=OPERAMS_HYBRID / "{assembly_id}.fa.gz",
        spades=SPADES_HYBRID / "{assembly_id}.fa.gz",
    output:
        report=COMPARE / "{assembly_id}/report.html",
    log:
        COMPARE / "{assembly_id}/quast.log",
    container:
        docker["metaquast"]
    threads: 8
    params:
        out_dir=lambda w: COMPARE / w.assembly_id,
    shell:
        r"""
        mkdir -p {params.out_dir}
        metaquast.py -t {threads} \
            <(zcat {input.flye}) \
            <(zcat {input.strainy}) \
            <(zcat {input.polypolish}) \
            <(zcat {input.operams}) \
            <(zcat {input.spades}) \
            -o {params.out_dir} > {log} 2>&1

        # metaQUAST schreibt report.html in out_dir
        cp {params.out_dir}/report.html {output.report}
        """

rule metaquast:
    """Collect all Strainy results"""
    input:
        [COMPARE / f"{assembly_id}/report.html" for assembly_id in ASSEMBLIES],
