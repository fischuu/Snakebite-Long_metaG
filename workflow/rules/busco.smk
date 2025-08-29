rule busco__compare:
    input:
        fa=COMPARE / "{assembly_id}/report.html",   # hängt an quast
    output:
        summary=COMPARE / "{assembly_id}/busco/short_summary.txt",
    log:
        COMPARE / "{assembly_id}/busco.log",
    container:
        docker["busco"]
    threads: 8
    params:
        out_dir=lambda w: COMPARE / f"{w.assembly_id}/busco",
        lineage=params["qc"]["busco"]["lineage"],   # z. B. bacteria_odb10
    shell:
        r"""
        busco \
            -i {COMPARE}/{wildcards.assembly_id}/../*.fa.gz \
            -l {params.lineage} \
            -o {params.out_dir} \
            -m genome \
            -c {threads} > {log} 2>&1
        """
        
rule busco:
    """Collect all Busco results"""
    input:
        [COMPARE / f"{assembly_id}/busco/short_summary.txt" for assembly_id in ASSEMBLIES],

