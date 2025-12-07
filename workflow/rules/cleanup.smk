rule long_read_assemblies_run:
    """
    Link polished assemblies into final_assemblies/long_reads/
    """
    input:
        fa=POLYPOLISH / "{assembly_id}.polypolish.fa.gz"
    output:
        fa=FINAL_LONG / "{assembly_id}.polypolish.fa.gz"
    shell:
        """
        ln -sf $(realpath {input.fa}) {output.fa}
        """


rule long_read_assemblies:
    """Collect all polished assemblies"""
    input:
        [FINAL_LONG / f"{assembly_id}.polypolish.fa.gz" for assembly_id in ASSEMBLIES],
