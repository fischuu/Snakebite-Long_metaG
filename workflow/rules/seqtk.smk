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
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["quitehighmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        partition=config["resources"]["partition"]["small"],
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
