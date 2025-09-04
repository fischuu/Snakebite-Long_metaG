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
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["quitehighmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        partition=config["resources"]["partition"]["small"],
    params:
        out=lambda w: MEDAKA / w.assembly_id,
        #out=MEDAKA / "{assembly_id}.medaka",
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
