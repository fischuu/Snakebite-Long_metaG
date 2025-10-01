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
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        partition=config["resources"]["partition"]["small"],
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
        
