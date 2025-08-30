# YOU CAN ALSO USE STRAINY OUTPUT HERE INSTEAD OF FLYE!!!

rule racon_run:
    """1st Racon polishing round on Flye (or Strainy) assembly"""
    input:
        asm=FLYE_LONG / "{assembly_id}.fa.gz",
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
    shell:
        r"""
        mkdir -p {params.tmp}
        minimap2 -t {threads} -x map-ont <(zcat {input.asm}) {input.long} \
           > {params.tmp}/aln.paf 2> {log}
        racon -t {threads} {input.long} {params.tmp}/aln.paf {input.asm} \
           > {output.fa} 2>> {log}
        """


rule racon:
    """Collect all Racon results"""
    input:
        expand(RACON / "{assembly_id}.racon.fa" , assembly_id=ASSEMBLIES),
        
