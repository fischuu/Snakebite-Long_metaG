# YOU CAN ALSO USE STRAINY OUTPUT HERE INSTEAD OF FLYE!!!

rule racon_run:
    """1st Racon polishing round on Flye (or Strainy) assembly"""
    input:
        asm=FLYE_LONG / "{assembly_id}.fa.gz",
        long=get_longreads_from_assembly_id,
    output:
        fa=RACON / "{assembly_id}.racon1.fa.gz",
    log:
        RACON / "{assembly_id}.racon1.log",
    container:
        docker["medaka"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    params:
        tmp=RACON / "{assembly_id}.racon1",
    shell:
        r"""
        mkdir -p {params.tmp}
        minimap2 -t {threads} -x map-ont <(zcat {input.asm}) {input.long} \
          | gzip -1 > {params.tmp}/aln.paf.gz 2> {log}
        racon -t {threads} {input.long} {params.tmp}/aln.paf.gz <(zcat {input.asm}) \
          | gzip -c > {output.fa} 2>> {log}
        """

rule medaka_run:
    """Medaka polishing after Racon"""
    input:
        fa=RACON / "{assembly_id}.racon1.fa.gz",
        long=get_longreads_from_assembly_id,
    output:
        fa=MEDAKA / "{assembly_id}.medaka.fa.gz",
    log:
        MEDAKA / "{assembly_id}.medaka.log",
    container:
        docker["medaka"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    params:
        out=MEDAKA / "{assembly_id}.medaka",
        model=params["assemble"]["medaka"]["model"],  # z.B. r1041_e82_400bps_sup_v4.2.0
    shell:
        r"""
        mkdir -p {params.out}
        # Medaka erwartet BAM/reads+ref; wir lassen es intern mappen:
        medaka_consensus \
          -t {threads} \
          -m {params.model} \
          -i {input.long} \
          -d <(zcat {input.fa}) \
          -o {params.out} 2> {log}
        gzip -c {params.out}/consensus.fasta > {output.fa}
        """
        
        
rule medaka:
    """Collect all Medaka results"""
    input:
        [MEDAKA / f"{assembly_id}.medaka.fa.gz" for assembly_id in ASSEMBLIES],
