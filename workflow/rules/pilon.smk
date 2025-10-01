# FOR PILON WE NEED TO ALIGN EACH SAMPLE AGAINST THE GENOME WITH BWA AND GET THE BAM
# WE CAN THEN FEED IT IN WITH MULTIPLE --frags OPTIONS
# SEE: https://github.com/broadinstitute/pilon/wiki/Requirements-&-Usage


rule pilon__run:
    """Short-read polishing with Pilon"""
    input:
        asm=MEDAKA / "{assembly_id}.medaka.fa.gz",
        bam=BWA / "{assembly_id}" / "aln.bam",
    output:
        fa=PILON / "{assembly_id}" / "{assembly_id}.pilon.fa.gz",
    log:
        PILON / "{assembly_id}.pilon.log",
    container:
        docker["pilon"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        partition=config["resources"]["partition"]["small"],
    params:
        memory=params["assemble"]["pilon"]["memory"],
        outdir=PILON / "{assembly_id}",
        pre="{assembly_id}",
    shell:
        r"""
        
        java -Xmx{params.memory} -jar pilon.jar \
             --genome {input.asm} \
             --frags short_sorted.bam \
             --out {params.pre} \
             --outdir {params.outdir} \
             --threads {threads} \
             2> {log} 1>&2
        """

rule pilon:
    """Collect all Medaka results"""
    input:
        [PILON / f"{assembly_id}" / f"{assembly_id}.pilon.fa.gz" for assembly_id in ASSEMBLIES],
