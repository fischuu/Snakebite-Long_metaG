# helper to collect BAMs for a given assembly
def bams_for_assembly(wildcards):
    seen = set()
    uniq_bams = []
    for a, s, l in ASSEMBLY_SAMPLE_LIBRARY:
        if a == wildcards.assembly_id:
            bam = str(BWA_PAIRED / a / f"{s}.{l}.bam")
            if bam not in seen:
                uniq_bams.append(bam)
                seen.add(bam)
    return uniq_bams



rule pilon_run:
    """Short-read polishing with Pilon"""
    input:
        asm = MEDAKA / "{assembly_id}.medaka.fa.gz",
        bams = bams_for_assembly,
    output:
        fa = PILON / "{assembly_id}" / "{assembly_id}.pilon.fa.gz",
    log:
        PILON / "{assembly_id}" / "{assembly_id}.pilon.log",
    container:
        docker["pilon"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["highmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        partition=config["resources"]["partition"]["small"],
    params:
        memory = params["assemble"]["pilon"]["memory"],
        outdir = lambda w: PILON / w.assembly_id,
        pre = lambda w: w.assembly_id,
    shell:"""
        # build --frags args from the input BAMs
        frags=""
        for bam in {input.bams}; do
            frags="$frags --frags $bam"
        done
      
        gunzip -c {input.asm} > {params.outdir}/asm.fa
  
        # Export Java memory options
        export _JAVA_OPTIONS="-Xmx{params.memory} -Xms4g"
        export JAVA_TOOL_OPTIONS="-Xmx{params.memory} -Xms4g"
      
        pilon --genome {params.outdir}/asm.fa \
              $frags \
              --output {params.pre} \
              --outdir {params.outdir} \
              --threads {threads} \
              2> {log}
        
        gzip -c {params.outdir}/{params.pre}.fasta > {output.fa}
        rm {params.outdir}/{params.pre}.fasta
        """


rule pilon:
    """Collect all Pilon-polished assemblies"""
    input:
        [PILON / f"{assembly_id}" / f"{assembly_id}.pilon.fa.gz" for assembly_id in ASSEMBLIES],
