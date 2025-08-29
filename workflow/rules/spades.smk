rule spades_hybrid__run:
    """Run SPAdes over one sample, merging all libraries in the process"""
    input:
        long=get_longreads_from_assembly_id,
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        fasta=SPADES_HYBRID / "{assembly_id}.fa.gz",
        tarball=SPADES_HYBRID / "{assembly_id}.tar.gz",
        concatenated_forwards=temp(SPADES_HYBRID / "{assembly_id}_R1_concat.fastq.gz"),
        concatenated_reverses=temp(SPADES_HYBRID / "{assembly_id}_R2_concat.fastq.gz"),
    log:
        log=SPADES_HYBRID / "{assembly_id}.log",
    container:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["veryhighmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        partition=config["resources"]["partition"]["highlong"],
    params:
        out_dir=lambda w: SPADES_HYBRID / w.assembly_id,
        kmer_size=params["assemble"]["spades"]["kmer_size"],
        additional_options=params["assemble"]["spades"]["additional_options"],
        forwards=aggregate_forwards_for_spades,
        reverses=aggregate_reverses_for_spades,
        assembly_id=lambda w: w.assembly_id,
        memory = int(config["resources"]["mem_per_cpu"]["veryhighmem"] / 1024),
    shell:
        """
        # Concatenate forward reads into a single file
        cat {params.forwards} > {output.concatenated_forwards}
        
        # Concatenate reverse reads into a single file
        cat {params.reverses} > {output.concatenated_reverses}
        
        spades.py \
            --meta \
            -t {threads} \
            --memory {params.memory} \
            --nanopore {input.long} \
            -1 {output.concatenated_forwards} \
            -2 {output.concatenated_reverses} \
            -o {params.out_dir} \
        2> {log} 1>&2
        """
        
rule spades_long__run:
    """Run SPAdes using only the long reads"""
    input:
        long=get_longreads_from_assembly_id,
    output:
        fasta=SPADES_LONG / "{assembly_id}.fa.gz",
        tarball=SPADES_LONG / "{assembly_id}.tar.gz",
        concatenated_forwards=temp(SPADES_LONG / "{assembly_id}_R1_concat.fastq.gz"),
        concatenated_reverses=temp(SPADES_LONG / "{assembly_id}_R2_concat.fastq.gz"),
    log:
        log=SPADES_LONG / "{assembly_id}.log",
    container:
        docker["assemble"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["veryhighmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["longrun"],
        partition=config["resources"]["partition"]["highlong"],
    params:
        out_dir=lambda w: SPADES_LONG / w.assembly_id,
        kmer_size=params["assemble"]["spades"]["kmer_size"],
        additional_options=params["assemble"]["spades"]["additional_options"],
        forwards=aggregate_forwards_for_spades,
        reverses=aggregate_reverses_for_spades,
        assembly_id=lambda w: w.assembly_id,
        memory = int(config["resources"]["mem_per_cpu"]["veryhighmem"] / 1024),
    shell:
        """
        # Concatenate forward reads into a single file
        cat {params.forwards} > {output.concatenated_forwards}
        
        # Concatenate reverse reads into a single file
        cat {params.reverses} > {output.concatenated_reverses}
        
        spades.py \
            --meta \
            -t {threads} \
            --memory {params.memory} \
            --nanopore {input.long} \
            -o {params.out_dir} \
        2> {log} 1>&2
        """        
        
rule spades_hybrid:
    """Rename all assemblies contigs to avoid future collisions"""
    input:
        [SPADES_HYBRID / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],

rule spades_long:
    """Rename all assemblies contigs to avoid future collisions"""
    input:
        [SPADES_LONG / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
