rule operams__run:
    """Run Opera-MS hybrid assembly (long + short reads)"""
    input:
        long=get_longreads_from_assembly_id,
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        fasta=OPERAMS_HYBRID / "{assembly_id}.fa.gz",
        tarball=OPERAMS_HYBRID / "{assembly_id}.tar.gz",
        concatenated_forwards=temp(OPERAMS_HYBRID / "{assembly_id}_R1_concat.fastq"),
        concatenated_reverses=temp(OPERAMS_HYBRID / "{assembly_id}_R2_concat.fastq"),
        uncompressed_long=temp(OPERAMS_HYBRID / "{assembly_id}_longreads.fastq"),
    log:
        log=OPERAMS_HYBRID / "{assembly_id}.log",
    container:
        docker["operams"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    resources:
        cpu_per_task=config["resources"]["cpu_per_task"]["multi_thread"],
        mem_per_cpu=config["resources"]["mem_per_cpu"]["quitehighmem"] // config["resources"]["cpu_per_task"]["multi_thread"],
        time=config["resources"]["time"]["verylongrun"],
        partition=config["resources"]["partition"]["longrun"],
    params:
        out_dir=lambda w: OPERAMS_HYBRID / w.assembly_id,
        additional_options=params["assemble"]["operams"]["additional_options"],
        forwards=aggregate_forwards_for_spades,
        reverses=aggregate_reverses_for_spades,
        cpus=config["resources"]["cpu_per_task"]["multi_thread"],
    shell:
        """
        # Concatenate and uncompress forward reads
        for f in {params.forwards}; do
            if [[ $f == *.gz ]]; then
                gunzip -c $f
            else
                cat $f
            fi
        done > {output.concatenated_forwards}

        # Concatenate and uncompress reverse reads
        for f in {params.reverses}; do
            if [[ $f == *.gz ]]; then
                gunzip -c $f
            else
                cat $f
            fi
        done > {output.concatenated_reverses}

        # Uncompress long reads
        if [[ {input.long} == *.gz ]]; then
            gunzip -c {input.long} > {output.uncompressed_long}
        else
            cp {input.long} {output.uncompressed_long}
        fi

        # Set Java ENV variables to avoid a heap memory crash
        export _JAVA_OPTIONS="-Xmx120g -Xms4g"
        export JAVA_TOOL_OPTIONS="-Xmx120g -Xms4g"

        # Run OPERA-MS
        perl /operams/OPERA-MS.pl \
            --short-read1 {output.concatenated_forwards} \
            --short-read2 {output.concatenated_reverses} \
            --long-read {output.uncompressed_long} \
            --no-ref-clustering \
            --num-processors {params.cpus} \
            --out-dir {params.out_dir} \
        2> {log} 1>&2

        # Compress final assembly
        gzip -c {params.out_dir}/final_assembly.fasta > {output.fasta}
        """



rule operams:
    """Collect Opera-MS hybrid assemblies"""
    input:
        [OPERAMS_HYBRID / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
