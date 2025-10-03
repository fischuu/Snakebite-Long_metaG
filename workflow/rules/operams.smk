rule operams__run:
    """Run Opera-MS hybrid assembly (long + short reads)"""
    input:
        long=get_longreads_from_assembly_id,
        forwards=get_forwards_from_assembly_id,
        reverses=get_reverses_from_assembly_id,
    output:
        fasta=OPERAMS_HYBRID / "{assembly_id}.fa.gz",
        concatenated_forwards=temp(OPERAMS_HYBRID / "{assembly_id}_R1_concat.fastq"),
        concatenated_reverses=temp(OPERAMS_HYBRID / "{assembly_id}_R2_concat.fastq"),
        uncompressed_long=temp(OPERAMS_HYBRID / "{assembly_id}_longreads.fastq"),
    log:
        log=OPERAMS_HYBRID / "{assembly_id}.log",
    container:
        docker["operams"]
    threads: esc("cpus", "operams__run")
    resources:
        runtime=esc("runtime", "operams__run"),
        mem_mb=esc("mem_mb", "operams__run"),
        cpus_per_task=esc("cpus", "operams__run"),
        slurm_partition=esc("partition", "operams__run"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'operams__run')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("operams__run"))
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
        export _JAVA_OPTIONS="-Xmx360g -Xms4g"
        export JAVA_TOOL_OPTIONS="-Xmx360g -Xms4g"

        # Run OPERA-MS
        perl /operams/OPERA-MS.pl \
            --short-read1 {output.concatenated_forwards} \
            --short-read2 {output.concatenated_reverses} \
            --long-read {output.uncompressed_long} \
            --no-ref-clustering \
            --num-processors {threads} \
            --out-dir {params.out_dir} \
        2> {log} 1>&2

        # Compress final assembly
        gzip -c {params.out_dir}/contigs.polished.fasta > {output.fasta}
        """



rule operams:
    """Collect Opera-MS hybrid assemblies"""
    input:
        [OPERAMS_HYBRID / f"{assembly_id}.fa.gz" for assembly_id in ASSEMBLIES],
