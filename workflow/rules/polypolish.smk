rule polypolish_filter:
    """Filter SAM alignments for Polypolish"""
    input:
        sam_1=BWA_SINGLE / "{assembly_id}" / "{sample_id}.{library_id}_1.sam",
        sam_2=BWA_SINGLE / "{assembly_id}" / "{sample_id}.{library_id}_2.sam",
    output:
        filt_1=POLYPOLISH / "{assembly_id}" / "{sample_id}.{library_id}_1.filtered.sam",
        filt_2=POLYPOLISH / "{assembly_id}" / "{sample_id}.{library_id}_2.filtered.sam",
    log:
        POLYPOLISH / "{assembly_id}" / "{sample_id}.{library_id}.filter.log",
    container:
        docker["polypolish"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    shell:
        r"""
        # mkdir -p $(dirname {output.filt_1})
        polypolish filter \
            --in1 {input.sam_1} --in2 {input.sam_2} \
            --out1 {output.filt_1} --out2 {output.filt_2} \
            2>> {log}
        """


rule polypolish_run:
    """Run Polypolish polishing on all filtered SAMs for an assembly"""
    input:
        asm=MEDAKA / "{assembly_id}.medaka.fa.gz",
        sams = filtered_sams_for_assembly,
    output:
        fa=POLYPOLISH / "{assembly_id}.polypolish.fa.gz",
    log:
        POLYPOLISH / "{assembly_id}.polypolish.log",
    container:
        docker["polypolish"]
    threads: config["resources"]["cpu_per_task"]["multi_thread"]
    params:
        tmp=POLYPOLISH / "{assembly_id}.tmp"
    shell:
        r"""
        polypolish polish {input.asm} {input.sams} \
            | gzip -c > {output.fa} 2>> {log}
        """


rule polypolish:
    """Collect all polished assemblies"""
    input:
        [POLYPOLISH / f"{assembly_id}.polypolish.fa.gz" for assembly_id in ASSEMBLIES],
