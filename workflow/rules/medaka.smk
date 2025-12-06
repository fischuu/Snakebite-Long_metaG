rule medaka__mini_align:
    """Medaka mini align step"""
    input:
        fa=RACON / "{assembly_id}.racon.fa.gz",
        long=get_longreads_from_assembly_id,
    output:
        bam=MEDAKA / "{assembly_id}.calls_to_draft.bam",
    log:
        MEDAKA / "{assembly_id}.medaka_mini_align.log",
    container:
        docker["medaka"]
    threads: esc("cpus", "medaka__mini_align")
    resources:
        runtime=esc("runtime", "medaka__mini_align"),
        mem_mb=esc("mem_mb", "medaka__mini_align"),
        cpus_per_task=esc("cpus", "medaka__mini_align"),
        slurm_partition=esc("partition", "medaka__mini_align"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'medaka__mini_align')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("medaka__mini_align"))
    params:
        out=lambda w: MEDAKA / f"{w.assembly_id}.calls_to_draft",
        model=params["assemble"]["medaka"]["model"],
    shell:
        r"""
          mini_align -i {input.long} \
                     -r {input.fa} -P -m \
                     -p {params.out} -t {threads} \
                     2> {log} 1>&2
        """

rule extract_regions:
    input:
        fa = RACON / "{assembly_id}.racon.fa.gz"
    output:
        regions = MEDAKA / "{assembly_id}/regions.txt"
    shell:
        r"""
        gunzip -c {input.fa} | grep '^>' | cut -d' ' -f1 | sed 's/>//' > {output.regions}
        """



checkpoint split_regions:
    input:
        regions = MEDAKA / "{assembly_id}/regions.txt"
    output:
        directory(MEDAKA / "{assembly_id}/chunks")
    params:
        size = params["assemble"]["medaka"]["group_size"]
    shell:
        r"""
        mkdir -p {output}
        split -l {params.size} {input.regions} {output}/chunk_
        """

#def get_medaka_chunks(wildcards):
#    checkpoint_output = checkpoints.split_regions.get(**wildcards).output[0]
#    # Get the chunk indices (numbers) dynamically
#    chunk_files = glob_wildcards(os.path.join(checkpoint_output, "chunk_{i}")).i
#    return chunk_files

def get_medaka_chunks(wildcards):
    checkpoint_output = checkpoints.split_regions.get(**wildcards).output[0]
    # Only match chunk files without any extension
    chunk_files = [f for f in glob_wildcards(os.path.join(checkpoint_output, "chunk_{i}")).i
                   if not f.endswith((".hdf", ".log"))]
    return chunk_files


rule medaka_inference_chunk:
    input:
        bam = MEDAKA / "{assembly_id}.calls_to_draft.bam",
        chunk = MEDAKA / "{assembly_id}/chunks/chunk_{chunk}",
        regions = MEDAKA / "{assembly_id}/regions.txt"
    output:
        hdf = MEDAKA / "{assembly_id}/chunks/chunk_{chunk}.hdf"
    log:
        MEDAKA / "{assembly_id}/chunks/{assembly_id}_chunk_{chunk}.medaka_inference_chunk.log",
    container:
        docker["medaka"]
    threads: 1
    resources:
        runtime=esc("runtime", "medaka_inference_chunk"),
        mem_mb=esc("mem_mb", "medaka_inference_chunk"),
        cpus_per_task=esc("cpus", "medaka_inference_chunk"),
        slurm_partition=esc("partition", "medaka_inference_chunk"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'medaka_inference_chunk')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("medaka_inference_chunk"))
    shell:
        r"""
        regions=$(tr '\n' ' ' < {input.chunk})
        medaka inference \
            {input.bam} \
            {output.hdf} \
            --region $regions \
            2> {log}
        """


#def get_medaka_chunk_hdfs(wildcards):
#    # Wait for split_regions checkpoint
#    chunks = checkpoints.split_regions.get(**wildcards).output[0]
#    chunk_ids = glob_wildcards(os.path.join(chunks, "chunk_{i}")).i
#    return expand(MEDAKA / "{assembly_id}/chunks/chunk_{chunk}.hdf",
#                  assembly_id=wildcards.assembly_id,
#                  chunk=chunk_ids)

def get_medaka_chunk_hdfs(wildcards):
    # Wait for split_regions checkpoint
    chunks_dir = checkpoints.split_regions.get(**wildcards).output[0]
    chunk_ids = [f for f in glob_wildcards(os.path.join(chunks_dir, "chunk_{i}")).i
                 if not f.endswith((".hdf", ".log"))]
    return expand(MEDAKA / "{assembly_id}/chunks/chunk_{chunk}.hdf",
                  assembly_id=wildcards.assembly_id,
                  chunk=chunk_ids)



rule medaka_sequence:
    input:
        hdfs = get_medaka_chunk_hdfs,
        chunks = directory(MEDAKA / "{assembly_id}/chunks"),
        bam=MEDAKA / "{assembly_id}.calls_to_draft.bam",
    output:
        fasta = MEDAKA / "{assembly_id}.polished.fasta"
    log:
        MEDAKA / "{assembly_id}.medaka_sequence.log",
    container:
        docker["medaka"]
    threads: 1
    resources:
        runtime=esc("runtime", "medaka_sequence"),
        mem_mb=esc("mem_mb", "medaka_sequence"),
        cpus_per_task=esc("cpus", "medaka_sequence"),
        slurm_partition=esc("partition", "medaka_sequence"),
        gres=lambda wc, attempt: f"{get_resources(wc, attempt, 'medaka_sequence')['nvme']}",
        attempt=get_attempt,
    retries: len(get_escalation_order("medaka_sequence"))
    shell:
        r"""
        medaka sequence {input.hdfs} {output.fasta} \
                     2> {log} 1>&2
        """


rule medaka:
    """Collect all Medaka results"""
    input:
        [MEDAKA / f"{assembly_id}.polished.fasta" for assembly_id in ASSEMBLIES],
