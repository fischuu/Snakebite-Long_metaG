rule reads__link_run:
    """Make a link to the original file, with a prettier name than default"""
    input:
        forward_=get_forward,
        reverse_=get_reverse,
    output:
        forward_=READS / "{sample}.{library}_1.fq.gz",
        reverse_=READS / "{sample}.{library}_2.fq.gz",
    log:
        READS / "{sample}.{library}.log"
    benchmark:
        READS / "benchmark/{sample}.{library}.tsv"
    container:
        docker["reads"]
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_} 2>  {log} 1>&2
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_} 2>> {log} 1>&2
        """

rule reads__link:
    input:
        expand(READS / "{sample}.{library}_{end}.fq.gz",
               sample=[s for s,_ in SAMPLE_LIBRARY],
               library=[l for _,l in SAMPLE_LIBRARY],
               end=["1","2"])
