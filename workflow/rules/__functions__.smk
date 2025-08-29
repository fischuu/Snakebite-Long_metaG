def get_sample_and_library_from_assembly_id(assembly_id):
    """Get all the sample and library ids for a given assembly_id"""
    samples_in_assembly = samples[samples.assembly_id == assembly_id][
        ["sample_id", "library_id"]
    ].values.tolist()
    return samples_in_assembly


def _get_reads_from_assembly_id(wildcards, end):
    """Get the pair-end files"""
    assert end in ["forward", "reverse"]
    end = 1 if end == "forward" else 2
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    return [
        READS / f"{sample_id}.{library_id}_{end}.fq.gz"
        for sample_id, library_id in samples_in_assembly
    ]

def get_forwards_from_assembly_id(wildcards):
    """Get the forward files for megahit"""
    return _get_reads_from_assembly_id(wildcards, end="forward")


def get_reverses_from_assembly_id(wildcards):
    """Get the forward files for megahit"""
    return _get_reads_from_assembly_id(wildcards, end="reverse")


def aggregate_forwards_for_spades(wildcards):
    """Put all the forwards together separated by a whitespace"""
    forwards = [
        str(forward_)
        for forward_ in _get_reads_from_assembly_id(wildcards, end="forward")
    ]
    return " ".join(forwards)


def aggregate_reverses_for_spades(wildcards):
    """Put all the reverses together separated by a whitespace"""
    reverses = [
        str(reverse_)
        for reverse_ in _get_reads_from_assembly_id(wildcards, end="reverse")
    ]
    return " ".join(reverses)
  
def _get_longreads_from_assembly_id(wildcards):
    """Get the long-read files for a given assembly_id"""
    assembly_id = wildcards.assembly_id
    # select only rows with non-missing longread_filename
    longreads = samples.loc[
        (samples.assembly_id == assembly_id) & samples.longread_filename.notna(),
        "longread_filename"
    ].tolist()
    return [Path(lr) for lr in longreads]


def get_longreads_from_assembly_id(wildcards):
    """Wrapper for input functions in rules"""
    return _get_longreads_from_assembly_id(wildcards)

