def _get_read_file(wildcards, end):
    assert end in ["forward_filename", "reverse_filename"]
    return samples[
        (samples["sample_id"] == wildcards.sample)
        & (samples["library_id"] == wildcards.library)
        & (samples["longread_filename"].isna()) 
    ][end].values[0]


def get_forward(wildcards):
    """Get the forward read for a given sample and library"""
    return _get_read_file(wildcards, end="forward_filename")


def get_reverse(wildcards):
    """Get the reverse read for a given sample and library"""
    return _get_read_file(wildcards, end="reverse_filename")
  

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

# fastp
def _get_adapter(wildcards, end):
    """Get the adapter of the en from a file"""
    assert end in ["forward", "reverse"]
    end = "forward_adapter" if end == "forward" else "reverse_adapter"
    adapter = samples[
        (samples["sample_id"] == wildcards.sample_id)
        & (samples["library_id"] == wildcards.library_id)
        & (samples["longread_filename"].isna()) 
    ][end].tolist()[0]
    return adapter


def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    return _get_adapter(wildcards, end="forward")


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    return _get_adapter(wildcards, end="reverse")

def get_fastp_reads_from_assembly_id(wildcards, end):
    """Get the file_end for megahit"""
    assert end in ["forward", "reverse"]
    end = 1 if end == "forward" else 2
    assembly_id = wildcards.assembly_id
    samples_in_assembly = get_sample_and_library_from_assembly_id(assembly_id)
    return [
            FASTP / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in samples_in_assembly
    ]

def get_forwards_from_assembly_id(wildcards):
    """Get the forward files for megahit"""
    return get_fastp_reads_from_assembly_id(wildcards, end="forward")


def get_reverses_from_assembly_id(wildcards):
    """Get the forward files for megahit"""
    return get_fastp_reads_from_assembly_id(wildcards, end="reverse")

def aggregate_forwards_for_bwa(wildcards):
    """Put all the forwards together separated by a comma"""
    forwards = [
        str(forward_)
        for forward_ in get_fastp_reads_from_assembly_id(wildcards, end="forward")
    ]
    return ",".join(forwards)


def aggregate_reverses_for_bwa(wildcards):
    """Put all the reverses together separated by a comma"""
    reverses = [
        str(reverse_)
        for reverse_ in get_fastp_reads_from_assembly_id(wildcards, end="reverse")
    ]
    return ",".join(reverses)


# helper to collect unique filtered SAMs for a given assembly (preserves order)
def filtered_sams_for_assembly(wildcards):
    paths = []
    for a, sample, library in ASSEMBLY_SAMPLE_LIBRARY:
        if a == wildcards.assembly_id:
            for read in (1, 2):
                # match the naming you showed: sample.library_read.filtered.sam
                p = POLYPOLISH / a / f"{sample}.{library}_{read}.filtered.sam"
                paths.append(str(p))
    # dedupe while preserving order
    seen = set()
    uniq = []
    for p in paths:
        if p not in seen:
            uniq.append(p)
            seen.add(p)
    return uniq
