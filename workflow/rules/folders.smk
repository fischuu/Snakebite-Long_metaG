# Generate more configuration keys
# Define the script_folder dynamically based on the pipeline_folder
SCRIPT_FOLDER = os.path.join(config["pipeline_folder"], "workflow", "scripts")

READS = Path("reads/")
WD = os.getcwd()

BWA = Path("results/bwa/")
BWA_INDEX = BWA / "index/"
BWA_PAIRED = BWA / "paired/"
BWA_SINGLE = BWA / "single/"
COMPARE = Path("results/metaquast/")
FASTP = Path("results/fastp")
FLYE_LONG = Path("results/flye/")
SEQTK = Path("results/seqtk/")
MEDAKA = Path("results/medaka/")
PILON = Path("results/pilon/")
OPERAMS_HYBRID = Path("results/operams/")
POLYPOLISH = Path("results/polypolish/")
RACON = Path("results/racon/")
SPADES_HYBRID = Path("results/spades/")
STRAINY = Path("results/strainy/")
