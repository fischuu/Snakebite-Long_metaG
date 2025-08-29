# Generate more configuration keys
# Define the script_folder dynamically based on the pipeline_folder
SCRIPT_FOLDER = os.path.join(config["pipeline_folder"], "workflow", "scripts")

READS = Path("reads/")
WD = os.getcwd()

COMPARE = Path("results/metaquast/")
FLYE_LONG = Path("results/flye/")
MEDAKA = Path("results/medaka/")
OPERAMS_HYBRID = Path("results/operams/")
POLYPOLISH = Path("results/polypolish/")
RACON = Path("results/polish/")
SPADES_HYBRID = Path("results/spades/hybrid/")
SPADES_LONG = Path("results/spades/long_only/")
STRAINY = Path("results/strainy/")
