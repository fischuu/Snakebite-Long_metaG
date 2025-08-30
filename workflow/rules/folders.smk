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
RACON = Path("results/racon/")
SPADES_HYBRID = Path("results/spades/")
STRAINY = Path("results/strainy/")
