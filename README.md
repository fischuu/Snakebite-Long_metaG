# Snakebite-Long_metaG

This Snakemake pipeline dedicated to Metagenomic data analysis with long and short reads consists out of several modules that cover
a) full hybrid b) long reads with short read polishing and c) existing workflows (opera-ms)

![alt text](https://github.com/fischuu/Snakebite-Long_metaG/blob/main/flowchart/flowchart.png?raw=true)


# Requirements
The pipeline requires version 9. Further, it is currently tested with the slurm executor and as such this one is also required to be installed. 

In essence, it can also run with smaller Snakemake versions, but it might cause some troubles with the slurm executor and the escalation part
of the resource allocation. 

Supports:
SLURM executor / local execution
conda environment (not tested)
docker/singularity/apptainer support

## Python dependencies
Since Snakemake 8, it is required to install an executor plugin to submit jobs to a queueing system of a HPC system. Please ensure you have installed the corresponding Snakemake plugins installed in case you want to submit your jobs to a queueing system

```
pip install snakemake-executor-plugin-cluster-generic
pip install snakemake-executor-plugin-slurm

```

Of course you can also install other, specific executor plugins, but this might need more adjustments to the existing files.

Depending on your Python version, you need to install a Pandas version > 2.1, there were errors when newer Python versions met older Pandas version. In case you run into obscure Pandas error, please make sure to install a newer pandas, e.g.

```
pip install pandas==2.2.3
```

# Installation

You can install the pipeline by cloning this repository

The recommended setup is to have a dedicated pipeline folder (the cloned repository), that carries the functionality and which should not require any changes. 

Then the project should have somewhere an own folder and the required configuration files are copied to it. The steps to perform are

```
# Go to the folder, to where you would like to clone the pipeline, e.g. 
 cd /users/fischerd/git

# First, clone the pipeline into that folder
  git clone git@github.com:fischuu/Snakebite-Long_metaG.git

# In case the previous steps fails with an error that contains
# git@github.com: Permission denied (publickey)
# it indicates that you do not have a ssh key exchanged with GitHub and you could clone the repository then instead like this
# git clone https://github.com/fischuu/Snakebite-Long_metaG.git
  
# Setting ENV variable to get downstream code more generic (so, this is the directory to where you cloned the pipeline)
  cd Snakebite-Long_metaG
  PIPELINEFOLDER=$(pwd)
  
# If previous doesn't work, you can set it also manually like for example this
  PIPELINEFOLDER="/users/fischerd/git/Snakebite-Long_metaG"
```

Next, we setup a project folder in our scratch space of the HPC, here we will run the pipeline

```
# Go to the project space of your HPC, e.g. 
  cd /scratch/project_2009831

# Create a folder for the new project
  mkdir My_metag_project
  cd My_metag_project   
  
# For convenience, we set again a ENV variable, so that the code will be more generic
  PROJECTFOLDER=$(pwd)
  
# Or manually the same thing:  
  PROJECTFOLDER="/scratch/project_2009831/My_metag_project"
```
