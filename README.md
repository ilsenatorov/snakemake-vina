# Vina Snakemake

## Installation

1. Clone the repository with
2. Create a conda environment with `conda env create -f workflow/envs/main.yml`
3. Activate the environment with `conda activate vina`
4. Edit the `config/config.yaml` file with your settings (see [configuration](##configuration))
5. Run the code with `snakemake --use-conda -j <num_cores>`

## Configuration

In the config file, the two most important parameters are `resources` and `results` - they define where the original files are located (receptors and ligands), and where to store the results.

The recommended structure of the projects is the following:
```
datasets
├── dataset1
│   └── resources
│       ├── ligands
│       └── receptors
└── dataset2
    └── resources
        ├── ligands
        └── receptors
```

Then if you want to run the workflow for dataset1, the `resources` entry in the config should be `datasets/dataset1/resources` and `results` should be `datasets/dataset1/results`.

The `results` folder will then be filled with the results of the workflow.