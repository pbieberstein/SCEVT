
![SCEVT LOGO](https://raw.githubusercontent.com/pbieberstein/SCEVT/master/scevt-logo.png)


SCEVT is a tool to easily visualize and analyze scaffolds during de-novo genome assembly.

Here is an sample output from the scaff_to_scaff_vis.py script which takes two scaffolds sequences and BLAT gene mappings,
and then visualizes the gaps in the scaffolds (in red), the unique genes on the scaffolds (in green) and the genes that are present
on both scaffolds (in purple)

Eventually, you can expect plots like this:
![sample_scaffold_vis](https://github.com/pbieberstein/SCEVT/blob/master/sample_output.png)


#### Currently still in development:
Progress:
* scaff_to_scaff_vis.py script is DONE ![#c5f015](https://placehold.it/15/c5f015/000000?text=+)
* gene_vis.py script is currently being worked on ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)
* snp_vis.py has not been started yet ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)

#### Getting Started

**Installing Python Dependencies**
This script was developed in Python 2.7

Use miniconda to create a new Python environment and install BioPython and matplotlib
(https://conda.io/docs/install/quick.html)
```bash
cd path/to/project
conda create --prefix ./scevt-env biopython matplotlib
# This creates a new environment with biopython and matplotlib installed inside the folder "scevt_env"
```
(Additional conda help: https://conda.io/docs/_downloads/conda-cheatsheet.pdf)

Once this is completed, you can run the script! Just be sure to edit the paths inside the scripts
for your relevant files (for now the specifics are inside the script headers)