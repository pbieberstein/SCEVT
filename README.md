
![SCEVT LOGO](https://github.com/pbieberstein/SCEVT/blob/master/documentation/scevt-logo.png)


SCEVT is a tool to easily visualize and analyze scaffolds during de-novo genome assembly.

SCEVT consists of two scripts:
* scaal.py
* scaphy.py

[Whats the difference?](https://github.com/pbieberstein/SCEVT/wiki/Overview) 






#### Progress:
* scaal.py script is DONE ![#c5f015](https://placehold.it/15/c5f015/000000?text=+)
* scaphy.py script is DONE ![#c5f015](https://placehold.it/15/c5f015/000000?text=+)
* Documentation is currently being worked on ![#f03c15](https://placehold.it/15/f03c15/000000?text=+)





#### Getting Started
**Getting the Files**
```bash
# Go to where you want to have this tool
cd path/to/Project/directory
git clone https://github.com/pbieberstein/SCEVT.git SCEVT
```

**Installing Python & Dependencies**

This script was developed in Python 2.7

The easiest way : Install anaconda for python 2.7 on your local machine and then install biopython via:

```bash
conda install biopython matplotlib==1.5.3 pandas
conda install --channel bioconda gffutils
```

Alternatively, if you want to stay organized we recommend you install miniconda and then create
a new virtual environment with the dependencies for this project.
(https://conda.io/docs/install/quick.html)
(Additional conda help: https://conda.io/docs/_downloads/conda-cheatsheet.pdf)


```bash
cd path/to/Project/directory
conda create --prefix ./scevt-env biopython matplotlib==1.5.3 pandas
# This creates a new environment with biopython and matplotlib installed inside the folder "scevt_env"
```

**It's important to use matplotlib 1.5.3 otherwise SCEVT will run very slowly

Now when you want to run SCEVT, you'll first have to activate this new python environment via:
```bash
source activate scevt-env/bin/activate
```

Now open up a new terminal window to update the PATHs and now you're ready to run scaal and scaphy


Then you can run the tools via
```bash
cd path/to/Project/directory/SCEVT
cd Scripts
python scaal.py
# or
python scaphy.py
```
