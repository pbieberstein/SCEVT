
![SCEVT LOGO](https://github.com/pbieberstein/SCEVT/blob/master/documentation/scevt-logo.png)


SCEVT is a tool to easily visualize and analyze scaffolds during de-novo genome assembly.

SCEVT consists of two scripts:
* scaal.py
* scaphy.py


**scaphy.py** (Scaffold to Physical Reference Mapping)
-- 
scaphy is a tool to visualize scaffolds in relation to a reference genome assembly. Specifically, it draws gaps within the scaffolds (esspecially helpful for BioNano assisted scaffolds) and draws mappings to a reference chromosome whenever the genes match. It also highlights when a gene is on the scaffold that is not on the specified chromosome on the reference genome (meaning you have probably anchored a new contig).

Here is an example output:
![example scaphy output](https://github.com/pbieberstein/SCEVT/blob/master/documentation/scaphy_example_output.png)

[How to use](https://github.com/pbieberstein/SCEVT/wiki/scaphy.py-User-Guide)


**scaco.py** (Scaffold Comparison) 
--
[How it Works]
scaco directly compares two scaffolds based on gene annotations. It highlights and maps which genes are similar on the two scaffolds, and also highlights which genes are present on one but not the other. Additionally, it also plots the gaps within the scaffolds. This is useful for comparing haplotype contigs of a de-novo assembly.

Here is an example output:
![example scaco output](https://github.com/pbieberstein/SCEVT/blob/master/documentation/scaal_example_output.png)


[How to use](https://github.com/pbieberstein/SCEVT/wiki/scaco.py-User-Guide)



**Installation**
--
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
conda create --prefix ./scevt-env biopython matplotlib==1.5.3 pandas biopython
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



**Progress:**
--
* scaal.py script is DONE ![#c5f015](https://placehold.it/15/c5f015/000000?text=+)
* scaphy.py script is DONE ![#c5f015](https://placehold.it/15/c5f015/000000?text=+)
* Documentation is DONE ![#c5f015](https://placehold.it/15/c5f015/000000?text=+)

*This tool was written to assist in a de-novo genome assembly project at ETH-Zurich*

It is not activily maintained but it should still be useful. If you have any questions/ideas/concerns, contact me.
