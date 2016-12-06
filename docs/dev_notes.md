# Notes on Project and Development Environment

** Go with Python >= 3.4**

 Install miniconda3 http://conda.pydata.org/miniconda.html

Worth taking a look at http://www-huber.embl.de/HTSeq/doc/index.html

http://dx.doi.org/10.1093/bioinformatics/btu638

Should use Bioconda channel to simplify install of tabix & bedtools.
See http://daler.github.io/pybedtools/main.html#quick-install-via-conda

Add this channel to our .condarc and upgrade pip to stop nags

```
conda config --add channels 'bioconda'
pip install --upgrade pip
```

Python version
-------------

- Pybedtools  version 2.7 or greater (Python 3 is supported)
- PyVCF  2.6-3.4
- Primer3-py 2.7-3.5


To simplify install of bedtools dependencies and Jupyter
initiated Python 3 env with these, then

```
conda create -y -n py3markers  bedtools jupyter pyvcf
source activate py3markers
pip install primer3-py pybedtools

```
Listed packages to export file
---------

```
conda list --export >  conda-package-list.txt
```

To regenerate Environment

```
conda create -n py3markers --file conda-package-list.txt
```

Fasta Access
------------

- could use https://github.com/mdshw5/pyfaidx  for cutting out slices

umelt
=====

- service at https://www.dna.utah.edu/db/services/cgi-bin/udesign.cgi
- See web app at https://www.dna.utah.edu/umelt/um.php

- Might be worth a look at http://biopython.org/DIST/docs/api/Bio.SeqUtils.MeltingTemp-module.html
