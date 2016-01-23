# Notes on Project and Development Environment

Should use Bioconda channel to simplify install of tabix & bedtools.
See http://daler.github.io/pybedtools/main.html#quick-install-via-conda

Add this channels to our .condarc

```
conda config --add channels 'bioconda'
```

Python version
-------------

- Pybedtools  version 2.7 or greater (Python 3 is supported)
- PyVCF  2.6-3.4
- Primer3-py 2.7-3.5


To simplify install of bedtools dependencies and Jupyter
initiated Python 3 env with these, then

```
conda create -n py3markers  bedtools jupyter python=3
source activate py3markers
pip install primer3-py pybedtools
pip install --upgrade pip
conda install pyvcf
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
