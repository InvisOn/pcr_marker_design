---
title: 'pcr-marker-design'
...

[![Latest PyPI version](https://img.shields.io/pypi/v/pcr-marker-design.svg)](https://pypi.python.org/pypi/pcr-marker-design)



# Tools for PCR assay design from NGS variant data

Galaxy-free refactoring of https://github.com/cfljam/galaxy-pcr-markers

Usage
=====

Installation
============

- create and activate  a clean Python 3xx env

e.g.
```

conda create -y -n Py3PCRtest python=3.5  numpy cython bedtools
source activate  Py3PCRtest
```


- optionally append jupyter notebook (or just ipykernel ) for usage in notebooks
```
conda install  jupyter notebook
```

- pull this repo

```
git clone https://github.com/PlantandFoodResearch/pcr_marker_design
```
- pip install it
```
pip install ./pcr_marker_design
```

## For development work
```
pip install -e ./pcr_marker_design
```
## Create an identical env

- in OSX
```
conda create --name <env> --file ./docs/Py3_environment.spec.txt
```

- other (may need some editing)
```
conda env create -f ./docs/Py3_environment.yml
```

Requirements
------------

Python 3.5

To Run Tests
-----------

from root directory
>pytest

Melt Prediction
---------------

- Relies on Univ of Utah Wittwer Lab Service
- See a test by visting https://www.dna.utah.edu/db/services/cgi-bin/udesign.cgi?seq=CTGATCGATCGTACGGCGCATCGTAGCTCWTAGCTACGCGCGTAGCTAGCTGCCGTAGC&rs=0&cation=20&mg=2&dmso=0


Compatibility
=============

Licence
=======

Authors
=======

pcr-marker-design was written by [John
McCallum](john.mccallum@plantandfood.co.nz) ,
[Susan Thomson](susan.thomson@plantandfood.co.nz) [Leshi Chen](), [Scout Liu]() [Brett Davis]()
