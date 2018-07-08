# BIological Data Analysis Library - Documentation
<img title="bidali logo" src="bidali_logo.svg" width="200">

## Requirements

- Install python3 from https://www.python.org/downloads/ (version 3.6.2)
- Install git from https://git-scm.com/downloads
- https://virtualenvwrapper.readthedocs.io/en/latest/

### netCDF4 output

- sudo apt-get install libnetcdf-dev libhdf5-dev

## Installation

### Mac OS X

Open `Terminal` and copy paste below line by line:

     mkdir -p ~/{repos,LSData/cache} && cd ~/repos
     git clone https://github.com/dicaso/bidali.git
     BIDADIR=~/repos/bidali
     mkvirtualenv -a $BIDADIR -i ipython -i rpy2 -i twine -i Sphinx \
                  -r $BIDADIR/requirements.txt bidali
     python setup.py test # runs all the tests
     deactivate

The bidali package can now be used. The following example illustrates preparing
the neuroblastoma 39 cell line data ready for use in R:

    . ~/repos/bidali/velsp/bin/activate
    python -m LSD --list LSD.dealer.external.celllines.get_NB39

To update the package, and get new available datasets:

   cd ~/repos/bidali
   git pull origin master
   . velsp/bin/activate

