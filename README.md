# Documentation
## Installation
### Mac OS X

- Install python3 from https://www.python.org/downloads/ (version 3.6.2)
- Install git from https://git-scm.com/downloads
- open `Terminal`

Copy paste below line by line:

     sudo pip3 install virtualenv # asks your password
     mkdir -p ~/{repos,Data/LSDpy} && cd ~/repos
     git clone https://github.ugent.be/cvneste/lospy.git # asks username/password
     cd ~/repos/lospy
     virtualenv velsp
     . velsp/bin/activate
     pip install .
     deactivate

The lospy package can now be used. The following example illustrates preparing
the neuroblastoma 39 cell line data ready for use in R:

    . ~/repos/lospy/velsp/bin/activate
    python -m LSD --list LSD.dealer.external.celllines.get_NB39

To update the package, and get new available datasets:

   cd ~/repos/lospy
   git pull origin master
   . velsp/bin/activate

