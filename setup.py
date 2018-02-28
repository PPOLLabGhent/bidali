from setuptools import setup

setup(name = 'bidali',
      version = '0.1.0',
      description = 'Biological Data Analysis Library',
      url = 'https://github.com/beukueb/bidali',
      author = 'Christophe Van Neste',
      author_email = 'christophe.vanneste@ugent.be',
      license = 'MIT',
      packages = ['bidali','LSD'],
      install_requires = [
          #Generated with `pipreqs .` and then moved here from requirements.txt
          'gffutils==0.8.7.1',
          'scipy==0.19.0',
          'lifelines==0.9.4',
          'seaborn==0.7.1',
          'requests==2.13.0',
          'pyliftover==0.3',
          'networkx==1.11',
          'rpy2==2.8.5',
          'matplotlib==2.0.0',
          'pandas==0.19.2',
          'setuptools==32.2.0',
          'numpy==1.12.0',
      ],
      zip_safe = False,
      #entry_points = {
      #    'console_scripts': ['getLSDataset=LSD.command_line:main'],
      #},
      test_suite = 'nose.collector',
      tests_require = ['nose']
)

#To install with symlink, so that changes are immediately available:
#pip install -e . 
