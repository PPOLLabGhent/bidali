from setuptools import setup, find_packages

setup(name = 'bidali',
      version = '0.1.11',
      description = 'Biological Data Analysis Library',
      url = 'https://github.com/beukueb/bidali',
      author = 'Christophe Van Neste',
      author_email = 'christophe.vanneste@ugent.be',
      license = 'MIT',
      packages = find_packages(),
      install_requires = [
          #Generated with `pipreqs .` and then moved here from requirements.txt
          'networkx',
          'biomart',
          'pandas',
          'gffutils',
          'scipy',
          'numpy',
          'pyliftover',
          'requests',
          'seaborn',
          'setuptools',
          'rpy2',
          'matplotlib',
          'plumbum',
          'lifelines'
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
