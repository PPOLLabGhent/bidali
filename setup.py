from setuptools import setup

setup(name = 'lospy',
      version = '0.1',
      description = 'Lab-related occult and secretive python package',
      url = 'https://github.ugent.be/cvneste/lospy',
      author = 'Christophe Van Neste',
      author_email = 'christophe.vanneste@ugent.be',
      license = 'MIT',
      packages = ['lospy','LSD'],
      install_requires = [
          'gseapy',
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
