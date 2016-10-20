
import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


def read(*paths):
    """Build a file path from *paths* and return the contents."""
    with open(os.path.join(*paths), 'r') as f:
        return f.read()

setup(name='FishTaco',
      version='1.1.1',
      classifiers=['License :: Free for non-commercial use'],
      description='FishTaco: a metagenomic computational framework, aiming to '
                  'identify the taxa that are driving functional shifts in microbiomes.',
      long_description=('\n\n' + ('README.rst') + '\n\n' + read('HISTORY.rst') +
                        '\n\n' + read('AUTHORS.rst') + '\n\n' + read('LICENSE') + '\n\n'),
      author='Ohad Manor',
      author_email='omanor@gmail.com',
      url='http://omanor.github.io/fishtaco/',
      packages=['fishtaco'],
      package_data={'fishtaco': ['data/*.tab', 'data/*.lst',
                                 'examples/*.tab', 'examples/output/*.tab']},
      install_requires=['NumPy >= 1.6.1', 'SciPy >= 0.9', 'scikit-learn >= 0.15.2',
                        'statsmodels >= 0.5.0', 'pandas >= 0.14', 'MUSiCC >= 1.0.1'],
      scripts=['scripts/run_fishtaco.py', 'tests/test_fishtaco.py'],
      )




