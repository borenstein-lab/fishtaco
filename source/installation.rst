Installation Instructions
=========================
.. index:: Installation

Prerequisites for installation
------------------------------

In order for FishTaco to run successfully, the following Python modules should be pre-installed on your system:

- Numpy >= 1.6.1 (`<http://www.numpy.org/>`_)
- Scipy >= 0.9 (`<http://www.scipy.org/>`_)
- Scikit-learn >= 0.15.2 (`<http://scikit-learn.org/stable/>`_)
- Pandas >= 0.14 (`<http://pandas.pydata.org/>`_)
- StatsModels >= 0.5.0 (`<http://statsmodels.sourceforge.net/>`_)
- MUSiCC >= 1.0.1 (`<https://github.com/omanor/MUSiCC>`_)

If you have *pip* installed, you can install these packages by running the following command:

``pip install -U numpy scipy scikit-learn pandas statsmodels MUSiCC``

Installing FishTaco
-------------------

To install FishTaco, download the package from `GitHub <https://github.com/omanor/fishtaco/archive/v1.0.1.tar.gz>`_.

After downloading FishTaco, you’ll need to unzip the file. If you’ve downloaded the release version, do this with the following command:

``tar -xzf fishtaco-1.0.tar.gz``

You’ll then change into the new FishTaco directory as follows:

``cd fishtaco-1.0``

and install using the following command:

``python setup.py install``

ALTERNATIVELY, you can install FishTaco directly from PyPI by running:

``pip install -U fishtaco``

Note for windows users: Under some windows installations, Scipy may fail when importing the Stats module. Workarounds may be found online, such
as `here <https://code.google.com/p/pythonxy/issues/detail?id=745>`_.

Testing the software package
----------------------------

After downloading and installing the software, we recommend testing it by running the following command:

``test_fishtaco.py``

This will invoke a series of tests. A correct output should end with:

``Ran 3 tests in X.XXXXs``

``OK``
