===============================
MUSiCC API via the command line
===============================
The MUSiCC module handles all calculations internally.
MUSiCC offers an interface to the MUSiCC functionality via the command line and the run_musicc script.

Usage:
------

``run_musicc.py input_file [options]``

Required arguments:
-------------------

**input_file**
    Input abundance file to correct

Optional arguments:
-------------------

**-h, --help**
    show help message and exit

**-o OUTPUT_FILE, --out OUTPUT_FILE**
    Output destination for corrected abundance (default: MUSiCC.tab)

**-if {tab,csv}, --input_format {tab,csv}**
    Option indicating the format of the input file (default: tab)

**-of {tab,csv}, --output_format {tab,csv}**
    Option indicating the format of the output file (default: tab)

**-n, --normalize**
    Apply MUSiCC normalization (default: false)

**-c {use_generic, learn_model}, --correct {use_generic,learn_model}**
    Correct abundance per-sample using MUSiCC (default: false)

**-perf, --performance**
    Calculate model performance on various gene sets (may add to running time) (default: false)

**-v, --verbose**
    Increase verbosity of module (default: false)


============================
MUSiCC API via python script
============================
MUSiCC can also be used directly inside a python script. Passing variables and flags to the MUSiCC script is done by
creating a dictionary and passing it to the function *correct_and_normalize*, as shown below.

Usage:
------

>>> from musicc.core import correct_and_normalize
>>> musicc_args = {'input_file': 'test_musicc/lib/python3.3/site-packages/musicc/examples/simulated_ko_relative_abundance.tab', 'output_file': 'MUSiCC.tab','input_format': 'tab', 'output_format': 'tab', 'musicc_inter': True, 'musicc_intra': 'learn_model','compute_scores': True, 'verbose': True}
>>> correct_and_normalize(musicc_args)

Required arguments:
-------------------

**input_file**
    Input abundance file to correct

Optional arguments:
-------------------

**output_file**
    Output destination for corrected abundance (default: MUSiCC.tab)

**input_format {'tab','csv'}**
    Option indicating the format of the input file (default: 'tab')

**output_format {'tab','csv'}**
    Option indicating the format of the output file (default: 'tab')

**musicc_inter {True, False}**
    Apply MUSiCC normalization (default: False)

**musicc_intra {'use_generic', 'learn_model', 'None'}**
    Correct abundance per-sample using MUSiCC (default: 'None')

**compute_scores {True, False}**
    Calculate model performance on various gene sets (may add to running time) (default: False)

**verbose {True, False}**
    Increase verbosity of module (default: False)

========
Examples
========
In the *musicc/examples* directory, the file *simulated_ko_relative_abundance.tab* contains simulated KO abundance measurements of 20 samples described in the
MUSiCC manuscript. Using this file as input for MUSiCC results in the following files:

- simulated_ko_MUSiCC_Normalized.tab (only normalization)
- simulated_ko_MUSiCC_Normalized_Corrected_use_generic.tab (normalize and correct using the generic model learned from HMP)
- simulated_ko_MUSiCC_Normalized_Corrected_learn_model.tab (normalize and correct learning a new model for each sample)

The commands used were the following (via command line):

``run_musicc.py musicc/examples/simulated_ko_relative_abundance.tab -n -perf -v -o musicc/examples/simulated_ko_MUSiCC_Normalized.tab``

``run_musicc.py musicc/examples/simulated_ko_relative_abundance.tab -n -c use_generic -perf -v -o musicc/examples/simulated_ko_MUSiCC_Normalized_Corrected_use_generic.tab``

``run_musicc.py musicc/examples/simulated_ko_relative_abundance.tab -n -c learn_model -perf -v -o musicc/examples/simulated_ko_MUSiCC_Normalized_Corrected_learn_model.tab``
