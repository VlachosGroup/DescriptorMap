Automated Descriptor Selection, Volcano Curve Generation, and Active Site Determination Using the DescMAP Software
========================================================

The **Desc**\ riptor-Based **M**\ icrokinetic **A**\ nalysis
**P**\ ackage (DescMAP) is a Python library developed by the Vlachos
Research Group at the University of Delaware. This code was developed to
automate descriptor selection and volcano curve generation for
heterogeneous catalysis using empirical and semi-empirical approaches
coupled with microkinetic modeling. Both electronic and geometric
descriptors are supported. Inputting data via spreadsheets and
controlling program behavior via template files increases flexibility
and supported capabilities.

.. image:: https://raw.githubusercontent.com/VlachosGroup/DescriptorMap/master/docs/logos/descmap_logo.png
    :width: 400px

Documentation
-------------

See our `documentation page <https://descriptormap.readthedocs.io/en/latest/>`__ for 
docstrings and more details.


Getting Started
---------------

1. Install using pip (see documentation for full instructions)

.. code:: bash

   $ pip install descmap

2. Look at the provided examples


Developers
----------

-  Xue Zong (xzong@udel.edu)
-  Jonathan Lym (jonathanlym@gmail.com)


Dependencies
------------

-  Python >= 3.9
-  `Numpy <http://www.numpy.org/>`__ >= 1.24.2: Used for vector and matrix operations
-  `Pandas <https://pandas.pydata.org/>`__ >= 1.5.3: Used to import data from Excel or CSV files
-  `Scipy <https://www.scipy.org/>`__ >= 1.10.0: Used for curve fitting
-  `Scikit-Learn <https://scikit-learn.org/stable/>`__ >= 1.2.1: Choose
   descriptors based on DFT data
-  `RDKit <https://www.rdkit.org/docs/Overview.html>`__ >= 2022.9.4: Used for 
   constructing feasible chemical structures required by pGradd
-  `Matplotlib <https://matplotlib.org/>`__ >= 3.6.3: Used for generating plots
-  `Plotly <https://plotly.com/>`__ >= 5.13.0: Plots interactive volcano curves
-  `Chart-studio <https://chart-studio.plotly.com/feed/#/>`__ >= 1.1.0: Provide utilities 
   for interfacing with Plotly's Chart Studio service
-  `pMuTT <https://vlachosgroup.github.io/pMuTT/>`__ >= 1.3.2: Generates input files for OpenMKM
-  `pGradd <https://github.com/VlachosGroup/PythonGroupAdditivity/>`__ >= 2.9.5: Estimate 
   thermodynamic properties of molecules
-  `VUnits <https://vlachosgroup.github.io/vunits/>`__ >= 0.0.4: Unit conversion and constants
-  `xlsxwriter <https://xlsxwriter.readthedocs.io/>`__ >= 3.0.8: Create Excel xlsx files
-  `pyDOE <https://pythonhosted.org/pyDOE/>`__ >= 0.3.8: Experimental design package to 
   provide sampling method


License
-------

This project is licensed under the MIT License - see the LICENSE.md file
for details.


Contributing
------------

If you have any suggestion or find a bug, please post to our Issues page
on GitHub.


Questions
---------

If you have any question or run into any issue, please post to our
Issues page on GitHub.


Funding
-------

This material is based upon work supported by the Department of Energy's
Office of Energy Efficient and Renewable Energy's Advanced Manufacturing
Office under Award Number DE-EE0007888-9.5.
