Descriptor-Based Microkinetic Analyais Package (DescMAP)
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

.. image:: https://github.com/VlachosGroup/DescriptorMap/blob/32b25da138e320cd3c4d38083b563f30e0d80d3f/docs/logos/descmap_logo.png
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

-  Jonathan Lym (jonathanlym@gmail.com)
-  Xue Zong (xzong@udel.edu)


Dependencies
------------

-  Python >= 3.9
-  `Numpy <http://www.numpy.org/>`__: Used for vector and matrix operations
-  `Pandas <https://pandas.pydata.org/>`__: Used to import data from Excel or CSV files
-  `Scipy <https://www.scipy.org/>`__: Used for curve fitting
-  `Scikit-Learn <https://scikit-learn.org/stable/>`__: Choose
   descriptors based on DFT data
-  `RDKit <https://www.rdkit.org/docs/Overview.html>`__: Used for 
   constructing feasible chemical structures required by pGradd
-  `Matplotlib <https://matplotlib.org/>`__: Used for generating plots
-  `Plotly <https://plotly.com/>`__: Plots interactive volcano curves
-  `Chart-studio <https://chart-studio.plotly.com/feed/#/>`__: Provide utilities 
   for interfacing with Plotly's Chart Studio service
-  `pMuTT <https://vlachosgroup.github.io/pMuTT/>`__: Generates input files for OpenMKM
-  `pGradd <https://github.com/VlachosGroup/PythonGroupAdditivity/>`__: Estimate 
   thermodynamic properties of molecules
-  `VUnits <https://vlachosgroup.github.io/vunits/>`__: Unit conversion and constants
-  `xlsxwriter <https://xlsxwriter.readthedocs.io/>`__: Create Excel xlsx files
-  `pyDOE <https://pythonhosted.org/pyDOE/>`__: Experimental design package to 
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
