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

.. image:: docs/logos/descmap_logo.png
    :width: 400px

Installation
------------

.. code:: bash

   $ pip install descmap

Developers
----------

-  Jonathan Lym (jonathanlym@gmail.com)
-  Xue Zong (xzong@udel.edu)

Dependencies
------------

-  Python3
-  `pMuTT <https://vlachosgroup.github.io/pMuTT/>`__: Generates input
   files for OpenMKM
-  `OpenMKM <https://vlachosgroup.github.io/openmkm/>`__: Runs
   Microkinetic models
-  `Scikit-Learn <https://scikit-learn.org/stable/>`__: Choose
   descriptors based on DFT data
-  `Plotly <https://plotly.com/>`__: Plots interactive volcano curves

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
