Descriptor-Based Microkinetic Analyais Package (DescMAP)
=============

The **Desc**\ riptor-Based **M**\ icrokinetic **A**\ nalysis **P**\ ackage (DescMAP) is a Python library developed by the Vlachos Research Group at the University of Delaware. This code was developed to automate feature selection and volcano curve generation for heterogeneous catalysis using empirical and semi-empirical approaches coupled with microkinetic modeling. Inputting data via spreadsheets and controlling program behavior via template files increases flexibility and supported capabilities.

.. image:: docs/DescMAP_logo.png
    :width: 400px
    :align: center

Developers
----------

-  Jonathan Lym (jlym@udel.edu)
-  Xue Zong (xzong@udel.edu)

Dependencies
------------

- Python3
- `pMuTT`_: Generates input files for OpenMKM
- `OpenMKM`_: Runs Microkinetic models
- `Scikit-Learn`_: Choose descriptors based on DFT data
- `Plotly`_: Plots interactive volcano curves

License
-------

This project is licensed under the MIT License - see the LICENSE.md file for
details.

Funding
-------

This material is based upon work supported by the Department of Energy's Office 
of Energy Efficient and Renewable Energy's Advanced Manufacturing Office under 
Award Number DE-EE0007888-9.5.

.. _`pMuTT`: https://vlachosgroup.github.io/pMuTT/
.. _`OpenMKM`: https://vlachosgroup.github.io/openmkm/
.. _`Scikit-Learn`: https://scikit-learn.org/stable/
.. _`Plotly`: https://plotly.com/
