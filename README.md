# Descriptor-Based Microkinetic Analyais Package (DescMAP)

The **Desc**riptor-Based **M**icrokinetic **A**nalysis **P**ackage
(DescMAP) is a Python library developed by the Vlachos Research Group at
the University of Delaware. This code was developed to automate feature
selection and volcano curve generation for heterogeneous catalysis using
empirical and semi-empirical approaches coupled with microkinetic
modeling. Inputting data via spreadsheets and controlling program
behavior via template files increases flexibility and supported
capabilities.

![image](docs/DescMAP_logo.png)

## Installation

```bash
$ pip install descmap
```

## Developers

-   Jonathan Lym (<jlym@udel.edu>)
-   Xue Zong (<xzong@udel.edu>)

## Dependencies

-   Python3
-   [pMuTT](https://vlachosgroup.github.io/pMuTT/): Generates input
    files for OpenMKM
-   [OpenMKM](https://vlachosgroup.github.io/openmkm/): Runs
    Microkinetic models
-   [Scikit-Learn](https://scikit-learn.org/stable/): Choose descriptors
    based on DFT data
-   [Plotly](https://plotly.com/): Plots interactive volcano curves

## License

This project is licensed under the MIT License - see the LICENSE.md file
for details.

## Funding

This material is based upon work supported by the Department of
Energy\'s Office of Energy Efficient and Renewable Energy\'s Advanced
Manufacturing Office under Award Number DE-EE0007888-9.5.
