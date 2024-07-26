# TWN-Region-Analysis
Source code for TWN Region Analysis. 
Anyone can test the TWN Region Analysis for fragment growing.

TWN region analysis is conducted through the following procedure.
## Installation
This code was tested in windows with Python 3.9.

A yaml file containing all requirements is provided.

You can set up a virtual environment in conda as follows:

    conda env create -f TWN-Region-Analysis.yaml
    conda activate TWN-Region-Analysis

## Preset
To run this code, you need the following three preset files.
1. Molecular Dynamics (MD) simulation trajectories for target protein
2. Topological Water Networks (TWN)
3. Boundary File

We used GROMACS version 5.1.4 and the CHARMM27 forcefield to get MD trajecotries. And you just need to prepare it with the trajectory in the DATA directory.

