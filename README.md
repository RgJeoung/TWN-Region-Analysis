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

###### MD trajecoties
We used GROMACS version 5.1.4 and the CHARMM27 forcefield to get MD trajecotries. And you just need to prepare it with the trajectory in the DATA directory.

###### TWNs
TWNs can be obtained by running the TWN_Extractor_v1.exe file.
The input file must contain a_input in the upper PDB directory (for example, 1NVR in the trajectory directory) and must be entered in the format of "trajectory number.pdb" (for example, 0.pdb).

For example:

    TWN_Extractor_v1.exe -path ./DATA/trajectory/1NVR

The output file is created as "a_output" in the same pdb directory by default. In detail, you can set the output directory name and whether to include duplicates of the upper ring ("True" by default which means include duplicates).

For example:

    TWN_Extractor_v1.exe -path ./DATA/trajectory/1NVR -opname "output directory name" -dup False
    
After the extraction process is completed, you must change the rfour directory created in the output directory to the "protein-name_pdb-code_ring-type" format to create a TWN directory like "CHK1_1NVR_R4" in the DATA/TWN directory.
###### Boundary file
