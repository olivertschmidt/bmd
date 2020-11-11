# Bispectral Mode Decomposition
`BMD()` Bispectral Mode Decomposition of multidimensional and multivariate snapshot data

`CBMD()` Cross-Bispectral Mode Decomposition of two or three fields

Triadic interactions are the fundamental mechanism of energy transfer in fluid flows. Bispectral mode decomposition (BMD) educes coherent flow structures that are associated with triadic interactions from experimental or numerical data. Triadic interactions are characterized by quadratic phase coupling which can be detected by the bispectrum. The proposed method maximizes an integral measure of this third-order statistic to compute modes associated with frequency triads, as well as a mode bispectrum that identifies resonant three-wave interactions. Unlike the classical bispectrum, the decomposition establishes a causal relationship between the three frequency components of a triad. This permits the distinction of sum- and difference-interactions, and the computation of interaction maps that indicate regions of nonlinear coupling. 

## Files
| File        |     Description     |
| ------------- |:-------------|
| bmd.m | Bispectral mode decomposition | 
| cbmd.m | Cross-bispectral mode decomposition | 
| example_1.m | BMD of single-variable, 2-D test case with default spectral estimation parameters | 
| example_2.m | BMD of 2-variable, 2-D test case with manual specification of spectral estimation parameters and interactive visualization of modes | 
| example_3.m | Demonstration of CBMD | 
| wake_Re500.mat | Snapshot data of the flow behind a cylinder at Re=500 | 
| LICENSE.txt | License | 

## Cite As
[1] Schmidt, O. T., Bispectral mode decomposition of nonlinear flows, Nonlinear Dynamics, 2020  
