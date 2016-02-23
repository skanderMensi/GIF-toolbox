# GIF-toolbox
toolbox to fit the Generalized-Integrate-and-Fire (GIF) model as described in: "Automated High-Throughput Characterization of Single Neurons by Means of Simplified Spiking Models"

For details, see:
"Automated High-Throughput Characterization of Single Neurons by Means of Simplified Spiking Models"
Christian Pozzorini, Skander Mensi, Olivier Hagens, Richard Naud, Christof Koch, Wulfram Gerstner
PLOS CB, 2015
DOI: 10.1371/journal.pcbi.1004275

6 files:

1. mainFullFit.m :        run the full fitting procedure
2. generateData.m :       generate surrogate data using a GIF model parameters obtain from L5 pyr neurons. This script generates a training set and a test set of input currents I and membrane potential V used to demonstrate the fitting procedure.
4. mainFitSubV.m :        parameter estimation of the subthreshold voltage dynamics
5. mainFitVThreshold.m :  parameter estimation of the threshold dynamics
6. mainTestMd.m :         compute performance Md* of the fitted model
7. mainVoltageTraces.m :  simulate the GIF model and generate figures

the /tools folder contains functions used by the main scripts.

Skander Mensi, 2015
