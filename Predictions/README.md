# Code for predicting full depletion voltage evolution

This code aims at predicting the full depletion voltage evolution with radiations and time. The leakage current evolution is also computed here but only for predicting temperature evolution. A more fine temperature and fluence history would be needed for a precise prediction of the leakage current.

Precise input files names and options in `ComputeVdepEvolution.C` before compiling and executing.

## include

*include* directory contains the main codes.
*HamburgModel* is the class containing all the methods used for computing leakage current and charge carriers evolution.
*HamburgModelFactory* is the main class. There are method to read the various inputs, call the HamburgModel methods to do predictions for a given module, draw and save curves.

## Inputs

Directory to store localy the inputs: fluence/temperature scenario.
The scenario file format must be:
`StepIndex Fluence Temperature IsModulePowered`

