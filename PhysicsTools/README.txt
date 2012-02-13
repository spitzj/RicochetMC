README for Neutrino MC Libraries
Adam Anderson
adama@mit.edu

These is the beginning of a simple library for simulating reactor neutrino experiments with short baseline neutrino oscillations.  ReactorExperiment contains methods for generating fake Monte Carlo data, and OscillationFitter contains methods for fitting the fake experimental data.  ReactorNuAnalysis.C contains the main method, which does a random sample analysis using the two classes.  Everything should compile easily by typing 'make', assuming that you have ROOT and the usual UNIX compilers installed.