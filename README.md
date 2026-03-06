FKM_ADI: Fast Krasnosel’skiĭ–Mann Accelerated ADI Method

This repository provides the MATLAB implementation and simulation codes associated with the paper: “An Fast Krasnosel’skiĭ–Mann-Based Policy Iteration Approach for Multi-Agent Control Systems.”

The repository contains the numerical codes used to evaluate the efficiency of the proposed FKM_ADI method and the simulation programs for the multi-agent consensus control experiments presented in the paper.

The goal of this repository is to ensure the reproducibility of the numerical results reported in the manuscript.
the repository contains two main folders:

1. FKM_ADI

This folder contains the MATLAB implementation of the proposed FKM_ADI algorithm and the numerical experiments used to evaluate its performance.

Main functionalities include:

Implementation of ADI, GADI, and FKM_ADI methods

Numerical experiments for solving Sylvester equations and linear matrix equations

Performance comparison in terms of

iteration numbers

CPU time

relative accuracy

These scripts reproduce the numerical results reported in Section 5 (Numerical Experiments) of the paper.
2. MASs_Simulation

This folder contains the simulation codes for the multi-agent control system (MASs) experiments.

The simulations demonstrate:

Optimal output consensus control

Policy iteration based control design

The effectiveness of the proposed accelerated framework in multi-agent systems

These scripts reproduce the simulation results presented in the control application section of the paper.












