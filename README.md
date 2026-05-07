# A Fast Krasnosel'skiĭ-Mann-based Policy Iteration (FKM-PI) for Multi-Agent Systems

[![MATLAB](https://img.shields.io/badge/MATLAB-R2021a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

This repository contains the official MATLAB implementation for the numerical simulations of the paper:  
**"A fast Krasnosel'skiĭ-Mann-based policy iteration approach for multi-agent cooperative control systems."**

The project provides an efficient numerical framework to solve cooperative control problems in Multi-Agent Systems (MASs) using an accelerated policy iteration (PI) algorithm based on the Krasnosel'skiĭ-Mann (KM) iteration.

## 📁 Repository Structure

*   **`FKM_ADI numerical test/`**: Contains the main numerical experiments and scripts to generate the results presented in the manuscript.
    *   `Example1.m`: Script to reproduce results for Example 1.
    *   `Example2.m`: Script to reproduce results for Example 2.
*   **`MASs/`**: A library of utility functions, solvers, and core modules specifically designed for Multi-Agent System dynamics and control.
*   **`LICENSE`**: MIT License.
*   **`README.md`**: Project documentation.

## 💻 Prerequisites

*   **MATLAB**: Developed and tested on **R2021a** or later.
*   **Toolboxes**: (Recommended) Control System Toolbox.

## 🚀 How to Reproduce the Results

Follow these steps to replicate the numerical simulations:

### 1. Clone the Repository
```bash
git clone https://github.com/yangbo116/FKM-PI-for-MAS.git
cd FKM-PI-for-MAS
