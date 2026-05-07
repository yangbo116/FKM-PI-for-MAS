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

---

## 💻 Prerequisites

*   **MATLAB**: Developed and tested on **R2021a** or later.
   
## 🚀 How to Reproduce the Results

Follow these steps to replicate the numerical simulations:

### 1. Clone the Repository
Open your terminal or command prompt and run:
```bash
git clone https://github.com/yangbo116/FKM-PI-for-MAS.git
cd FKM-PI-for-MAS

### 2. Open MATLAB
Launch MATLAB.
Navigate to the FKM-PI-for-MAS folder.
Add the folder and its subfolders to the MATLAB path:
code
Matlab
addpath(genpath(pwd));

#### 3. Run the Simulation
To reproduce the numerical results, please follow these steps:
📂 Navigate to the test folder
In the MATLAB file browser, go to the directory:
FKM_ADI numerical test/

#### 🧪 Select the experiment **
Run Example1.m to replicate the results for Example 1.
Run Example2.m to replicate the results for Example 2.

#### 4. 💻 Execute the code
You can run the scripts by:
Typing the filename in the Command Window.
Or by opening the file and pressing F5.

### 📧 Contact
If you have any questions regarding the code, numerical methods, or the paper, please feel free to:
Open an Issue: Submit a bug report or a question directly on this GitHub repository.
Email the Author: Send an email to Yang Bo at bbo_yang@163.com.


