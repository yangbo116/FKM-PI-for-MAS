Numerical Simulation of FKM using ADI Method

This repository contains the MATLAB source code for the numerical simulation of An fast Krasnosel'skiĭ-Mann-based policy iteration approach for multi-agent
cooperative control systems

📁 Repository Structure
FKM_ADI numerical test/: Contains the main numerical experiments, convergence tests, and scripts to generate primary results.
MASs/: Contains utility functions, solvers, or modules related to Multi-Agent Systems (if applicable) used by the main scripts.
LICENSE: The legal terms for using this code.
README.md: This documentation file.

💻 Prerequisites
MATLAB: Developed and tested on version R2021a or later.

🚀 How to Reproduce the Results
To replicate the numerical results presented in this project, follow these steps:
1. Clone the Repository
code
Bash
git clone https://github.com/yangbo116/FKM-PI-for-MAS.git
cd FKM-PI-for-MAS

3. Set Up MATLAB Path
Open MATLAB and add the project folders to your path to ensure all functions in the MASs folder are recognized:
code
Matlab

% Run this in the MATLAB Command Window
addpath(genpath(pwd));
savepath;

5. Run the Simulation
Navigate to the FKM_ADI numerical test folder in the MATLAB file browser.
Locate the main entry script (The files Example1.m and Example2.m correspond to Example 1 and Example 2 in the manuscript, respectively).
Run the script by typing its name in the Command Window or pressing F5.

📊 Expected Output
Once the script finishes execution:
Figures: MATLAB will generate plots showing the evolution of the system, error analysis, or convergence orders.
Data: Numerical results may be saved in a .mat file or printed directly to the console.
(Optional: You can add a screenshot of a result plot here to show users what the output should look like)

📧 Contact
For any questions regarding the code or the numerical methods, please feel free to:
Open an Issue on this repository.
Contact the author at: bbo_yang@163.com














