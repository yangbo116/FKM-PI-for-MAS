## 📁 Repository Structure

- **`FKM_ADI numerical test/`**: Contains the main numerical experiments and the scripts for reproducing the results in Examples 1 and 2 of the manuscript.

- **`MASs/`**: Contains the main numerical experiments and the scripts for reproducing the results in Examples 3 and 4 of the manuscript.
- **`LICENSE`**: MIT License.
- **`README.md`**: Project documentation.

---

## 💻 Prerequisites

- **MATLAB**: Developed and tested on **R2021a** or later.

---

## 🚀 How to Reproduce the Results

Follow these steps to replicate the numerical experiments:

### 1. Clone the Repository

Open your terminal or command prompt and run:

```bash
git clone https://github.com/yangbo116/FKM-PI-for-MAS.git
cd FKM-PI-for-MAS
```

### 2. Open MATLAB

- Launch MATLAB.
- Navigate to the `FKM-PI-for-MAS` or `MASs` folder.
- Add the folder and its subfolders to the MATLAB path:

```matlab
addpath(genpath(pwd));
```

### 3. Run the numerical experiments

To reproduce the numerical results, please follow these steps:

#### 📂 Step 1 — Navigate to the test folder

In the MATLAB file browser, go to the directory: 
- FKM_ADI numerical test/
- MASs/
#### 🧪 Step 2 — Select the experiment

- Run `Example1.m` to replicate the results for **Example 1**.
- Run `Example2.m` to replicate the results for **Example 2**.
- Run `example_mas.m` to replicate the results for **Example 3**.
- Run `example_mas_mechanical.m` to replicate the results for **Example 4**.

#### ▶️ Step 3 — Execute the code

You can run the scripts by:

- Entering the script name in the MATLAB Command Window and pressing Enter, or
- Opening the script file and pressing F5 (or clicking Run).

---

## 📧 Contact

If you have any questions regarding the code, numerical methods, or the paper, please feel free to:

- **Open an Issue**: Submit a bug report or question directly on this GitHub repository.
- **Email the Author**: Send an email to [bbo_yang@163.com](mailto:bbo_yang@163.com).





