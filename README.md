# Estimating Photodissociation Resources

This repository estimates the logical-qubit and Toffoli-gate resources required to simulate photodissociation using the end-to-end quantum algorithm introduced in [1], *End-to-End Simulation of Chemical Dynamics on a Quantum Computer* ([arXiv:2603.19007](https://arxiv.org/abs/2603.19007)). Starting from a molecular structure, the workflow computes:

- equilibrium geometry and Hessian
- electronic MPS preparation resources
- nuclear initial-state preparation resources
- nuclear initial-state preparation error via Monte Carlo
- Hamiltonian simulation resources
- total logical-qubit and Toffoli-gate counts

The directory `MolecularData` contains a subdirectory for each molecule listed in Table 4 of [1]. Each subdirectory contains the data used to generate the corresponding resource estimates. These estimates can be reproduced by following the steps below using the data provided for each molecule. 

All equations and variables used throughout this respository have the meanings defined in [1].

## Citation

If you use this software, please cite:

[1] Elliot C. Eklund, Arkin Tikku, Patrick Sinnott, William J. Huggins, Guang Hao Low, Dominic W. Berry, Ivan Kassal, *End-to-End Simulation of Chemical Dynamics on a Quantum Computer*, (2026) ([arXiv:2603.19007](https://arxiv.org/abs/2603.19007)).

## 0. Preliminaries

### Software and Dependencies

The following software and packages are required:

- [Psi4](https://psicode.org/) (version 1.9.1) to compute equilibrium geometry and Hessian
- [PySCF](https://pyscf.org/) (version 2.9) to compute molecular orbitals for MPS encoding
- [Julia](https://julialang.org/downloads/) (version 1.10.4) for MPS encoding
    - [ITensors.jl](https://github.com/ITensor/ITensors.jl) (version 0.9.7)
    - [ITensorMPS.jl](https://github.com/ITensor/ITensorMPS.jl) (version 0.3.18)
    - [Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl) (version 4.1.0)
    - [SpecialPolynomials.jl](https://github.com/JuliaMath/SpecialPolynomials.jl) (version 0.5.0)
- [Mathematica](https://www.wolfram.com/mathematica/) (version 13.2) for Monte Carlo simulation
- [Jupyter Notebook](https://jupyter.org/) (version 6.4.12) to improve readability of Python code
- [Python](https://www.python.org/downloads/) (version 3.9.13) to compute resource estimates
    - [NumPy](https://numpy.org/) (version 1.21.5)
    - [SciPy](https://scipy.org/) (version 1.9.1)
    - [Matplotlib](https://matplotlib.org/) (version 3.5.2)

### Units

All quantities below, except in the Psi4 input files, are in atomic units.

### Example

The instructions below use the hydrogen molecule, H2, as an example.


## 1. Directory Setup and Input Files

The first step is to create the directories and files that are used throughout the calculation.

- Navigate to the `MolecularData` directory.

- Create a new directory named after the molecule. For the H2 example, `MolecularData/H2`.

- Inside `H2`, create the following input files:
  - `h2_psi4.in`
  - `h2_error_budget.in`
  - `h2_parameters.in`

- Create a new subdirectory named `Output`.

**Important:** All input file names must begin with `h2_`, matching the molecule directory name in lowercase. This naming convention is required and should be followed for every molecule.

### Populating `h2_psi4.in`

The `h2_psi4.in` file contains an initial guess for the ground-state equilibrium molecular geometry, along with the Psi4 commands required to optimize the structure and compute its vibrational frequencies and the Hessian. A minimal working example for H2 is provided below. The commands following the molecular geometry are required to output the Hessian. Further details on constructing Psi4 input files can be found in the [Psi4 documentation](https://psicode.org/psi4manual/master/).

```
molecule h2 {

H         -0.55410         0.00000         0.00000
H          0.55410         0.00000         0.00000
    
}

set basis cc-pvdz
optimize('mp2')
set hessian_write on
E, wfn = frequencies('mp2', return_wfn=True)
frequencies('mp2', return_wfn=True)
wfn.hessian().print_out()
```

### Populating `h2_error_budget.in`

The `h2_error_budget.in` file specifies the simulation error budget. The error-budget input file must contain the nine parameters shown below. `ErrorBudgetBreakdown.md` explains the complete breakdown of the overall error budget given in [1] and shows how the parameters in the example below fit within this breakdown.

```
total_budget=0.095
isp_budget=0.015
elec_isp_budget=1e-4
pkb_budget=1e-12
trim_budget=1e-4
ham_budget=0.00124
qsp_budget=1e-11
basis_trans_budget=10e-16
qae_qft_budget=1e-20
```

The nine parameters correspond to the following error terms in [1]:

| Parameter | Mathematical notation |
| :--- | :--- |
| `total_budget` | $\epsilon$ |
| `isp_budget` | $\epsilon_{\mathrm{ISP}}^{(\mathrm{ISP})}$ |
| `elec_isp_budget` | $\epsilon_{\mathrm{ISP}}^{(e)}$ |
| `pkb_budget` | $\epsilon_{\mathrm{PK}}$ |
| `trim_budget` | $\epsilon_{\mathrm{trim}}$ |
| `ham_budget` | $\epsilon_H$ |
| `qsp_budget` | $\epsilon_{\mathrm{QSP}}$ |
| `basis_trans_budget` | $\epsilon_B$ |
| `qae_qft_budget` | $\epsilon_{\mathrm{QAE-QFT}}$ |


### Populating `h2_parameters.in`

A parameters file must contain the 10 high-level simulation parameters listed below. An example for `h2_parameters.in` follows.

```
nuclei=2
basis=3-21g
spin=0
electrons=2
L_guess=100
E_guess=1400
time=30
is_linear=true
num_bonds=1
b_R=8
```

| Parameter | Description |
| :--- | :--- |
| `nuclei` | Number of nuclei in the molecule. |
| `basis` | Basis set used to calculate the molecular orbitals for electronic MPS construction. |
| `spin` | Molecular spin parameter, defined in PySCF as $N_\alpha-N_\beta=2S$. |
| `electrons` | Number of electrons in the molecule. |
| `L_guess` | Initial estimate of the real-space simulation-cell length required to satisfy the specified error bounds. Given in Bohr radii. |
| `E_guess` | Initial estimate of the maximum plane-wave kinetic-energy cutoff required to satisfy the specified error bounds. Given in Hartrees.|
| `time` | Total simulation time. Given in atomic units of time.|
| `is_linear` | Boolean parameter that is `true` when the molecule has a linear ground-state equilibrium geometry and `false` otherwise. |
| `num_bonds` | Number of bonds that are checked for dissociation. |
| `b_R` | Number of precision bits used for rotation angles when preparing uniform superposition states |


## 2. Compute Electronic Hessian

The second step is to run a Psi4 calculation that computes the Hessian of the electronic ground state. The output of this calculation will be used to compute the normal modes of the molecule.

- Navigate to `MolecularData/H2`.

- From the command line, run: `psi4 h2_psi4.in Output/h2_psi4.out`.

- If the calculation runs successfully, the directory `MolecularData/H2/Output` will contain with the following four files:
    - `h2_psi4.h2.####.hess`		
    - `h2_psi4.log`
    - `h2_psi4.h2.####.vibrec`
    - `h2_psi4.out`
    
- Here, `####` is a number that will be different for each calculation.
- The only file that needs to be kept is `h2_psi4.out`. The other three may be deleted.
- A file called `time.dat` will be created in `MolecularData/H2`. This file may also be deleted.



## 3. Electronic Matrix Product State Calculation

The next step is to compute the resources needed to prepare an MPS approximation to the electronic state. 

Portions of `elec_mps.jl` are derived from the [`fq-mps-stateprep`](https://github.com/wjhuggins/fq-mps-stateprep) repository developed by William Huggins and collaborators. The original software is distributed under the MIT License. The corresponding copyright and license notice is included in `THIRD_PARTY_LICENSES/fq-mps-stateprep-LICENSE`.

- Navigate to the directory `InitialStatePrep`.
- Run `julia elec_mps.jl H2 h2` in the terminal. Here, `H2` specifies the molecular-data directory, while `h2` specifies the filename prefix associated with the molecule.
- If the calculation completes successfully, it generates `h2_elec_mps.out` and `h2_elec_mps.log` in `MolecularData/H2/Output`. An example of each file is shown below.
- PySCF is called during `elec_mps.jl` to compute the MOs that are encoded as MPS.

Example `h2_elec_mps.out` file:

```
elec_isp_eps=6.771492802281801e-6
E_elec=200
L_elec=80.26769229921922
n_qubits_elec=9
elec_bond_dim_max=26
```

| Parameter | Mathematical notation | Description |
| :--- | :--- | :--- |
| `elec_isp_eps` | $\epsilon^{(e,c)} = \displaystyle\sum_{a=1}^{N_{\mathrm{MOB}}}\epsilon_a^{(e,c)}$ | Total error incurred by projecting the electronic wavefunction onto a truncated plane-wave basis (PWB) and encoding the resulting state as an MPS. This error is evaluated by explicitly performing the projection and MPS encoding. |
| `E_elec` | — | Maximum plane-wave kinetic energy used when calculating `elec_isp_eps`. Given in Hartrees. |
| `L_elec` | — | Real-space simulation-cell length used when calculating `elec_isp_eps`. Given in Bohr radii.|
| `n_qubits_elec` | — | Number of qubits per electron per spatial degree of freedom required to represent the electronic wavefunction using `E_elec` and `L_elec`. |
| `elec_bond_dim_max` | $\displaystyle\max_a M_a$ | Maximum MPS bond dimension across all molecular orbitals. |


Example `h2_elec_mps.log` file:

```
H2
h2
Running SCF
converged SCF energy = -0.826953624558593
Restarting SCF with second order method
converged SCF energy = -0.826953624558593
SCF complete
HOMO index: 1
Number of MOs found: 4
Error tol: 0.0001
Min Energy: 1400
Max Energy: 9800.0

Iteration Error: 1
Trying with energy: 1400
Thread 1 working on MO #1
0.99999999999652
Thread 1 working on MO #2
1.0000000000034777
Total State Error: 5.275473998687958e-6
Total Rotational Error: 3.590781520195377e-17
Total Error: 5.2754739987238655e-6
Total T gates: 6.5990518e7
Maximum Bond Dimension: 26

The computed error is smaller than the error budget. Try decreasing the energy to increase the error.


Finished
```

The file `elec_mps.jl` contains three hard-coded parameters found at the top of the file and chosen to make the resulting MPS compression error negligible relative to the target error budget, `total_budget`, used in [1].

```
ds_trunc = 1e-20
final_trunc = 1e-12
rot_qubits = 70
```

| Parameter | Mathematical notation | Description |
| :--- | :--- | :--- |
| `ds_trunc` | — | Truncation threshold used when compressing intermediate MPS representations. Singular values smaller than `ds_trunc` are discarded. |
| `final_trunc` | — | Truncation threshold used when compressing the final MPS representing a molecular orbital. The current description states that singular values smaller than `final_trunc` are discarded. |
| `rot_qubits` | $b_{\mathrm{rot}}^{(e)}$ | Number of rotation-synthesis precision qubits used during electronic MPS preparation. |


## 4. Compute Initial-State-Preparation Resources

The fourth step is to compute resource estimates for the initial-state preparation (ISP).

- Navigate to the directory `InitialStatePrep` 
- Open `compute_isp_resources.ipynb`. 
- Step-by-step instructions to execute the notebook and the definitions of the hardcoded parameters are contained within.

Executing `compute_isp_resources.ipynb` generates two output files: `h2_isp_final_data.out` and `h2_lambda_mat.out`. Both files are needed for upcoming steps. `h2_isp_final_data.out` contains final ISP resource estimates and `h2_lambda_mat.out` stores the diagonal matrix elements of the matrix that defines the nuclear Gaussian initial state. 


## 5. Perform Monte Carlo Error Estimation

The next step is to verify that the padded grid is large enough that the trimming error $\epsilon_{\mathrm{trim}}$ is below the corresponding error budget. This is done by performing a Monte Carlo simulation. 

- Navigate to the directory `InitialStatePrep`.
- Open `monte_carlo.nb`. 
- Step-by-step instructions to execute the notebook are given inside. 


## 6. Compute End-to-End Resource Estimation

The final step is to compute resource estimates for running the entire end-to-end algorithm. 

- Navigate to the directory `TimeEvolution`.
- Open `main.ipynb`. 
- Step-by-step instructions to execute the notebook are given inside. 

`main.ipynb` uses three hardcoded parameters. These parameters are the default inputs of the `get_error_budget()` function, in `resource_estimate.ipynb`. They are:

```
error_rot=1e-22,
error_gamma=1e-25,
error_phi=1e-25
```

`main.ipynb` generates `h2_quantum_yield_resources.out` in `MolecularData/H2/Output`. The beginning of the file includes a list of output parameters that should look like this:

```
Generated: 2026-10-01 17:05:34 AEST (UTC+1000)

------------
Total Resources
Total Qubits=3615
State Vector Qubits=192
Electron Qubits=96
Nuclear Qubits=96
Ancilla Qubits=3423
Toffoli Count=1.94e+11
Lambda=1387.1770868078356
delta_L=0.15707963267948966

------------
Electronic ISP Breakdown
Electronic:
Ancilla Qubits=277
Toffoli Count=8.56e+06

------------
Nuclear ISP Breakdown
Nuclear:
Ancilla Qubits=601
Toffoli Count=2.48e+09

------------
Time Evolution Breakdown
Qubits per dimension=16
Simulation cell length=1.03e+04
Ancilla Qubits=3412
Toffoli Count=8.94e+09

------------
Errors
Initial State Preparation Error=0.015
Hamiltonian Simulation Error=0.00124
QPE Error=0.0625
Total Quantum Yield Error=0.095
```

### Total Resources

| Variable | Description |
| :--- | :--- |
| `Total Qubits` | Total number of logical qubits required, including the state-vector and ancilla registers. |
| `State Vector Qubits` | Number of logical qubits used to encode the complete electron–nuclear wavefunction. |
| `Electron Qubits` | Total number of state-vector qubits used to represent the electrons. |
| `Nuclear Qubits` | Total number of state-vector qubits used to represent the nuclei. |
| `Ancilla Qubits` | Peak number of ancilla required during the calculation. |
| `Toffoli Count` | Total number of Toffoli gates required by the calculation. |
| `Lambda` | Hamiltonian normalization parameter used in the block encoding. Given in Hartrees. |
| `delta_L` | Real-space grid spacing associated with the simulation cell. Given in Bohr radii. |

### Electronic ISP

| Variable | Description |
| :--- | :--- |
| `Ancilla Qubits` | Peak number of ancilla required to prepare the electronic initial state. |
| `Toffoli Count` | Number of Toffoli gates required to prepare the electronic initial state. |

### Nuclear ISP

| Variable | Description |
| :--- | :--- |
| `Ancilla Qubits` | Peak number of ancilla qubits required to prepare the nuclear initial state. |
| `Toffoli Count` | Number of Toffoli gates required to prepare the nuclear initial state. |

### Time Evolution

| Variable | Description |
| :--- | :--- |
| `Qubits per Dimension` | Number of qubits used per particle per spatial dimension. |
| `Simulation Cell Length` | Length of the real-space simulation cell. |
| `Ancilla Qubits` | Peak number of ancilla qubits required to implement time evolution. |
| `Toffoli Count` | Number of Toffoli gates required to implement time evolution. |

### Errors

| Variable | Description |
| :--- | :--- |
| `Initial State Preparation Error` | Error allocated to preparing the initial electron–nuclear state. |
| `Hamiltonian Simulation Error` | Error arising from the approximation and simulation of the Hamiltonian. |
| `QPE Error` | Error allocated to quantum phase estimation. |
| `Total Quantum Yield Error` | Target upper bound on the total error in the estimated quantum yield. |


Following this output, all input files that were used to generate the output are copied into `h2_quantum_yield_resources.out` for easy record keeping.
