import nbimporter
import numpy as np
from scipy import linalg
import sys
from pathlib import Path
from math import copysign, hypot
import math
import matplotlib.pyplot as plt
from itertools import product

# Configure project directories.
current_dir = str(Path.cwd())
root_dir = str(Path.cwd().parent.parent)
src_dir = root_dir + "/src"

# Add the shared source directory to the import path.
sys.path.append(src_dir)

import physical_consts
import parse_psi4
import io_functions as io
import isp_resources

def get_ssct_error(Sigma, delta, num_points, eta_n):
    """Compute the state-preparation error introduced by SSCT.

    Parameters
    ----------
    Sigma : array_like
        One-dimensional array containing the Gaussian coefficients.
    delta : float
        Momentum-space grid spacing.
    num_points : int
        Number of grid points.
    eta_n : int
        Number of nuclei, used to determine the number of degrees of freedom.

    Returns
    -------
    float
        SSCT error.
    """
    Sigma = np.asarray(Sigma, dtype=np.float64)

    dof = int(3 * eta_n)

    delta_sq = delta * delta

    m = num_points // 2

    l = np.arange(1, m, dtype=np.float64)

    l_minus = l * (l - 0.5)
    l_plus = l * (l + 0.5)

    m_minus = m * (m - 0.5)
    m_plus = m * (m + 0.5)

    log_ratio = 0.0
    

    for sigma in Sigma:

        c = delta_sq * sigma

        # Evaluate the numerator sum without constructing the full grid.
        s_num = (
            np.exp(-c * l_minus).sum()
            + np.exp(-c * l_plus).sum()
        )

        G_num = 1.0 + s_num

        # The normalization differs only through its two endpoint terms.
        extra = (
            np.exp(-c * m_minus)
            + np.exp(-c * m_plus)
        )

        # Work in log space to avoid cancellation when the ratio is near one.
        log_G_ratio = -np.log1p(extra / G_num)

        # Account for the exp(-c / 8) factor contributed by the numerator.
        log_ratio += log_G_ratio - c / 8.0

    value_squared = -2.0 * np.expm1(log_ratio)

    return np.sqrt(value_squared)

def get_single_point_parameters(Sigma, eta_n, n_qubits, L):
    """Evaluate the SSCT error for a given grid size and length.

    Parameters
    ----------
    Sigma : array_like
        Diagonal matrix of Gaussian coefficients.
    eta_n : int
        Number of nuclei.
    n_qubits : int
        Number of qubits used for each grid coordinate.
    L : float
        Grid length.

    Returns
    -------
    float
        SSCT error for the specified grid.
    """
    # Convert the diagonal coefficient matrix into a one-dimensional vector.
    sigma_vec = np.diag(Sigma)
    num_points = np.power(2, n_qubits)
    delta = 2 * np.pi / L
    current_ssct_error = get_ssct_error(sigma_vec, delta, num_points, eta_n)
    
    return current_ssct_error

def get_smallest_L_nuc_grid_parameters(
    Sigma,
    eta_n,
    n_low,
    n_high,
    L_low,
    L_high,
    ssct_budget,
    step_rate,
):
    """Find the smallest tested grid length satisfying the SSCT error budget.

    The search is repeated for every qubit count in the inclusive interval
    ``[n_low, n_high]``.

    Parameters
    ----------
    Sigma : array_like
        Diagonal matrix of Gaussian coefficients.
    eta_n : int
        Number of nuclei.
    n_low : int
        Smallest number of qubits to test.
    n_high : int
        Largest number of qubits to test.
    L_low : float
        Smallest grid length to test.
    L_high : float
        Largest grid length to test.
    ssct_budget : float
        Maximum allowed SSCT error.
    step_rate : float
        Increment between consecutive tested grid lengths.

    Returns
    -------
    numpy.ndarray
        Smallest successful grid length for each qubit count. An entry remains
        zero when no tested length satisfies the budget.
    numpy.ndarray
        Corresponding SSCT errors. An entry remains at its initialized value
        when no tested length satisfies the budget.
    """
    # Store one candidate length and error for each tested qubit count.
    best_lengths = np.zeros(int(n_high - n_low + 1))
    error_vec = np.zeros(int(n_high - n_low + 1))
    total_steps = int((L_high - L_low) / step_rate)

    # Convert the diagonal coefficient matrix into a one-dimensional vector.
    sigma_vec = np.diag(Sigma)

    # Search the length interval independently for each qubit count.
    for n in range(n_low, n_high + 1):
        indx = n - n_low
        num_points = np.power(2, n)
        best_L = 0
        best_error = 1e6
        L_not_found = True
        j_indx = 0
        
        # Stop at the first tested length that satisfies the error budget.
        while L_not_found:
            L = L_low + (step_rate * j_indx)
            delta = 2 * np.pi / L
            current_ssct_error = get_ssct_error(
                sigma_vec,
                delta,
                num_points,
                eta_n,
            )

            if current_ssct_error <= ssct_budget:
                best_L = L
                best_error = current_ssct_error
                L_not_found = False

            if j_indx >= total_steps:
                # Give up once the requested length interval is exhausted.
                print("Exceeded max steps for n = ", n)
                L_not_found = False

            j_indx += 1

        best_lengths[indx] = best_L
        error_vec[indx] = best_error

    return best_lengths, error_vec


def get_nuc_grid_parameters(
    Sigma,
    eta_n,
    n_low,
    n_high,
    L_low,
    L_high,
    ssct_budget,
    nuc_c_budget,
    K_max_elec,
    L_elec,
    step_rate,
):
    """Select nuclear grid parameters that satisfy the supplied constraints.

    Parameters
    ----------
    Sigma : array_like
        Diagonal matrix of nuclear Gaussian coefficients.
    eta_n : int
        Number of nuclei.
    n_low : int
        Smallest nuclear qubit count to test.
    n_high : int
        Largest nuclear qubit count to test.
    L_low : float
        Smallest nuclear grid length to test.
    L_high : float
        Largest nuclear grid length to test.
    ssct_budget : float
        Maximum allowed SSCT error.
    nuc_c_budget : float
        Maximum allowed nuclear cutoff error.
    K_max_elec : float
        Electronic momentum cutoff that the nuclear grid must support.
    L_elec : float
        Electronic grid length that the nuclear grid must cover.
    step_rate : float
        Increment between consecutive tested nuclear grid lengths.

    Returns
    -------
    numpy.integer
        Selected number of qubits per nuclear coordinate.
    float
        Selected nuclear grid length.
    float
        Selected nuclear momentum cutoff.
    """
    # First find the smallest length meeting the SSCT budget at each qubit count.
    L_nuc_vec, ssct_eps_vec = get_smallest_L_nuc_grid_parameters(
        Sigma,
        eta_n,
        n_low,
        n_high,
        L_low,
        L_high,
        ssct_budget,
        step_rate,
    )

    # Check that at least one candidate length satisfies the SSCT budget.
    no_L_found = np.all(L_nuc_vec == 0)

    if no_L_found:
        print(
            "ERROR: No L_nuc found that satisfies ssct_budget. "
            "Try increasing n_high and L_high."
        )

    # Retain candidates beginning with the first successful length.
    length_index_error_satisfied = np.flatnonzero(L_nuc_vec > 0)[0]
    qubits_error_satisfied = np.arange(n_low + length_index_error_satisfied, n_high + 1)
    L_nuc_satisfied = L_nuc_vec[length_index_error_satisfied:]

    # Compare the cutoff supported by each discrete grid with the analytical
    # cutoff required by the nuclear error bound.
    K_max_nuc_asymp_vec = np.zeros(len(L_nuc_satisfied))
    K_max_nuc_compute_vec = np.zeros(len(L_nuc_satisfied))

    indx = 0

    for L in L_nuc_satisfied:

        K_max_nuc_asymp_vec[indx] = isp_resources.get_nuc_k_max(
            nuc_c_budget,
            np.diag(Sigma),
            L,
            1,
        )
        K_max_nuc_compute_vec[indx] = (
            np.pi * np.power(2, qubits_error_satisfied[indx]) / L
        )

        indx += 1

    # Determine which candidates satisfy the nuclear and electronic cutoffs.
    K_asymp_satisfied = K_max_nuc_compute_vec >= K_max_nuc_asymp_vec
    K_elec_satisfied = K_max_nuc_compute_vec >= K_max_elec
    K_both_satisfied = K_asymp_satisfied == K_elec_satisfied

    # Check that at least one candidate satisfies both cutoff conditions.
    no_K_both_found = np.all(K_both_satisfied == False)

    if no_K_both_found:
        print("ERROR: No K_max_nuc was found that satisfies all conditions.")

    indx_both_satisfied = np.flatnonzero(K_both_satisfied)[0]

    # Set the final nuclear cutoff to the electronic cutoff, then recompute the
    # nuclear length implied by the selected qubit count.

    K_max = K_max_nuc_compute_vec[indx_both_satisfied] #sets K_max to largest K 
    n_nuc = qubits_error_satisfied[indx_both_satisfied]
    L_nuc = np.pi * np.power(2, n_nuc) / K_max

    # Verify that the recomputed length covers both prior length requirements.
    L_nuc_temp = L_nuc_satisfied[indx_both_satisfied]
    not_larger = (L_nuc < L_nuc_temp) or (L_nuc < L_elec)

    if not_larger:
        print("ERROR: L_nuc found is not larger than L_elec, or L needed to satisfy error bounds.")
        if L_nuc < L_nuc_temp:
            print("L_nuc < L_nuc_temp")
        else:
            print("L_nuc < L_elec")
  
    return n_nuc, L_nuc, K_max


def get_padded_grid_params(n_isp, R_mat):
    """Compute the additional grid qubits required after a linear transform.

    Parameters
    ----------
    n_isp : int
        Original number of qubits per grid coordinate.
    R_mat : array_like
        Transformation matrix used to determine the required padding.

    Returns
    -------
    numpy.float64
        Number of additional qubits required for the padded grid.
    """
    
    # Bound the transformed grid using the matrix infinity norm.
    num_pw = np.power(2, n_isp)
    n_pad = np.ceil(np.log2(num_pw * linalg.norm(R_mat, ord=np.inf) + 1)) - n_isp
    
    return n_pad
