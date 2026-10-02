
import nbimporter
import numpy as np
from scipy import linalg
import sys
from pathlib import Path
from math import copysign, hypot
import math

# Configure project directories.
current_dir = str(Path.cwd())
root_dir = str(Path.cwd().parent.parent)
src_dir = root_dir + "/src"

# Add the shared source directory to the import path.
sys.path.append(src_dir)

import physical_consts
import parse_psi4
import io_functions as io
from scipy.linalg import ldl


def get_mass_weight_mat(molec_masses, eta_n):
    """Construct the mass-weighting matrix and its inverse in atomic units.

    Parameters
    ----------
    molec_masses : dict
        Mapping whose values are the nuclear masses in atomic mass units.
    eta_n : int
        Number of nuclei.

    Returns
    -------
    numpy.ndarray
        Diagonal mass-weighting matrix with entries ``1 / sqrt(m_i)`` repeated
        for the three Cartesian coordinates of each nucleus.
    numpy.ndarray
        Inverse mass-weighting matrix with entries ``sqrt(m_i)``.
    """
    mass_vec = np.zeros(eta_n)

    for i, val in enumerate(molec_masses.values()):
        mass_vec[i] = physical_consts.convert_mass(val)

    dof = 3 * eta_n
    mass_weight_mat = np.zeros((dof, dof))
    mass_weight_mat_inv = np.zeros((dof, dof))

    # Repeat each nuclear mass along its three Cartesian coordinates.
    for i in range(eta_n):
        mass_weight_mat[i * 3, i * 3] = 1.0 / np.sqrt(mass_vec[i])
        mass_weight_mat[i * 3 + 1, i * 3 + 1] = 1.0 / np.sqrt(mass_vec[i])
        mass_weight_mat[i * 3 + 2, i * 3 + 2] = 1.0 / np.sqrt(mass_vec[i])

    # Construct the inverse explicitly from the same converted masses.
    for i in range(eta_n):
        mass_weight_mat_inv[i * 3, i * 3] = np.sqrt(mass_vec[i])
        mass_weight_mat_inv[i * 3 + 1, i * 3 + 1] = np.sqrt(mass_vec[i])
        mass_weight_mat_inv[i * 3 + 2, i * 3 + 2] = np.sqrt(mass_vec[i])

    return mass_weight_mat, mass_weight_mat_inv

def get_mass_weight_hess(pre_weighted_hessian, mass_weight_mat):
    """Apply mass weighting to a molecular Hessian.

    Parameters
    ----------
    pre_weighted_hessian : numpy.ndarray
        Molecular Hessian before mass weighting.
    mass_weight_mat : numpy.ndarray
        Diagonal mass-weighting matrix.

    Returns
    -------
    numpy.ndarray
        Mass-weighted Hessian.
    """
    He = np.matmul(mass_weight_mat, np.matmul(pre_weighted_hessian, mass_weight_mat))

    return He

def check_freqs(molec_masses, eta_n, pre_weighted_hessian, is_linear=False):
    """Compute and print vibrational frequencies from a molecular Hessian.

    Parameters
    ----------
    molec_masses : dict
        Mapping whose values are the nuclear masses in atomic mass units.
    eta_n : int
        Number of nuclei.
    pre_weighted_hessian : numpy.ndarray
        Molecular Hessian before mass weighting.
    is_linear : bool, optional
        Whether molecule has linear or non-linear geometry

    Returns
    -------
    numpy.ndarray
        All sorted frequencies, including translational and rotational modes,
        in inverse centimetres.
    numpy.ndarray
        Vibrational frequencies after removing the translational and rotational
        modes, in inverse centimetres.
    """
    mass_weight_mat, _ = get_mass_weight_mat(molec_masses, eta_n)
    He = get_mass_weight_hess(pre_weighted_hessian, mass_weight_mat)
    omega_sq_mw, O = linalg.eig(He)

    freqs = physical_consts.hartree2wavenum(np.sqrt(np.abs(omega_sq_mw)))
    freqs_contam = np.sort(freqs)
    freqs = freqs_contam[6:]

    # Linear molecules have five translational and rotational modes.
    if is_linear:
        freqs = freqs_contam[5:]

    print("Frequencies in cm^-1:", freqs, "\n")

    return freqs_contam, freqs

def clean_mat(A, thresh):
    """Set matrix entries below a magnitude threshold to zero in place.

    Parameters
    ----------
    A : numpy.ndarray
        Matrix to clean. This array is modified in place.
    thresh : float
        Magnitude below which entries are set to zero.

    Returns
    -------
    numpy.ndarray
        The modified input matrix.
    """
    A[np.abs(A) < thresh] = 0

    return A

def get_RQ_decomp(A):
    """Compute the RQ decomposition of a square matrix.

    Parameters
    ----------
    A : numpy.ndarray
        Square matrix to decompose.

    Returns
    -------
    numpy.ndarray
        Upper-triangular matrix ``R`` satisfying ``A = R @ Q``.
    numpy.ndarray
        Orthogonal matrix ``Q``.
    """
    R, Q = linalg.rq(A)

    return R, Q

def get_DR_decomp(R):
    """Factor an upper-triangular matrix into diagonal and triangular parts.

    Parameters
    ----------
    R : numpy.ndarray
        Upper-triangular matrix to decompose.

    Returns
    -------
    numpy.ndarray
        Scaled upper-triangular factor ``r`` satisfying ``R = D @ r``.
    numpy.ndarray
        Diagonal factor ``D`` containing the absolute diagonal entries of
        ``R``.
    """
    D = np.diag(np.abs(np.diag(R)))
    D_inv = np.diag(1.0 / np.abs(np.diag(R)))
    r = np.matmul(D_inv, R)

    return r, D

def get_DrQ_decomp(A):
    """Compute the ``D @ r @ Q`` decomposition of a square matrix.

    Parameters
    ----------
    A : numpy.ndarray
        Square matrix to decompose.

    Returns
    -------
    numpy.ndarray
        Orthogonal factor ``Q``.
    numpy.ndarray
        Scaled upper-triangular factor ``r``.
    numpy.ndarray
        Diagonal factor ``D``.
    """
    R, Q = get_RQ_decomp(A)
    r, D = get_DR_decomp(R)

    return Q, r, D

def get_nuclear_wfxn_parameters(
    molec_masses,
    eta_n,
    pre_weighted_hessian,
    round_tol,
    is_linear=False,
):
    """Compute transformation data and Gaussian widths for the nuclear state.

    Parameters
    ----------
    molec_masses : dict
        Mapping whose values are the nuclear masses in atomic mass units.
    eta_n : int
        Number of nuclei.
    pre_weighted_hessian : numpy.ndarray
        Molecular Hessian before mass weighting.
    round_tol : float
        Magnitude below which transformation-matrix entries are set to zero.
    is_linear : bool, optional
        Whether molecule has linear or non-linear geometry

    Returns
    -------
    numpy.ndarray
        Cleaned orthogonal factor from the ``D @ r @ Q`` decomposition.
    numpy.ndarray
        Cleaned upper-triangular factor from the decomposition.
    numpy.ndarray
        Gaussian coefficients for all nuclear degrees of freedom.
    """
    MW, IMW = get_mass_weight_mat(molec_masses, eta_n)
    He = get_mass_weight_hess(pre_weighted_hessian, MW)
    omega_sq_mw, O = linalg.eig(He)

    # Transform the normal-mode eigenvectors back to unweighted coordinates.
    W = np.matmul(O, IMW)
    W = clean_mat(W, round_tol)

    O_bar, SU, d = get_DrQ_decomp(W)

    O_bar_clean = clean_mat(O_bar, round_tol)
    SU_clean = clean_mat(SU, round_tol)
    d_clean = clean_mat(d, round_tol)

    omega_mw = np.sqrt(np.real(np.abs(omega_sq_mw)))
    # The squared diagonal factor follows from the expression above Eq. (127).
    gamma_dirty = np.multiply(np.square(np.diag(d_clean)), omega_mw)

    # Sort the modes before separating vibrations from translations/rotations.
    sorted_indices = np.argsort(omega_sq_mw)

    omega_sq_clean = []
    F_clean = []
    gamma = []

    # A linear molecule has three translations and two rotations.
    if is_linear:
        gamma = gamma_dirty[sorted_indices][5:]

        # Use representative vibrational widths for the removed modes.
        gamm_trans = np.min(gamma)
        gamma_rot = np.max(gamma)

        # Append three translational widths.
        for i in range(3):
            gamma = np.append(gamma, gamm_trans)

        # Append two rotational widths.
        for i in range(2):
            gamma = np.append(gamma, gamma_rot)

    # A nonlinear molecule has three translations and three rotations.
    else:
        gamma = gamma_dirty[sorted_indices][6:]

        # Use representative vibrational widths for the removed modes.
        gamm_trans = np.min(gamma)
        gamma_rot = np.max(gamma)

        # Append three translational widths.
        for i in range(3):
            gamma = np.append(gamma, gamm_trans)

        # Append three rotational widths.
        for i in range(3):
            gamma = np.append(gamma, gamma_rot)

    return O_bar_clean, SU_clean, gamma

def get_ssct_parameters(
    molec_masses,
    eta_n,
    pre_weighted_hessian,
    round_tol,
    is_linear,
):
    """Compute the Gaussian and shear matrices required by the SSCT.

    Parameters
    ----------
    molec_masses : dict
        Mapping whose values are the nuclear masses in atomic mass units.
    eta_n : int
        Number of nuclei.
    pre_weighted_hessian : numpy.ndarray
        Molecular Hessian before mass weighting.
    round_tol : float
        Magnitude below which transformation-matrix entries are set to zero.
    is_linear : bool
        Whether molecule has linear or non-linear geometry

    Returns
    -------
    numpy.ndarray
        Block-diagonal factor from the LDL decomposition of ``Lambda``.
    numpy.ndarray
        Upper-triangular shear matrix obtained by transposing the LDL lower
        factor.
    """
    O_bar, SU, gamma_ps = get_nuclear_wfxn_parameters(
        molec_masses,
        eta_n,
        pre_weighted_hessian,
        round_tol,
        is_linear,
    )

    # Convert position-space Gaussian coefficients to momentum-space values.
    gamma_ms = linalg.inv(np.diag(gamma_ps))

    # Form the composite coordinate transformation A.
    A = np.matmul(SU, O_bar)

    # Transform the momentum-space Gaussian matrix into Cartesian coordinates.
    A_inv = linalg.inv(A)
    A_inv_trans = np.transpose(A_inv)
    Lambda = np.matmul(A_inv, np.matmul(gamma_ms, A_inv_trans))

    # The KW procedure is implemented here as an LDL decomposition.
    L_mat, Sigma, _ = ldl(Lambda, lower=True)
    R_mat = np.transpose(L_mat)

    return Sigma, R_mat
