
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

def get_nuc_k_vec(eps_c, gamma, L, N_SMB):
    """Compute the nuclear momentum cutoff for each mode using Eq. (50).

    Parameters
    ----------
    eps_c : float
        Error budget for preprocessing the nuclear state.
    gamma : numpy.ndarray
        One-dimensional array of transformed Gaussian coefficients.
    L : float
        Nuclear grid length.
    N_SMB : int
        Number of single-modal basis functions.

    Returns
    -------
    numpy.ndarray
        Maximum nuclear momentum for each mode.
    """
    n = gamma.size
    k_vec = np.zeros(n)

    # The multi-basis correction vanishes for a single-modal basis.
    F = 0

    if N_SMB > 1:
        F = (N_SMB - 1) * np.log(4 * (N_SMB - 1))

    ind = 0

    for gam in gamma:
        log_term = np.log(
            1 / np.sqrt(2) + (2 * np.sqrt(np.pi) / (L * gam))
        )
        k_term = np.sqrt(2 * gam) * np.sqrt(
            2 * np.log(1 / eps_c) + (2 * log_term) + F + 5.4
        )

        k_vec[ind] = k_term
        ind += 1

    return k_vec

def get_nuc_k_max(eps_c, gamma, L, N_SMB):
    """Compute the largest nuclear momentum cutoff across all modes.

    Parameters
    ----------
    eps_c : float
        Error budget for preprocessing the nuclear state.
    gamma : numpy.ndarray
        One-dimensional array of transformed Gaussian coefficients.
    L : float
        Nuclear grid length.
    N_SMB : int
        Number of single-modal basis functions.

    Returns
    -------
    numpy.floating
        Largest nuclear momentum cutoff obtained from Eq. (50).
    """
    k_per_mode = get_nuc_k_vec(eps_c, gamma, L, N_SMB)
    nuc_k_max = np.max(k_per_mode)

    return nuc_k_max

def get_nuc_bond_dim_max(eps_c, gamma, L, N_SMB):
    """Compute the maximum nuclear MPS bond dimension using Eq. (55).

    Parameters
    ----------
    eps_c : float
        Error budget for preprocessing the nuclear state.
    gamma : numpy.ndarray
        One-dimensional array of transformed Gaussian coefficients.
    L : float
        Nuclear grid length.
    N_SMB : int
        Number of single-modal basis functions.

    Returns
    -------
    numpy.float64
        Maximum nuclear MPS bond dimension.
    """
    k_per_mode = get_nuc_k_vec(eps_c, gamma, L, N_SMB)
    bond_dim_per_mode = np.square(np.e) * np.square(k_per_mode) / gamma
    nuc_bond_dim_max = 2 * np.ceil(np.max(bond_dim_per_mode)) + 3

    return nuc_bond_dim_max


def get_shear_ortho_resources(num_atoms, n):
    """Estimate shear and orthogonal-transform resources from Eqs. (159)-(161).

    Parameters
    ----------
    num_atoms : int
        Number of nuclei.
    n : int
        Number of qubits per nucleus per spatial dimension.

    Returns
    -------
    int
        Number of required ancilla qubits.
    int
        Number of required Toffoli gates.
    """
    ancilla = int(np.ceil(4 * n - 4))

    T1 = (9 / 2) * np.square(num_atoms) * (8 * np.square(n) + (39 * n) - 8)
    T2 = (3 / 2) * num_atoms * (8 * np.square(n) + (35 * n) - 8)

    Toffolis = int(np.ceil(T1 - T2 - n))

    return ancilla, Toffolis



def get_prep_nuc_mps_resources(N_SMB, n, b_rot, nuc_bond_dim_max, num_atoms):
    """Estimate resources for nuclear MPS preparation using Eq. (136).

    Parameters
    ----------
    N_SMB : int
        Number of single-modal basis functions.
    n : int
        Number of qubits per nucleus per spatial dimension.
    b_rot : int
        Rotation-synthesis precision in qubits.
    nuc_bond_dim_max : int
        Maximum nuclear MPS bond dimension.
    num_atoms : int
        Number of nuclei.

    Returns
    -------
    int
        Number of required ancilla qubits.
    int
        Number of required Toffoli gates.
    """
    T1 = N_SMB * n
    T2 = (
        32
        * (1 + np.sqrt(2))
        * np.sqrt(b_rot + 1)
        * pow(nuc_bond_dim_max, 3 / 2)
    )
    T3 = (8 * b_rot - 15) * nuc_bond_dim_max * np.log(2 * nuc_bond_dim_max)

    Toffolis = int(np.ceil(3 * num_atoms * (T1 + 2 * n * (T2 + T3))))

    X1 = (1 / 2) * np.log(nuc_bond_dim_max)
    X2 = 2 * np.sqrt(nuc_bond_dim_max) * b_rot / np.sqrt(b_rot + 1)
    X3 = (1 / 2) * np.log(b_rot + 1)
    X4 = 3 * b_rot

    ancilla = int(np.ceil(X1 + X2 + X3 + X4))

    return ancilla, Toffolis



def get_b_rot_nuc(eps_q, nuc_bond_dim_max, n, num_atoms):
    """Compute the nuclear MPS rotation-synthesis precision.

    Parameters
    ----------
    eps_q : float
        Target unitary-synthesis error.
    nuc_bond_dim_max : int
        Maximum nuclear MPS bond dimension.
    n : int
        Number of qubits per nucleus per spatial dimension.
    num_atoms : int
        Number of nuclei.

    Returns
    -------
    int
        Number of rotation-synthesis qubits.
    """
    T1 = 3 * num_atoms * (
        pow(2, 7 / 2)
        * n
        * nuc_bond_dim_max
        * np.log2(2 * nuc_bond_dim_max)
    )
    b_rot = int(np.ceil(np.log2(T1 / eps_q)))

    return b_rot


def get_pk_resources(b_grad, n, eps_pk, num_atoms):
    """Estimate phase-kickback resources using Eqs. (163) and (164).

    Parameters
    ----------
    b_grad : int
        Number of gradient-precision qubits used for phase kickback.
    n : int
        Number of qubits per nucleus per spatial dimension.
    eps_pk : float
        Target phase-kickback error.
    num_atoms : int
        Number of nuclei.

    Returns
    -------
    int
        Number of required ancilla qubits.
    int
        Number of required Toffoli gates.
    """
    ancilla = int(np.ceil(4 * b_grad + n - 1))

    T1 = (3 * num_atoms) * (2 * n * b_grad + b_grad - n)
    T2 = 1.149 * (
        np.log2(3 * num_atoms)
        + np.log2(b_grad)
        + np.log2(1 / eps_pk)
    ) + 9.2

    Toffolis = int(np.ceil(T1 + (b_grad * T2 / 2)))

    return ancilla, Toffolis


def get_nuc_resources(
    N_SMB,
    n_p_nuc,
    n_p_nct,
    b_rot_nuc,
    b_grad,
    nuc_bond_dim_max,
    eta_n,
    eps_pk,
):
    """Combine all nuclear initial-state-preparation resource estimates.

    Parameters
    ----------
    N_SMB : int
        Number of single-modal basis functions.
    n_p_nuc : int
        Number of qubits per nuclear coordinate before transformation.
    n_p_nct : int
        Number of qubits per nuclear coordinate during the coordinate
        transformation.
    b_rot_nuc : int
        Rotation-synthesis precision for nuclear MPS preparation.
    b_grad : int
        Gradient-precision qubits used for phase kickback.
    nuc_bond_dim_max : int
        Maximum nuclear MPS bond dimension.
    eta_n : int
        Number of nuclei.
    eps_pk : float
        Target phase-kickback error.

    Returns
    -------
    int
        Nuclear MPS ancilla cost.
    int
        Nuclear MPS Toffoli cost.
    int
        Ancilla cost.
    int
        Toffoli cost.
    int
        Phase-kickback ancilla cost.
    int
        Phase-kickback Toffoli cost.
    int or float
        Total nuclear ISP Toffoli cost.
    numpy.integer
        Peak nuclear ISP ancilla cost.
    int
        Number of nuclear data qubits.
    """
    # Estimate the resources for preparing the nuclear MPS.
    W_nuc_anc, W_nuc_Toff = get_prep_nuc_mps_resources(
        N_SMB,
        n_p_nuc,
        b_rot_nuc,
        nuc_bond_dim_max,
        eta_n,
    )

    # Estimate the resources for applying A^{-T}.
    U_A_inv_trans_anc, U_A_inv_trans_Toff = get_shear_ortho_resources(eta_n, n_p_nct)

    # Estimate the phase-kickback resources.
    PK_anc, PK_Toff = get_pk_resources(b_grad, n_p_nct, eps_pk, eta_n)

    # Add sequential Toffoli costs and take the peak reusable ancilla count.
    nuc_isp_Toff = W_nuc_Toff + U_A_inv_trans_Toff + PK_Toff
    nuc_isp_anc = np.max((W_nuc_anc, U_A_inv_trans_anc, PK_anc))
    nuc_data_qubits = 3 * eta_n * n_p_nuc

    return (
        W_nuc_anc,
        W_nuc_Toff,
        U_A_inv_trans_anc,
        U_A_inv_trans_Toff,
        PK_anc,
        PK_Toff,
        nuc_isp_Toff,
        nuc_isp_anc,
        nuc_data_qubits,
    )


def get_asym_resources(n_p, eta_e):
    """Estimate antisymmetrization resources using Eqs. (80) and (91).

    Parameters
    ----------
    n_p : int
        Number of qubits per electron per spatial dimension.
    eta_e : int
        Number of electrons.

    Returns
    -------
    int
        Number of required ancilla qubits.
    numpy.floating
        Number of required Toffoli gates.
    """
    n_p_bar = pow(2, np.ceil(np.log2(n_p + 1)))

    T1 = eta_e * np.log2(eta_e)
    T2 = n_p_bar * np.log2(n_p_bar) * (1 + np.log2(n_p_bar)) / 4
    T3 = 2 * (eta_e - 1)

    ancillas = int(np.ceil(T1 + T2 + T3))

    X1 = 2 * (eta_e - 1) * (np.log2(n_p_bar) + 1)
    X2 = (
        n_p_bar
        * np.log2(n_p_bar)
        * (1 + np.log2(n_p_bar))
        * (6 * np.log2(n_p_bar) + n_p + 1)
        / 4
    )

    Toffolis = X1 + X2

    return ancillas, Toffolis


def prep_elec_mps_resources(N_MOB, n, b_rot, elec_bond_dim_max, num_elec):
    """Estimate resources for electronic MPS preparation using Eq. (90).

    Parameters
    ----------
    N_MOB : int
        Number of molecular-orbital basis functions.
    n : int
        Number of qubits per electron per spatial dimension.
    b_rot : int
        Rotation-synthesis precision in qubits.
    elec_bond_dim_max : int
        Maximum electronic MPS bond dimension.
    num_elec : int
        Number of electrons.

    Returns
    -------
    int
        Number of required ancilla qubits.
    int
        Number of required Toffoli gates.
    """
    T1 = num_elec * N_MOB * n
    T2 = (
        32
        * (1 + np.sqrt(2))
        * np.sqrt(b_rot + 1)
        * pow(elec_bond_dim_max, 3 / 2)
    )
    T3 = (8 * b_rot - 15) * elec_bond_dim_max * np.log(2 * elec_bond_dim_max)

    Toffolis = int(np.ceil(T1 + 2 * num_elec * n * (T2 + T3)))

    X1 = (1 / 2) * np.log(elec_bond_dim_max)
    X2 = 2 * b_rot * np.sqrt(elec_bond_dim_max) / np.sqrt(b_rot + 1)
    X3 = (1 / 2) * np.log(b_rot + 1)
    X4 = 3 * b_rot

    ancilla = int(np.ceil(X1 + X2 + X3 + X4))

    return ancilla, Toffolis


def get_elec_resources(n_p, eta_e, N_MOB, b_rot_elec, elec_bond_dim_max):
    """Combine all electronic initial-state-preparation resource estimates.

    Parameters
    ----------
    n_p : int
        Number of qubits per electron per spatial dimension.
    eta_e : int
        Number of electrons.
    N_MOB : int
        Number of molecular-orbital basis functions.
    b_rot_elec : int
        Rotation-synthesis precision for electronic MPS preparation.
    elec_bond_dim_max : int
        Maximum electronic MPS bond dimension.

    Returns
    -------
    int
        Antisymmetrization ancilla cost.
    numpy.floating
        Antisymmetrization Toffoli cost.
    int
        Electronic MPS ancilla cost.
    int
        Electronic MPS Toffoli cost.
    int or float
        Total electronic ISP Toffoli cost.
    numpy.integer
        Peak electronic ISP ancilla cost.
    int
        Number of electronic data and spin qubits.
    """
    ASYM_anc, ASYM_toff = get_asym_resources(n_p, eta_e)
    W_elec_anc, W_elec_toff = prep_elec_mps_resources(
        N_MOB,
        n_p,
        b_rot_elec,
        elec_bond_dim_max,
        eta_e,
    )

    # Add sequential Toffoli costs and take the peak reusable ancilla count.
    elec_isp_Toff = W_elec_toff + ASYM_toff
    elec_isp_anc = np.max((W_elec_anc, ASYM_anc))
    # Include one spin qubit for each electron in addition to spatial data.
    elec_data_qubits = 3 * eta_e * n_p + eta_e

    return (
        ASYM_anc,
        ASYM_toff,
        W_elec_anc,
        W_elec_toff,
        elec_isp_Toff,
        elec_isp_anc,
        elec_data_qubits,
    )

def get_total_resources(
    elec_isp_Toff,
    elec_isp_anc,
    elec_data_qubits,
    nuc_isp_Toff,
    nuc_isp_anc,
    nuc_data_qubit,
):
    """Combine the nuclear and electronic ISP resource estimates.

    Parameters
    ----------
    elec_isp_Toff : int or float
        Electronic initial-state-preparation Toffoli cost.
    elec_isp_anc : int
        Electronic initial-state-preparation ancilla cost.
    elec_data_qubits : int
        Number of electronic data qubits.
    nuc_isp_Toff : int or float
        Nuclear initial-state-preparation Toffoli cost.
    nuc_isp_anc : int
        Nuclear initial-state-preparation ancilla cost.
    nuc_data_qubit : int
        Number of nuclear data qubits.

    Returns
    -------
    int or float
        Total initial-state-preparation Toffoli cost.
    numpy.integer
        Peak initial-state-preparation ancilla cost.
    int
        Total number of data qubits.
    """
    isp_Toff = elec_isp_Toff + nuc_isp_Toff
    isp_anc = np.max((elec_isp_anc, nuc_isp_anc))
    isp_data = elec_data_qubits + nuc_data_qubit

    return isp_Toff, isp_anc, isp_data
