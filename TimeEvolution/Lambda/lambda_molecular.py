import numpy as np
import lookup_table as lam_nu
import sys
from pathlib import Path

current_dir = str(Path.cwd())
root_dir = str(Path.cwd().parent.parent)
src_dir = root_dir + "/src"
lambda_dir = current_dir + "/Lambda"
data_dir = root_dir + "/MolecularData/"

sys.path.append(src_dir)
sys.path.append(lambda_dir)
sys.path.append(data_dir)

import physical_consts

def convert_nuc_mass_to_molecular_mass(nuc_masses, n_elec):
    
    """
    Convert dictionary of nuclear masses to a list of the full system masses (including electrons)

    Input:
        nuc_masses (dict) - dictionary of nuclear masses given in AMU
        n_elec (int) - number of electrons

    Output:
        mass_list (1d array) - array of particle masses in au

    """
    
    #handle nuclear masses
    molecular_masses = list(nuc_masses.values()) 
    molecular_masses = [physical_consts.convert_mass(nucleon) for nucleon in molecular_masses]    
    
    #handle electron masses
    for electron in range(n_elec):
        molecular_masses.insert(0,1) #append mass of 1 au to front of list
        
    return molecular_masses

def inverse_sum(array):
    
    """
    Compute sum of the inverses of a dictionary of values

    Input:
        array (array) - array of float values

    Output:
        inverse_sum (float) - inverse sum of the values

    """
    
    inv_sum = 0
    
    for i in range(len(array)):
        inv_sum += 1.0/array[i]

    return inv_sum

def lambda_T(num_elec, num_q, vol, masses):
    
    """
    Compute normalization parameter lambda for kinetic energy in a molecular Hamiltonian.
    
    Input:
        num_elec (int) - number of electrons in molecule
        num_q (int) - number of qubits used to represent a single dimension of a single particle
        vol (float) - volume of simulation cell in units
        masses (dict) - dictionary of masses of particles

    Output:
        lambda_t (float) - normalization parameter lambda for kinetic energy
        
    """

    pi_sq = np.square(np.pi)
    masses = convert_nuc_mass_to_molecular_mass(masses,num_elec)
    inverse_mass_sum = inverse_sum(masses)

    pre_fac =  (6 * pi_sq) / pow(vol,2.0/3.0)
    term_1 = np.square(np.power(2, num_q - 1) - 1)
    lambda_t = pre_fac * term_1 * inverse_mass_sum

    return lambda_t

def lambda_V(num_elec, lambda_nu, vol, charges):
    
    """
    Compute normalization parameter lambda for potential energy in a molecular Hamiltonian.
    
    Input:
        vol (float) - volume of simulation cell in units
        charges (array) - List of nuclear charges
        num_elec (int) - number of electrons in molecule
        lambda_nu (float) - lattice normalization parameter

    Output:
        lambda_v (float) - normalization parameter lambda for potential energy

    """

    numerator = (4*num_elec**2 - num_elec - sum(np.square(charges)))* lambda_nu
    denominator = 2 * np.pi * pow(vol, 1.0/3.0)
    lambda_v = numerator/denominator

    return lambda_v

def rotation_angle(num_states, rot_qubits):
    
    """
    Compute rotation angle used for preparing equal superposition states.
    Eq. (60) [Su, 2021]
    
    Input: 
        num_states (int) - number of non-zero coefficients
        rot_qubits (int) - number of qubits used for representing rotation angles
        
    Output:
        theta (float) - rotation angle
    """
    
    states_sq = np.power(2, np.ceil(np.log2(num_states)))
    dtheta = 2.0 * np.pi/(np.power(2, rot_qubits))
    arg = np.sqrt(states_sq / (4.0 * num_states))
    
    theta = dtheta * np.round(np.arcsin(arg)/dtheta)

    return theta

def equal_sup_sucess_prob(num_states, rot_qubits):
    
    """
    Compute success probability of preparing an equal superposition state.
    Eq. (59) [Su, 2021]
     
    Input:
        num_states (int) - number of non-zero coefficients
        rot_qubits (int) - number of qubits used for representing rotation angles
        
    Output:
        prob (float) - probability of successfully preparing an equal superposition of states
    """
    
    theta  = rotation_angle(num_states,rot_qubits)
    states_sq = np.power(2, np.ceil(np.log2(num_states)))
    
    prefac = num_states/states_sq
    term1 = 1.0 + (2 - 4.0 * num_states /states_sq) * np.square(np.sin(theta))
    term2 = np.square(np.sin(2.0 * theta ))
    
    prob = prefac * (np.square(term1) + term2)
    
    return prob

def total_sup_sucess_prob(rot_qubits, num_nuc, num_elec):
    """
    Compute total success probability of preparing an equal superposition state over all particles
    Eq. (63) [Su, 2021]
    
    Input:
        rot_qubits (int) - number of qubits used for representing rotation angles  
        num_nuc (int) - number of nuclei in the molecular system
        num_elec (int) - number of electrons in the molecular system
        
    Output:
        tot_prob (float) - total success probability for preparing an equal superposition state over all particles
    """
    
    num_particles = num_nuc + num_elec
    
    p1 = equal_sup_sucess_prob(3,8)
    p2 = equal_sup_sucess_prob(num_particles,rot_qubits)
    
    tot_prob = p1 * np.power(p2, 3)
    
    return tot_prob

def lambda_ham(vol, num_q, num_elec, num_nuc,
                molec_masses, lambda_nu, charges, b_R):
    """
    Compute normalization parameter lambda for a molecular Hamiltonian 

    Input: 
        vol (float) - volume of simulation cell in units
        num_q (int) - number of qubits used to represent a single dimension of a single particle
        num_elec (int) - number of electrons in molecule
        num_nuc (int) - number of nuclei in molecule
        molec_masses (array) - List of masses of particles in system
        lambda_nu (float) - lattice normalization parameter
        charges (array) - List of charges of particles in system
        b_R (int) - number of qubits used for representing rotation angles

    Output:
        lambda_molec (float) - lambda of molecule
    """

    #Compute normalization factors
    lambda_t = lambda_T(num_elec, num_q, vol, molec_masses)
    lambda_v = lambda_V(num_elec, lambda_nu, vol ,charges)
    success_probability = total_sup_sucess_prob(b_R, num_nuc, num_elec)

    lambda_molec = (lambda_t + lambda_v) / success_probability


    return lambda_molec

def lambda_alt(vol, b_R, num_nuc, num_elec, lambda_nu, charges, p_nu, a, p_Q):
    
    """
    Compute alternate normalization parameter lambda for a molecular Hamiltonian 

    Input: 
        vol (float) - volume of simulation cell in units
        b_R (int) - number of qubits used for representing rotation angles
        num_elec (int) - number of electrons in molecule
        num_nuc (int) - number of nuclei in molecule
        lambda_nu (float) - lattice normalization parameter
        charges (array) - List of charges of particles in system
        p_nu (float) - Success probability of creating nu state
        a (int) - 3 if amplitude amplification is performed on the nu state preparation, 1 otherwise.
        p_Q (float) - Success probability of creating chrage distribtution state

    Output:
        lambda_alt (float) - alternate lambda of molecule
    """
    
    #Check for amplitude amplification
    if a==3:
        p_nu = np.square(np.sin(3 * np.arcsin(np.sqrt(p_nu))))
    
    success_probability = total_sup_sucess_prob(b_R, num_nuc, num_elec)
    
    lambda_alt = lambda_V(num_elec, lambda_nu, vol, charges)/(p_nu * p_Q * success_probability)

    return lambda_alt
