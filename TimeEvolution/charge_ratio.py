
import numpy as np
import math
import nbimporter
import sys

import physical_consts
import parse_psi4
import io_functions as io

def compute_charge_ratio(sim_parameters,molec_labels):
    """
    Compute charge ratio for molecule specified by parameter_file and psi4_out_file.
    
    Input:
        sim_parameters - Parameters for simulation of molecule
        molec_labels - List of labels for molecules
        
    Output:
        charge_ratio (float) - charge ratio of given molecule
    """

    #Unpack data
    num_elec = sim_parameters["electrons"]
    num_nuc = sim_parameters["nuclei"]
    
    #Compute denomenator
    denom = 0
    
    #nuc-nuc sum: all molecules are neutral, sum of nuclear charges = num_elec
    denom += np.square(num_elec)
            
    #nuc-elec sum
    denom += 2*np.square(num_elec)
    
    #elec-elec sum
    denom += np.square(num_elec)
    
        
    #Compute numerator
    num = 0
    
    #nuc-nuc sum
    for i in range(num_nuc):
        
        nucleon_i = molec_labels[i]
        charge_i = physical_consts.nuclear_charge(nucleon_i)
        
        for j in range(num_nuc):
        
            nucleon_j = molec_labels[j]
            charge_j = physical_consts.nuclear_charge(nucleon_j)
            
            if i != j:
                num += np.abs(charge_i) * np.abs(charge_j)
    
    #nuc-elec sum
    num += 2*np.square(num_elec)
    
    #elec-elec sum
    num += num_elec * (num_elec - 1)
    
    charge_ratio = num/denom
    
    return charge_ratio

def compute_charge_list(sim_parameters,molec_labels):
    """
    Compute charge ratio for molecule specified by parameter_file and psi4_out_file.
    
    Input:
        sim_parameters - Parameters for simulation of molecule
        molec_labels - List of labels for molecules
        
    Output:
        charge_list (array) - List of charges of nuclei in a given molecule
    """


    #Unpack data
    num_nuc = sim_parameters["nuclei"]



    charge_list = []

    for i in range(num_nuc):
        
        nucleon_i = molec_labels[i]
        charge_i = physical_consts.nuclear_charge(nucleon_i)
        charge_list.append(charge_i)
    
    return charge_list
