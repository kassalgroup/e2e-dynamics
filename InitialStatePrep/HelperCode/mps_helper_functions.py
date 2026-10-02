
import numpy as np


def get_pyscf_molecule(molec_geom,molec_labels):
    """
    Return a pyscf molecular geometry files.

    Input:
        molec_geom (num_atoms x 3) - numpy array containing equilibrium molecular geometry of molecule
        molec_labes (dict) - dictionary containing labels of molecule

    Output:
        pyscf_molec_string (string) - molecular geometry and labels formatted as a string that can be passed
                                      into pyscf 
    """
    
    # Extract values from the dictionary
    labels = np.array(list(molec_labels.values()))

    # Reshape the vector to make it a column vector (n x 1)
    vector_column = labels.reshape(-1, 1)
    
    # Concatenate the vector as the first column
    pyscf_molec = np.hstack((vector_column, molec_geom))

    pyscf_molec_string = ''
    
    for row in pyscf_molec:
        pyscf_molec_string += f"{row[0]} {row[1]} {row[2]} {row[3]}\n"
    
    return pyscf_molec_string

def get_param_file_extn(molec_dir,molec_tag,root_dir):
    """
    Create proper file extension for the parameter file of a given molecule.
    
    Input:
        molec_dir (string) - directory where info about the molecule is stored, 
                             i.e. .photodissociation/Moleculardata/molec_dir/
        molec_tag (string) - tag used before all file names within molec_dir,
                             i.e. .photodissociation/Moleculardata/molec_dir/molec_tag...
        root_dir (string) - /path/to/photodissociation/
        
    Output:
        parameter_extn (string) - root_dir/Moleculardata/molec_dir/molec_tag_parameters.in
    """
    
    parameter_extn = root_dir + "/MolecularData/" + molec_dir + "/" + molec_tag + "_parameters.in"
    
    return parameter_extn

def get_psi4_file_extn(molec_dir,molec_tag,root_dir):
    """
    Create proper file extension for the psi4 file of a given molecule.
    
    Input:
        molec_dir (string) - directory where info about the molecule is stored, 
                             i.e. .photodissociation/Moleculardata/molec_dir/
        molec_tag (string) - tag used before all file names within molec_dir,
                             i.e. .photodissociation/Moleculardata/molec_dir/molec_tag...
        root_dir (string) - /path/to/photodissociation/
                             
    Output:
        psi4_extn (string) - root_dir/Moleculardata/molec_dir/molec_tag_psi4.in
    """    
    
    psi4_extn = root_dir + "/MolecularData/" + molec_dir + "/Output/" + molec_tag + "_psi4.out"
    
    return psi4_extn

def get_output_file_extn(molec_dir, molec_tag, root_dir, particle):
    """
    Create proper file extension for the output file of a given molecule.
    
    Input:
        molec_dir (string) - directory where info about the molecule is stored, 
                             i.e. .photodissociation/Moleculardata/molec_dir/
        molec_tag (string) - tag used before all file names within molec_dir,
                             i.e. .photodissociation/Moleculardata/molec_dir/molec_tag...
        root_dir (string) - /path/to/photodissociation/  
        particle (string) - must be either "nuc" or "elec".
        
    Output:
        output_extn (string) - /Moleculardata/molec_dir/molec_tag_particle.in
    """    
    
    output_extn = root_dir + "/MolecularData/" + molec_dir + "/Output/" + molec_tag + "_" + particle + "_mps.out"
    
    return output_extn
    
def get_elec_mps_file_extn(molec_dir,molec_tag,root_dir):
    """
    Create proper file extension for the output file of a given molecule.
    
    Input:
        molec_dir (string) - directory where info about the molecule is stored, 
                             i.e. .photodissociation/Moleculardata/molec_dir/
        molec_tag (string) - tag used before all file names within molec_dir,
                             i.e. .photodissociation/Moleculardata/molec_dir/molec_tag...
        root_dir (string) - /path/to/photodissociation/  
        
    Output:
        output_extn (string) - /Moleculardata/molec_dir/molec_tag_elec_mps.in
    """    
    
    output_extn = root_dir + "/MolecularData/" + molec_dir + "/Output/" + molec_tag + "_elec_mps.out"
    
    return output_extn


def get_error_budget_file_extn(molec_dir,molec_tag,root_dir):
    """
    Create proper file extension for the error budget of a given molecule.

    Input:
        molec_dir (string) - directory where info about the molecule is stored,
                             i.e. .photodissociation/Moleculardata/molec_dir/
        molec_tag (string) - tag used before all file names within molec_dir,
                             i.e. .photodissociation/Moleculardata/molec_dir/molec_tag...
        root_dir (string) - /path/to/photodissociation/

    Output:
        error_budget_extn (string) - /Moleculardata/molec_dir/molec_tag_error_budget.in
    """

    error_budget_extn = root_dir + "/MolecularData/" + molec_dir + "/" + molec_tag + "_error_budget.in"

    return error_budget_extn


def get_log_file_extn(molec_dir, molec_tag, root_dir):
    """
    Create proper file extension for the error budget of a given molecule.

    Input:
        molec_dir (string) - directory where info about the molecule is stored,
                             i.e. .photodissociation/Moleculardata/molec_dir/
        molec_tag (string) - tag used before all file names within molec_dir,
                             i.e. .photodissociation/Moleculardata/molec_dir/molec_tag...
        root_dir (string) - /path/to/photodissociation/

    Output:
        log_file_extn (string) - /Moleculardata/molec_dir/elec_mps.log
    """

    log_file_extn = root_dir + "/MolecularData/" + molec_dir + "/Output/" + molec_tag + "_elec_mps.log"

    return log_file_extn



