import numpy as np

def find_last_instance(file_name, search_phrase):
    """
    Find and return the line number corresponding to the last instance of search_phrase
    in file_name.

    Input:
        file_name (string) - name of the file that is to be searched
        search_phrase (string) - phrase that is being searched

    Output:
        line_number (int) - line number that search_phrase occurs at
    """

    line_number = None
    line_content = None

    with open(file_name, 'r') as file:
        for num, line in enumerate(file, start=1):
            if search_phrase in line:
                line_number = num
                line_content = line.strip()
    
    if line_number is not None:
        return line_number
    else:
        print(f"'{search_phrase}' not found in the file.")
        return None


def read_lines_range(file_name, start_line, stop_line):
    """
    Read in the lines starting at start_line and ending at stop_line
    of my_file. Store data in an array and return the array.

    Input:
        file_name (string) - file name
        start_line (int) - starting line of file that is to be read in
        stop_line (int) - stopping line of file that is to be read in

    Output:
        my_data (array) - array containing all the lines read in from file_name
    """
    
    my_data = []
    
    with open(file_name, 'r') as file:
        for line_number, line in enumerate(file, start=1):
            if start_line <= line_number <= stop_line:
                my_data.append(line.split())
    
    return my_data

def prase_hess_data(raw_data,num_atoms):
    """
    Parse Hessian data and return as a numpy array.

    Input:
        raw_data (array) - array of Hessian data that is to be parsed
        num_atoms (int) - number of atoms in the molecule

    Output:
        hess_np (array) - numpy array of Hessian data
    """

    #Convert each row of string data to a row of float data
    clean_data = []
    dof = 3 * num_atoms
    hess_np = np.zeros((dof,dof))
    
    for n ,row in enumerate(raw_data):
        float_row = [float(fl_row) for fl_row in row]
        
        if np.mod(n,dof+3)-3 >= 0:
            clean_data.append(float_row)
            
    #Remove first element
    for n ,row in enumerate(clean_data):
        row.pop(0)

    for n_row in range(dof):
        for n_col in range(dof):
            cd_col = int(np.mod(n_col,5))
            cd_row = int(n_row + (np.floor(n_col/5)*dof))
            clean_data[cd_row][cd_col]
            hess_np[n_row,n_col] = clean_data[cd_row][cd_col]
        
    return hess_np

def parse_molecular_data(dirty_data,num_atoms):
    """
    Parse molecular geometry data and return molecular massess, labels, and geometry

    Input:
        dirty_data (array) - array containing relevant molecular data
        num_atoms (int) - number of atoms in the molecule
    Output:
        molec_geometry (num_atoms x 3 array) - array of molecular geometry
        molec_masses (dict) - dictionary of molecular masses
        molec_labels (dict) - dictionary of labels (atomic symbol) of each atom in molecule
    """

    molec_geometry = np.zeros((num_atoms,3))
    molec_masses = {}
    molec_labels = {}

    for i in range(num_atoms):

        label_i = dirty_data[i][0]
        mass_i = float(dirty_data[i][4])
        
        molec_labels[i] = label_i        
        molec_masses[i] = mass_i

        for j in range(3):
            molec_geometry[i,j] = float(dirty_data[i][1+j])

    return molec_geometry, molec_masses, molec_labels

def get_hessian(file_name,num_atoms):
    """
    Return the hessian from a psi4.out file.

    Input:
        file_name (string) - psi4.out file with hessian data
        num_atoms (int) - number of atoms in molecule

    Output:
        my_hess (nxn array) - hessian from file_name 
    """
    
    #Define stride parameters
    num_dof = 3 * num_atoms 
    num_big_rows = int(np.ceil(num_dof/5))
    big_row_size = num_dof + 3

    #Find hessian data in file
    search_phrase = "## New Matrix (Symmetry 0) ##"
    line_num = find_last_instance(file_name, search_phrase)
    start_line = line_num + 3
    stop_line = start_line + (num_big_rows * big_row_size) - 1

    #Extract hesian data
    dirty_data = read_lines_range(file_name, start_line, stop_line)
    dirty_data.insert(0,[])

    #Parse hessian data
    my_hess = prase_hess_data(dirty_data,num_atoms)

    return my_hess

def get_molec_data(file_name,num_atoms):
    """
    Return molecular data (geometry, labels, and masses) from a psi4.out file.

    Input:
        file_name (string) - psi4.out file with molecular data
        num_atoms (int) - number of atoms in molecule

    Output:
        molec_geometry (num_atoms x 3 array) - array of molecular geometry
        molec_masses (dict) - dictionary of molecular masses
        molec_labels (dict) - dictionary of labels (atomic symbol) of each atom in molecule
    """

    #Find molecular data in file
    search_phrase = "Geometry (in Bohr),"
    line_number = find_last_instance(file_name, search_phrase)
    start_line = line_number + 4
    stop_line = start_line + num_atoms

    #Extract molecular data
    dirty_data = read_lines_range(file_name, start_line, stop_line)

    #Parse molecular data
    molec_geom, molec_masses, molec_labels = parse_molecular_data(dirty_data,num_atoms)

    return molec_geom, molec_masses, molec_labels
