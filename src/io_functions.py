import numpy as np

def read_file(file_name):
    data_dict = {}

    with open(file_name, "r") as file:
        for line in file:
            line = line.strip()
            if not line or "=" not in line:
                continue

            key, value_str = line.split("=", 1)
            key = key.strip()
            value_str = value_str.strip()

            # Check if comma-separated list
            if "," in value_str:
                items = value_str.split(",")
                parsed_items = []
                for item in items:
                    item = item.strip()
                    try:
                        parsed_item = float(item) if "." in item or "e" in item.lower() else int(item)
                    except ValueError:
                        parsed_item = item
                    parsed_items.append(parsed_item)
                data_dict[key] = parsed_items
            else:
                try:
                    parsed_value = float(value_str) if "." in value_str or "e" in value_str.lower() else int(value_str)
                except ValueError:
                    parsed_value = value_str
                data_dict[key] = parsed_value

    return data_dict

def get_param_file_extn(molec_dir,molec_tag):
    """
    Create proper file extension for the parameter file of a given molecule.
    
    Input:
        molec_dir (string) - directory where info about the molecule is stored, 
                             i.e. .photodissociation/Moleculardata/molec_dir/
        molec_tag (string) - tag used before all file names within molec_dir,
                             i.e. .photodissociation/Moleculardata/molec_dir/molec_tag...
                             
    Output:
        parameter_extn (string) - /Moleculardata/molec_dir/molec_tag_parameters.in
    """
    
    parameter_extn = "/MolecularData/" + molec_dir + "/" + molec_tag + "_parameters.in"
    
    return parameter_extn

def get_psi4_file_extn(molec_dir,molec_tag):
    """
    Create proper file extension for the psi4 file of a given molecule.
    
    Input:
        molec_dir (string) - directory where info about the molecule is stored, 
                             i.e. .photodissociation/Moleculardata/molec_dir/
        molec_tag (string) - tag used before all file names within molec_dir,
                             i.e. .photodissociation/Moleculardata/molec_dir/molec_tag...
                             
    Output:
        psi4_extn (string) - /Moleculardata/molec_dir/molec_tag_psi4.in
    """    
    
    psi4_extn = "/MolecularData/" + molec_dir + "/" + molec_tag + "_psi4.out"
    
    return psi4_extn
