
def get_param_file_extn(molec_name):
    """Construct the path suffix for a molecule's parameter input file.

    Parameters
    ----------
    molec_name : str
        Molecule name. The original spelling is used as the directory name and
        its lowercase form is used as the filename tag.

    Returns
    -------
    str
        Path suffix of the form
        ``/MolecularData/<molecule>/<tag>_parameters.in``.
    """
    molec_dir = molec_name
    molec_tag = molec_dir.lower()

    parameter_extn = "/MolecularData/" + molec_dir + "/" + molec_tag + "_parameters.in"

    return parameter_extn

def get_psi4_file_extn(molec_name):
    """Construct the path suffix for a molecule's Psi4 output file.

    Parameters
    ----------
    molec_name : str
        Molecule name. The original spelling is used as the directory name and
        its lowercase form is used as the filename tag.

    Returns
    -------
    str
        Path suffix of the form
        ``/MolecularData/<molecule>/Output/<tag>_psi4.out``.
    """
    molec_dir = molec_name
    molec_tag = molec_dir.lower()

    psi4_extn = "/MolecularData/" + molec_dir + "/Output/" + molec_tag + "_psi4.out"

    return psi4_extn

def get_output_file_extn(molec_name):
    """Construct the path suffix for a molecule's final ISP output file.

    Parameters
    ----------
    molec_name : str
        Molecule name. The original spelling is used as the directory name and
        its lowercase form is used as the filename tag.

    Returns
    -------
    str
        Path suffix of the form
        ``/MolecularData/<molecule>/Output/<tag>_isp_final_data.out``.
    """
    molec_dir = molec_name
    molec_tag = molec_dir.lower()

    output_extn = "/MolecularData/" + molec_dir + "/Output/" + molec_tag + "_isp_final_data.out"

    return output_extn

def get_elec_mps_file_extn(molec_name):
    """Construct the path suffix for a molecule's electronic-MPS output.

    Parameters
    ----------
    molec_name : str
        Molecule name. The original spelling is used as the directory name and
        its lowercase form is used as the filename tag.

    Returns
    -------
    str
        Path suffix of the form
        ``/MolecularData/<molecule>/Output/<tag>_elec_mps.out``.
    """
    molec_dir = molec_name
    molec_tag = molec_dir.lower()

    elec_mps = "/MolecularData/" + molec_dir + "/Output/" + molec_tag + "_elec_mps.out"

    return elec_mps

def get_error_budget_extn(molec_name):
    """Construct the path suffix for a molecule's error-budget input file.

    Parameters
    ----------
    molec_name : str
        Molecule name. The original spelling is used as the directory name and
        its lowercase form is used as the filename tag.

    Returns
    -------
    str
        Path suffix of the form
        ``/MolecularData/<molecule>/<tag>_error_budget.in``.
    """
    molec_dir = molec_name
    molec_tag = molec_dir.lower()

    error_budget_extn = "/MolecularData/" + molec_dir + "/" + molec_tag + "_error_budget.in"

    return error_budget_extn

def get_error_lambda_mat_extn(molec_name):
    """Construct the path suffix for a molecule's lambda-matrix output.

    Parameters
    ----------
    molec_name : str
        Molecule name. The original spelling is used as the directory name and
        its lowercase form is used as the filename tag.

    Returns
    -------
    str
        Path suffix of the form
        ``/MolecularData/<molecule>/<tag>_lambda_mat.out``.
    """
    molec_dir = molec_name
    molec_tag = molec_dir.lower()

    lambda_mat_extn = "/MolecularData/" + molec_dir + "/Output/" + molec_tag + "_lambda_mat.out"

    return lambda_mat_extn
