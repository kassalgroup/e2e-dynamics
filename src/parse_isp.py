import numpy as np

def parse_isp_file(file_name):
    """
    Read in an ISP file and store/return the data in a dictionary.

    Input:
        file_name (string) - name of isp file that is to be read in

    Output:
        data_dict (dict) - dictionary containing isp data
    """

    # Initialize an empty dictionary to store the data
    data_dict = {}

    # Open the file and process it line by line
    with open(file_name, "r") as file:
        for line in file:
            # Strip whitespace and skip empty lines
            line = line.strip()
            if line:
                # Split the line into key and value
                key, value = line.split("=")
                try:
                    # Attempt to convert value to a float or int
                    value = float(value) if "." in value or "e" in value.lower() else int(value)
                except ValueError:
                    # Keep as string if conversion fails
                    pass
                # Add the key-value pair to the dictionary
                data_dict[key] = value

    return data_dict
