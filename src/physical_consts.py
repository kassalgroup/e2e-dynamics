import numpy as np

def convert_mass(mass):
    """
    Convert mass given in atomic mass units (AMU) to 
    atomic units (au). This corresponds to multiples of the mass
    of an electron
    
    Input:
        mass (float) - mass in AMU
    
    Output:
        mass_au (float) - the same input mass, now in au
    """
    
    electron_mass_AMU = 5.485799090441/np.power(10,4) 
    mass_au = mass/electron_mass_AMU
    
    return mass_au

def nuclear_mass(nucleon):
    """
    Return the average atomic mass of an element in atomic units (au).

    Input:
        nucleon (string) - label of nucleon

    Output:
        my_mass (float) - mass of nucleon in au
    
    """

    # Average atomic masses in AMU, using CIAAW's abridged 2024 weights:
    # https://www.ciaaw.org/abridged-atomic-weights.htm
    # Entries marked * have no standard atomic weight. Their values are
    # isotope mass numbers (approximate masses), from the IUPAC 2022 table:
    # https://iupac.org/wp-content/uploads/2022/07/IUPAC_Periodic_Table-04May22_CRA.pdf
    masses = {
        "H": 1.0080,
        "He": 4.0026,
        "Li": 6.94,
        "Be": 9.0122,
        "B": 10.81,
        "C": 12.011,
        "N": 14.007,
        "O": 15.999,
        "F": 18.998,
        "Ne": 20.180,
        "Na": 22.990,
        "Mg": 24.305,
        "Al": 26.982,
        "Si": 28.085,
        "P": 30.974,
        "S": 32.06,
        "Cl": 35.45,
        "Ar": 39.95,
        "K": 39.098,
        "Ca": 40.078,
        "Sc": 44.956,
        "Ti": 47.867,
        "V": 50.942,
        "Cr": 51.996,
        "Mn": 54.938,
        "Fe": 55.845,
        "Co": 58.933,
        "Ni": 58.693,
        "Cu": 63.546,
        "Zn": 65.38,
        "Ga": 69.723,
        "Ge": 72.630,
        "As": 74.922,
        "Se": 78.971,
        "Br": 79.904,
        "Kr": 83.798,
        "Rb": 85.468,
        "Sr": 87.62,
        "Y": 88.906,
        "Zr": 91.222,
        "Nb": 92.906,
        "Mo": 95.95,
        "Tc": 97.0,  # *
        "Ru": 101.07,
        "Rh": 102.91,
        "Pd": 106.42,
        "Ag": 107.87,
        "Cd": 112.41,
        "In": 114.82,
        "Sn": 118.71,
        "Sb": 121.76,
        "Te": 127.60,
        "I": 126.90,
        "Xe": 131.29,
        "Cs": 132.91,
        "Ba": 137.33,
        "La": 138.91,
        "Ce": 140.12,
        "Pr": 140.91,
        "Nd": 144.24,
        "Pm": 145.0,  # *
        "Sm": 150.36,
        "Eu": 151.96,
        "Gd": 157.25,
        "Tb": 158.93,
        "Dy": 162.50,
        "Ho": 164.93,
        "Er": 167.26,
        "Tm": 168.93,
        "Yb": 173.05,
        "Lu": 174.97,
        "Hf": 178.49,
        "Ta": 180.95,
        "W": 183.84,
        "Re": 186.21,
        "Os": 190.23,
        "Ir": 192.22,
        "Pt": 195.08,
        "Au": 196.97,
        "Hg": 200.59,
        "Tl": 204.38,
        "Pb": 207.2,
        "Bi": 208.98,
        "Po": 209.0,  # *
        "At": 210.0,  # *
        "Rn": 222.0,  # *
        "Fr": 223.0,  # *
        "Ra": 226.0,  # *
        "Ac": 227.0,  # *
        "Th": 232.04,
        "Pa": 231.04,
        "U": 238.03,
        "Np": 237.0,  # *
        "Pu": 244.0,  # *
        "Am": 243.0,  # *
        "Cm": 247.0,  # *
        "Bk": 247.0,  # *
        "Cf": 251.0,  # *
        "Es": 252.0,  # *
        "Fm": 257.0,  # *
        "Md": 258.0,  # *
        "No": 259.0,  # *
        "Lr": 262.0,  # *
        "Rf": 267.0,  # *
        "Db": 268.0,  # *
        "Sg": 269.0,  # *
        "Bh": 270.0,  # *
        "Hs": 269.0,  # *
        "Mt": 277.0,  # *
        "Ds": 281.0,  # *
        "Rg": 282.0,  # *
        "Cn": 285.0,  # *
        "Nh": 286.0,  # *
        "Fl": 290.0,  # *
        "Mc": 290.0,  # *
        "Lv": 293.0,  # *
        "Ts": 294.0,  # *
        "Og": 294.0,  # *
    }

    return convert_mass(masses[nucleon.capitalize()])


def nuclear_charge(nucleon):
    """
    Return the nuclear charge of an element in atomic units (au).

    Input:
        nucleon (string) - label of nucleon

    Output:
        my_charge (float) - charge of nucleon in au
    
    """

    # Nuclear charges in atomic units (elementary charges), in atomic-number order.
    charges = {
        "H": 1.0,
        "He": 2.0,
        "Li": 3.0,
        "Be": 4.0,
        "B": 5.0,
        "C": 6.0,
        "N": 7.0,
        "O": 8.0,
        "F": 9.0,
        "Ne": 10.0,
        "Na": 11.0,
        "Mg": 12.0,
        "Al": 13.0,
        "Si": 14.0,
        "P": 15.0,
        "S": 16.0,
        "Cl": 17.0,
        "Ar": 18.0,
        "K": 19.0,
        "Ca": 20.0,
        "Sc": 21.0,
        "Ti": 22.0,
        "V": 23.0,
        "Cr": 24.0,
        "Mn": 25.0,
        "Fe": 26.0,
        "Co": 27.0,
        "Ni": 28.0,
        "Cu": 29.0,
        "Zn": 30.0,
        "Ga": 31.0,
        "Ge": 32.0,
        "As": 33.0,
        "Se": 34.0,
        "Br": 35.0,
        "Kr": 36.0,
        "Rb": 37.0,
        "Sr": 38.0,
        "Y": 39.0,
        "Zr": 40.0,
        "Nb": 41.0,
        "Mo": 42.0,
        "Tc": 43.0,
        "Ru": 44.0,
        "Rh": 45.0,
        "Pd": 46.0,
        "Ag": 47.0,
        "Cd": 48.0,
        "In": 49.0,
        "Sn": 50.0,
        "Sb": 51.0,
        "Te": 52.0,
        "I": 53.0,
        "Xe": 54.0,
        "Cs": 55.0,
        "Ba": 56.0,
        "La": 57.0,
        "Ce": 58.0,
        "Pr": 59.0,
        "Nd": 60.0,
        "Pm": 61.0,
        "Sm": 62.0,
        "Eu": 63.0,
        "Gd": 64.0,
        "Tb": 65.0,
        "Dy": 66.0,
        "Ho": 67.0,
        "Er": 68.0,
        "Tm": 69.0,
        "Yb": 70.0,
        "Lu": 71.0,
        "Hf": 72.0,
        "Ta": 73.0,
        "W": 74.0,
        "Re": 75.0,
        "Os": 76.0,
        "Ir": 77.0,
        "Pt": 78.0,
        "Au": 79.0,
        "Hg": 80.0,
        "Tl": 81.0,
        "Pb": 82.0,
        "Bi": 83.0,
        "Po": 84.0,
        "At": 85.0,
        "Rn": 86.0,
        "Fr": 87.0,
        "Ra": 88.0,
        "Ac": 89.0,
        "Th": 90.0,
        "Pa": 91.0,
        "U": 92.0,
        "Np": 93.0,
        "Pu": 94.0,
        "Am": 95.0,
        "Cm": 96.0,
        "Bk": 97.0,
        "Cf": 98.0,
        "Es": 99.0,
        "Fm": 100.0,
        "Md": 101.0,
        "No": 102.0,
        "Lr": 103.0,
        "Rf": 104.0,
        "Db": 105.0,
        "Sg": 106.0,
        "Bh": 107.0,
        "Hs": 108.0,
        "Mt": 109.0,
        "Ds": 110.0,
        "Rg": 111.0,
        "Cn": 112.0,
        "Nh": 113.0,
        "Fl": 114.0,
        "Mc": 115.0,
        "Lv": 116.0,
        "Ts": 117.0,
        "Og": 118.0,
    }

    return charges[nucleon.capitalize()]

def num_electrons(molecule):
    """
    Return the number of electrons (or total nuclear charge in a neutral molecule) for a given molecule.

    Input: 
        molecule (1d array) - array containing labels of nuclei in molecule

    Output:
        n_elec (float) - number of electrons in molecule
    """
    n_elec = 0

    for nucleon in molecule:
        charge = nuclear_charge(nucleon)
        n_elec += charge

    return n_elec


def num_nuclei(molecule):
    """
    Return the number of nuclei for a given molecule.

    Input: 
        molecule (1d array) - array containing labels of nuclei in molecule

    Output:
        num_nuc (float) - number of nuclei in molecule
    """

    num_nuc = int(len(molecule))

    return num_nuc

def femto2au(t):
    """
    Convert time in femtosecond to atomic units.

    Input:
        t (float) - time in femtoseconds

    Output:
        time_au (float) - time in atomic units
    """

    femto2sec = 1e15
    au_time = 2.4188843265864/np.power(10,17)
    time_au = (t/femto2sec)/au_time
    
    return time_au

def hartree2wavenum(energy):
    """
    Convert an energy from au to wavenumbers.

    Input:
        energy (float) - energy in au

    Output:
        energy_wn (float) - energy in wavenumbers.
    """
    
    #speed of light (sol) in cm/second
    sol = 2.998*np.power(10,10) 

    #conversion factor for frequency in au to SI units (1/seconds)
    au_freq_SI = 4.13 * np.power(10,16) 

    energy_wn = energy * au_freq_SI/(sol * 2 * np.pi)

    return energy_wn

def wavenum2hartree(energy_wn):
    """
    Convert an energy wavenumbers to au.

    Input:
        energy_wn (float) - energy in wavenumbers

    Output:
        energy (float) - energy in au.
    """
    
    #speed of light (sol) in cm/second
    sol = 2.998*np.power(10,10) 

    #conversion factor for frequency in au to SI units (1/seconds)
    au_freq_SI = 4.13 * np.power(10,16) 

    energy = energy_wn * (sol * 2 * np.pi) / au_freq_SI

    return energy

def angstrom2bohr(l):
    """
    Convert a length given in Angstroms to bohr.
    
    Input:
        l (float) - length in angstroms
        
    Output:
        l_bohr (float) - length in bohrs
    """
    
    l_bohr = l / 0.529177
    
    return l_bohr
