using ITensorMPS
using ITensors
using LinearAlgebra
using Random
using Polynomials
using SpecialPolynomials

function bond_bar(Mmo)
    """
    Equation F86 of Bill's paper.

    Input:
        Mmo (mps) - matrix product state?

    Output: 
        m_bar (int) - Eq. F86 evaluated for Mmo
    """
         
    bond_dims = linkdims(Mmo)
    m_bar = zeros(length(bond_dims))
    
    m_bar[1] = exp2(ceil(log2(bond_dims[1])))
    
    for j in 2:length(bond_dims)
        mj = bond_dims[j]
        mj_m1 =  bond_dims[j-1]
        
        vj = exp2(ceil(log2(mj)))
        vj_m1 = exp2(ceil(log2(mj_m1)))
        
        m_bar[j] = max(vj,vj_m1)
    end
    
    return m_bar
end

function error_rotation(b,Mmo)
    """
    Compute upper bound to error caused from rotations in multiplexors.

    Input:
        b (int) - number of bit of precision used for multiplexors
        Mmo (mps) - matrix product state?
    Output:
        error_final (float) - error resulting from multiplexors
    """

    bond_dims = linkdims(Mmo)
    m_bar = bond_bar(Mmo)
    error_tot = 0.0
    
    for j in 1:length(bond_dims)
        mj = bond_dims[j]
        m_bar_j = m_bar[j]
        
        v = mj*log2(2*m_bar_j)
        error_tot+= v
    end
    
    error_final = exp2((3.5)-b)*error_tot
    return error_final
end

function estimate_toff_gates(b,Mmo)
    """
    Compute estimated cost of preparing Mmo (MPS) in terms of Toffoli gates. 
    This uses Eq. F94 of Bill's paper.

    Input:
        b (int) - number of bit of precision used for multiplexors
        Mmo (mps) - matrix product state?

    Output:
        toff_gates (float) - number of Toffoli gates needed
    """
    
    toff_gates = 0
    bond_dims = linkdims(Mmo)
    bond_dims_bar = bond_bar(Mmo)
    
    pre_fac_1 = 32*(1 + np.sqrt(2))*np.sqrt(b+1)
    pre_fac_2 = (8*b-15)
    
    for k in 1:length(bond_dims)
        
        term_1 = pre_fac_1 * bond_dims[k] * np.sqrt(bond_dims_bar[k]) 
        term_2 = pre_fac_2 * bond_dims[k] * np.log2(2*bond_dims_bar[k])
        term_k = term_1 + term_2
        toff_gates += term_k
    end
    return toff_gates
end

function total_T_gates(num_qubits_dim, num_electron, toff_gate_sum)
    """
    Compute total cost of preparing electronic state in terms of T gates. 
    This uses Eq. 122 of Bill's paper.

    Input:
        num_qubits_dim (int) - number of bit of precision used for multiplexors
        num_electrons (int) - number of electrons
        toff_gate_sum (int) - sum of toffoli gates for preparing individual MPS

    Output:
        T_gate_tot (int) - total cost of preparing electronic state
    """    
    
    term1 = 3 * num_qubits_dim * num_electron^2
    term2 = 2 * num_electron * toff_gate_sum

    T_gate_tot = ceil(4 * (term1 + term2))

    return T_gate_tot
end
