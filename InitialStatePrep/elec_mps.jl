using PyCall
using DelimitedFiles
using Printf
using LinearAlgebra

BLAS.set_num_threads(1)

# Hard-coded parameters
ds_trunc = 1e-20
final_trunc = 1e-12
rot_qubits = 70 #rotational qubits used for unitary synthesis

verbose = false

molec_dir = ARGS[1]
molec_tag = ARGS[2]

# Directory containing elec_mps.jl
current_dir = @__DIR__
root_dir = dirname(current_dir)

# Helper-code directories
helper_dir = joinpath(current_dir, "HelperCode")
src_dir = joinpath(root_dir, "src")

# Add Python helper directories
sys_path = pyimport("sys")."path"
sys_path.insert(0, helper_dir)
sys_path.insert(1, src_dir)

# Include Julia helper files
include(joinpath(helper_dir, "elec_libfqmps.jl"))
include(joinpath(helper_dir, "resource_estimates.jl"))

# Load Python modules
np = pyimport("numpy")
scipy = pyimport("scipy")
pf = pyimport("pyscf_funcs")
pyscf = pyimport("pyscf")

gto = pyscf["gto"]
scf = pyscf["scf"]

special = pyimport("scipy.special")
factorial2 = special.factorial2

psi4_parse = pyimport("parse_psi4")
io = pyimport("io_functions")
helper_fxn = pyimport("mps_helper_functions")

psi4_file = helper_fxn.get_psi4_file_extn(
    molec_dir,
    molec_tag,
    root_dir,
)

params_file = helper_fxn.get_param_file_extn(
    molec_dir,
    molec_tag,
    root_dir,
)

output_file = helper_fxn.get_output_file_extn(
    molec_dir,
    molec_tag,
    root_dir,
    "elec",
)

error_budget_file = helper_fxn.get_error_budget_file_extn(
    molec_dir,
    molec_tag,
    root_dir,
)


log_file = helper_fxn.get_log_file_extn(
    molec_dir, 
    molec_tag, 
    root_dir
)

log_io = open(log_file, "w")
redirect_stdout(log_io)
redirect_stderr(log_io)

println(molec_dir)
println(molec_tag)
flush(stdout)

# Read data from files
parameters = io.read_file(params_file)
error_budget = io.read_file(error_budget_file)

num_atoms = parameters["nuclei"]
L = parameters["L_guess"]
E = parameters["E_guess"]
error_tol = error_budget["elec_isp_budget"]

molec_geom, molec_masses, molec_labels =
    psi4_parse.get_molec_data(psi4_file, num_atoms)

pyscf_molec = helper_fxn.get_pyscf_molecule(
    molec_geom,
    molec_labels,
)

# Unpack parameters and set values
basis = parameters["basis"]
spin = parameters["spin"]
num_electron = parameters["electrons"]

homo_num = Int(floor(num_electron / 2))
lumo_num = homo_num + 1


# Run RHF in PySCF and obtain data
a_vec, l_vec, gamma_vec, C_prim, norm_vec, G, nmo =
    pf.run_pyscf(pyscf_molec, basis, spin)

println("Number of MOs found: $(nmo)")

max_bond_dim = 0
state_error_total = 0.0
rot_error_total = 0.0
toff_gate_sum = 0

# Set lower bound on allowed error tolerance
error_tol_delta = error_tol / 10
error_tol_lower_bound = 9.5 * error_tol_delta

current_error = 1.0

# Set minimum and maximum energies to search through
min_E = E
max_E = 7.0 * E

total_error = 0
T_gates = 0
rot_error_total = 0
toff_gate_sum = 0
state_error_total = 0
current_E = E

println("Error tol: $(error_tol)")
println("Min Energy: $(min_E)")
println("Max Energy: $(max_E)")
println(" ")

# Outer loop over error tolerance
counter = 0

while counter < 1
    global current_E
    global min_E
    global max_E
    global current_error
    global counter
    global max_bond_dim
    global total_error
    global rot_error_total
    global toff_gate_sum
    global state_error_total
    global T_gates

    counter += 1

    # Reset after each iteration
    rot_error_total = 0
    state_error_total = 0
    toff_gate_sum = 0

    println("Iteration Error: $(counter)")
    println("Trying with energy: $(current_E)")

    # Preallocate arrays for thread-safe accumulation
    toff_gate_sums = zeros(Float64, lumo_num)
    rot_errors = zeros(Float64, lumo_num)
    state_errors = zeros(Float64, lumo_num)
    bond_dims = zeros(Int, lumo_num)

    # Loop over molecular orbitals
    Threads.@threads for mo_number = 1:lumo_num
        # Choose MO number: column of C_prim
        c_col = C_prim[:, mo_number]

        println(
            "Thread $(Threads.threadid()) working on MO #$(mo_number)"
        )

        # Get the MPS representation of the selected Gaussian-basis
        # molecular orbital in the plane-wave basis.
        MPS_mo, N_pw_local, num_qubits_dim_local, sv_fidelity =
            StateVecMPS(
                a_vec,
                l_vec,
                gamma_vec,
                c_col;
                L=L,
                E=current_E,
                verbose=verbose,
                ds_trunc=ds_trunc,
                final_trunc=final_trunc,
            )

        println(sv_fidelity)

        bond_dim = maximum(linkdims(MPS_mo))
        toff_gates = estimate_toff_gates(rot_qubits, MPS_mo)
        rot_error = error_rotation(rot_qubits, MPS_mo)
        trace_dist = TraceDist(sv_fidelity)

        toff_gate_sums[mo_number] = toff_gates

        rot_errors[mo_number] =
            mo_number < homo_num ? 2.0 * rot_error : rot_error

        state_errors[mo_number] =
            mo_number < homo_num ? 2.0 * trace_dist : trace_dist

        bond_dims[mo_number] = bond_dim
    end

    toff_gate_sum = sum(toff_gate_sums)
    rot_error_total = sum(rot_errors)
    state_error_total = sum(state_errors)
    max_bond_dim = maximum(bond_dims)

    total_error = rot_error_total + state_error_total
    current_error = total_error

    # Get the number of qubits per electron per dimension
    K_max = sqrt(2.0 * current_E)
    p_max = L * K_max / (2.0 * pi)
    N_pw_init = 2 * p_max + 1
    n_qb = ceil(Int, log2(N_pw_init))

    # Find the largest L for the given n_qb and E
    p_max_actual = (2^n_qb - 1) / 2
    L_opt = p_max_actual * pi * sqrt(2) / sqrt(current_E)

    T_gates = total_T_gates(
        n_qb,
        num_electron,
        toff_gate_sum,
    )

    println("Total State Error: $(state_error_total)")
    println("Total Rotational Error: $(rot_error_total)")
    println("Total Error: $(total_error)")
    println("Total T gates: $(T_gates)")
    println("Maximum Bond Dimension: $(max_bond_dim)")
    println(" ")

     # Write data to file
     open(output_file, "w") do io
         println(io, "elec_isp_eps=", total_error)
         println(io, "E_elec=", current_E)
         println(io, "L_elec=", L_opt)
         println(io, "n_qubits_elec=", n_qb)
         println(io, "elec_bond_dim_max=", max_bond_dim)
     end

    if current_error > error_tol
        println("The computed error is larger than the error budget. Try increasing the energy to lower the error.")
        min_E = current_E
    end

    if current_error < error_tol_lower_bound
        max_E = current_E
        println("The computed error is smaller than the error budget. Try decreasing the energy to increase the error.")
    end

    if counter > 1
        println("Exceeded max iterations")
        current_error = error_tol
    end

    println(" ")
    println(" ")
end

println("Finished")
flush(stdout)
flush(stderr)
