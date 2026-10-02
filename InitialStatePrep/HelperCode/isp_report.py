from typing import Any, Optional

def _format_section(
    title: str,
    rows: list[tuple[str, Any, Optional[str]]],
) -> str:
    """Format a titled section of aligned report rows.

    Parameters
    ----------
    title : str
        Heading displayed above the section.
    rows : list of tuple of (str, Any, str or None)
        Row labels, values, and optional format specifications.

    Returns
    -------
    str
        Formatted section containing the heading and aligned rows.
    """

    label_width = max(len(label) for label, _, _ in rows)

    lines = [
        title,
        "-" * len(title),
    ]

    for label, value, format_spec in rows:
        formatted_value = (
            format(value, format_spec)
            if format_spec is not None
            else str(value)
        )

        lines.append(f"{label:<{label_width}} : {formatted_value}")

    return "\n".join(lines)


def format_isp_summary(
    *,
    L: float,
    nuc_k_max: float,
    n_p: int,
    nuc_data_qubits: int,
    delta_L: float,
    nuc_bond_dim_max: int,
    W_nuc_Toff: float,
    W_nuc_anc: int,
    U_A_inv_trans_Toff: float,
    U_A_inv_trans_anc: int,
    PK_Toff: float,
    PK_anc: int,
    elec_k_max: float,
    elec_data_qubits: int,
    ASYM_toff: float,
    ASYM_anc: int,
    elec_bond_dim_max: int,
    W_elec_toff: float,
    W_elec_anc: int,
    isp_toff: float,
    isp_anc: int,
    isp_data: int,
) -> str:
    """Build a formatted initial-state-preparation resource summary.

    Parameters
    ----------
    L : float
        Final grid length.
    nuc_k_max : float
        Maximum nuclear momentum.
    n_p : int
        Number of qubits per spatial coordinate.
    nuc_data_qubits : int
        Number of nuclear data qubits.
    delta_L : float
        Nuclear grid spacing, ``L / N``.
    nuc_bond_dim_max : int
        Maximum nuclear MPS bond dimension.
    W_nuc_Toff : float
        Toffoli cost of nuclear MPS preparation.
    W_nuc_anc : int
        Ancilla cost of nuclear MPS preparation.
    U_A_inv_trans_Toff : float
        Toffoli cost of applying the inverse-transpose transformation.
    U_A_inv_trans_anc : int
        Ancilla cost of applying the inverse-transpose transformation.
    PK_Toff : float
        Toffoli cost of phase kickback.
    PK_anc : int
        Ancilla cost of phase kickback.
    elec_k_max : float
        Maximum electronic momentum.
    elec_data_qubits : int
        Number of electronic data qubits.
    ASYM_toff : float
        Toffoli cost of antisymmetrization.
    ASYM_anc : int
        Ancilla cost of antisymmetrization.
    elec_bond_dim_max : int
        Maximum electronic MPS bond dimension.
    W_elec_toff : float
        Toffoli cost of electronic MPS preparation.
    W_elec_anc : int
        Ancilla cost of electronic MPS preparation.
    isp_toff : float
        Total initial-state-preparation Toffoli cost.
    isp_anc : int
        Total initial-state-preparation ancilla cost.
    isp_data : int
        Total number of initial-state-preparation data qubits.

    Returns
    -------
    str
        Multi-section plain-text resource summary.
    """

    nuclear_rows = [
        ("Length", L, ".2e"),
        ("Nuclear K_max", nuc_k_max, ".2f"),
        ("n_p", n_p, None),
        ("Nuclear data qubits", nuc_data_qubits, None),
        ("Delta L (L/N)", delta_L, ".4f"),
        ("Nuclear bond dim", nuc_bond_dim_max, None),
        ("W_n Toff", W_nuc_Toff, ".2e"),
        ("W_n anc", W_nuc_anc, None),
        ("U_A^{-T} Toff", U_A_inv_trans_Toff, ".1e"),
        ("U_A^{-T} anc", U_A_inv_trans_anc, None),
        ("PK Toff", PK_Toff, ".1e"),
        ("PK anc", PK_anc, None),
    ]

    electronic_rows = [
        ("Electronic K_max", elec_k_max, ".2f"),
        ("Electronic data qubits", elec_data_qubits, None),
        ("ASYM Toff", ASYM_toff, ".1e"),
        ("ASYM anc", ASYM_anc, None),
        ("Electronic bond dim", elec_bond_dim_max, None),
        ("W_elec Toff", W_elec_toff, ".2e"),
        ("W_elec anc", W_elec_anc, None),
    ]

    total_rows = [
        ("Toffolis", isp_toff, ".2e"),
        ("Ancilla", isp_anc, None),
        ("Data qubits", isp_data, None),
    ]

    sections = [
        _format_section(
            "Nuclear Initial State Preparation",
            nuclear_rows,
        ),
        _format_section(
            "Electronic Initial State Preparation",
            electronic_rows,
        ),
        _format_section(
            "Total Initial State Preparation",
            total_rows,
        ),
    ]

    return "\n\n".join(sections)
