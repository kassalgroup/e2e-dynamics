# Error Budget Breakdown

Below is a breakdown of the total simulation error into its fundamental contributions, following [1]. <br>
Errors shown in <span style="color:blue">blue</span> correspond to error tolerances specified by the user in the `_error_budget.in` file. <br>
Errors shown in <span style="color:green">green</span> are hard-coded. These errors can typically be exponentially supressed, and we have set them to very small values. <br>
The remaining errors are solved for such that all errors sum to $\epsilon$.

$$
\color{blue}{\epsilon}
\le
2\left(
\color{blue}{\epsilon_{\mathrm{ISP}}^{(\mathrm{sep})}}
+
\epsilon_{\mathrm{prop}}
+
\color{blue}{\epsilon_B}
\right)
+
\epsilon_{\mathrm{meas}}.
$$

- **Initial-state preparation:** 
$\color{blue}{\epsilon_{\mathrm{ISP}}^{(\mathrm{sep})}} 
\leq 
\color{blue}{\epsilon_{\mathrm{ISP}}^{(e)}} + 
\epsilon_{\mathrm{ISP}}^{(n)}$

    - **Electronic ISP:** 
    $\color{blue}{\epsilon_{\mathrm{ISP}}^{(e)} }
    \leq 
    \color{green}{\epsilon_{\mathrm{ASP}}^{(e)}} + \sum_{a=1}^{N_{\mathrm{MOB}}}\epsilon_a^{(e)}$
    
        - $\color{green}{\epsilon_{\mathrm{ASP}}^{(e)}}$: electronic arbitrary-state-preparation error
        
        - $\epsilon_a^{(e)} \leq \epsilon_a^{(e,c)} + \color{green}{\epsilon_a^{(e,q)}}$
            - $\epsilon_a^{(e,c)}$: MO preprocessing / continuous approximation error
            - $\color{green}{\epsilon_a^{(e,q)}}$: unitary-synthesis / MPS-preparation error
            
    - **Nuclear ISP:**
     $\epsilon_{\mathrm{ISP}}^{(n)} \leq \color{green}{\epsilon_{\mathrm{ASP}}^{(n)}} +
    2^{3/2}\sum_{i=1}^{3\eta_n} \bigg(\epsilon_{i0}^{(n,c)} +  \color{green}{\epsilon_{i0}^{(n,q)}}\bigg) + \epsilon_{\mathrm{shear}} + \color{blue}{\epsilon_{\mathrm{PK}}} + \color{blue}{\epsilon_{\mathrm{trim}}} $

        - $\color{green}{\epsilon_{\mathrm{ASP}}^{(n)}}$: nuclear arbitrary-state-preparation error
        - $\epsilon_{i0}^{(n,c)}$: SMB preprocessing / continuous approximaiton error
        - $\color{green}{\epsilon_{i0}^{(n,q)}}$: unitary-synthesis / MPS-preparation error
        - $\epsilon_{\mathrm{shear}}$: non-orthogonal shear error
        - $\color{blue}{\epsilon_{\mathrm{PK}}}$: phase-kickback error
        - $\color{blue}{\epsilon_{\mathrm{trim}}}$: trimming error

- **Propagation:** 
$\epsilon_{\mathrm{prop}} 
\leq 
\color{blue}{\epsilon_H} t + 
\color{blue}{\epsilon_{\mathrm{QSP}}}$
    - **Hamiltonian block encoding:** 
    $\color{blue}{\epsilon_H} = \epsilon_T + \epsilon_V + \epsilon_\theta$
        - $\epsilon_T$: kinetic-energy block-encoding error
        - $\epsilon_V$: potential-energy block-encoding error
        - $\epsilon_\theta$: $T/V$ weighting-state preparation error
    - **QSP implementation:** 
    $\color{blue}{\epsilon_{\mathrm{QSP}}} = \epsilon_{\tilde d} + 
    (\tilde d+1)( \color{green}{\epsilon_\phi} + 
    \color{green}{\epsilon_\gamma} + 
    \color{green}{\epsilon_{\mathrm{rot}}})$
        - $\epsilon_{\tilde d}$: Jacobi--Anger truncation error
        - $\color{green}{\epsilon_\phi}$: QSP phase precision
        - $\color{green}{\epsilon_\gamma}$: QSP phase precision
        - $\color{green}{\epsilon_{\mathrm{rot}}}$: rotation-synthesis error

- **Basis transformation:** $\color{blue}{\epsilon_B}$

- **Measurement:** 
$\epsilon_{\mathrm{meas}} = 
\epsilon_{\mathrm{QAE}} + \color{blue}{\epsilon_{\mathrm{QAE-QFT}}}$
    - $\epsilon_{\mathrm{QAE}}$: finite-precision quantum amplitude-estimation error
    - $\color{blue}{\epsilon_{\mathrm{QAE-QFT}}}$: QFT error during QAE
