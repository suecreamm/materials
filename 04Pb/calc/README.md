## Pb EPW Example

### Workflow
[![EPW workflow](EPW_workflow.png)](EPW_workflow.png)

### Phonon Dispersion
<p align="center">
  <a href="phonon_dispersion.png">
    <img src="phonon_dispersion.png" alt="Phonon dispersion" width="600">
  </a>
</p>

### Phonon Linewidth & Electron–Phonon Coupling
Mode-resolved phonon linewidth and electron–phonon coupling.

[![Phonon linewidth and electron-phonon coupling](4epw/linewidth_lambda_T0000.075K.png)](4epw/linewidth_lambda_T0000.075K.png)

### Eliashberg Spectral Function
The Eliashberg spectral function, α²F(ω), was calculated using EPW, giving an electron–phonon coupling constant of λ = 0.686.

<p align="center">
  <a href="4epw/a2f.png">
    <img src="4epw/a2f.png" alt="Eliashberg spectral function" width="600">
  </a>
</p>

### Electron Self-Energy
`99elself.py` post-processes `linewidth.elself.0.075K` and analyzes the electron self-energy along the selected k-path.

For states within 50 meV of the Fermi level:
- Mean linewidth: 16.198 meV
- Maximum linewidth: 26.862 meV
- Mean lifetime: 48.19 fs

<p align="center">
  <a href="4epw/elself_bands.png">
    <img src="4epw/elself_bands.png" alt="Electron self-energy" width="600">
  </a>
</p>

### References

- [EPW School / Tutorial 01](https://docs.epw-code.org/tutorials/tutorial_01/index.html)
- [EPW FCC Lead Tutorial](https://docs.epw-code.org/tutorials/FCC-lead.html)
- [EPW Input Variables](https://docs.epw-code.org/doc/Inputs.html)
- [EPW: Electron–phonon coupling using Wannier functions](https://arxiv.org/abs/1604.03525)
