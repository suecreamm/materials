# Electron–Phonon Scattering in Pb

**Mode-resolved coupling, linewidths, self-energy, and quasiparticle lifetimes**

## Question

**How does electron–phonon scattering vary across phonon modes and electronic states?**

A linewidth can contain contributions from several scattering channels,

```math
\Gamma_{\mathrm{tot}}
\approx
\Gamma_{\mathrm{e-ph}}
+
\Gamma_{\mathrm{ph-ph}}
+
\Gamma_{\mathrm{defects}}
+\cdots .
```

This project focuses on the electron–phonon contribution in fcc Pb and examines how it changes with phonon mode, momentum, and electronic state.

---

## Key Result: Mode-Resolved Phonon Scattering

The phonon linewidth $\gamma_{\mathbf q\nu}$ and mode-resolved coupling $\lambda_{\mathbf q\nu}$ show how strongly the electron–phonon interaction varies across the phonon spectrum.

- $\gamma_{\mathbf q\nu}$ measures the electron–phonon contribution to phonon broadening.
- $\lambda_{\mathbf q\nu}$ gives the dimensionless coupling strength of an individual phonon mode.

The figure below maps both quantities onto the phonon branches.

The blue envelope is plotted as

```math
\omega_{\mathbf q\nu}
\pm
6\gamma_{\mathbf q\nu},
```

so a broader blue region indicates a larger phonon linewidth.

The red envelope is plotted as

```math
\omega_{\mathbf q\nu}
\pm
2\lambda_{\mathbf q\nu},
```

to show the relative mode-resolved coupling strength. Because $\lambda_{\mathbf q\nu}$ is dimensionless, the red width is only a visual guide and is not a physical spectral linewidth.

<p align="center">
  <a href="4epw/linewidth_lambda_T0000.075K.png">
    <img src="4epw/linewidth_lambda_T0000.075K.png" alt="Phonon linewidth and electron-phonon coupling" width="650">
  </a>
</p>

The strong variation along the phonon branches shows that electron–phonon scattering depends strongly on both mode and momentum. A single averaged coupling constant therefore cannot capture the full structure of the interaction.

---

## Phonon Dispersion

The phonon dispersion provides the lattice-dynamical basis for the mode-resolved analysis.

The selected high-symmetry path is

```math
\Gamma \rightarrow X \rightarrow W \rightarrow L \rightarrow \Gamma \rightarrow K .
```

<p align="center">
  <a href="phonon_dispersion.png">
    <img src="phonon_dispersion.png" alt="Phonon dispersion" width="650">
  </a>
</p>

---

## Eliashberg Spectral Function

The Eliashberg spectral function $\alpha^2F(\omega)$ shows how phonons at different frequencies contribute to the electron–phonon interaction.

The cumulative coupling is

```math
\lambda(\omega)
=
2\int_0^\omega
\frac{\alpha^2F(\omega')}{\omega'}
\,d\omega' ,
```

and the total coupling constant is

```math
\lambda
=
2\int_0^\infty
\frac{\alpha^2F(\omega)}{\omega}
\,d\omega .
```

For the present calculation,

```math
\lambda \approx 0.686 .
```

While the total $\lambda$ gives an overall measure of the interaction strength, $\alpha^2F(\omega)$ shows which phonon-frequency regions contribute most strongly.

<p align="center">
  <a href="4epw/a2f.png">
    <img src="4epw/a2f.png" alt="Eliashberg spectral function" width="650">
  </a>
</p>

---

## Electron Self-Energy, Quasiparticle Linewidth, and Lifetime

Electron–phonon scattering also modifies the electronic states through the electron self-energy,

```math
\Sigma_{n\mathbf{k}}(\omega)
=
\mathrm{Re}\,\Sigma_{n\mathbf{k}}(\omega)
+
i\,\mathrm{Im}\,\Sigma_{n\mathbf{k}}(\omega) .
```

The real part shifts the quasiparticle energy, while the imaginary part gives a finite spectral linewidth.

In the present post-processing (`calc/4epw/99elself.py`), the plotted on-shell linewidth is defined as

```math
\Gamma^{\mathrm{plot}}_{n\mathbf{k}}
=
2\left|
\mathrm{Im}\,\Sigma_{n\mathbf{k}}
\right| ,
```

with the corresponding lifetime estimate

```math
\tau_{n\mathbf{k}}
\approx
\frac{\hbar}
{\Gamma^{\mathrm{plot}}_{n\mathbf{k}}} .
```

This estimate neglects the quasiparticle renormalization factor $Z$.

<p align="center">
  <a href="4epw/elself_bands.png">
    <img src="4epw/elself_bands.png" alt="Electron self-energy and quasiparticle linewidth" width="650">
  </a>
</p>

For states within $50\,\mathrm{meV}$ of the Fermi level along the selected high-symmetry path:

- Mean linewidth: $16.198\,\mathrm{meV}$
- Maximum linewidth: $26.862\,\mathrm{meV}$
- Mean quasiparticle lifetime: $48.19\,\mathrm{fs}$

These values give the characteristic electron–phonon scattering timescale of the sampled low-energy states and show how it varies across bands and crystal momentum.

---

## Workflow

<p align="center">
  <a href="EPW_workflow.png">
    <img src="EPW_workflow.png" alt="EPW workflow" width="850">
  </a>
</p>

The main calculation stages are:

| Stage | Purpose |
|---|---|
| `1scf` | Ground-state electronic structure |
| `2ph` | DFPT phonons and perturbation potentials |
| `3nscf` | Dense electronic states for Wannier interpolation |
| `4.1epw` | Wannierization and coarse electron–phonon matrix elements |
| `4.2epw` | Interpolation check |
| `4.3epw` | Mode- and q-resolved linewidths and $\lambda_{\mathbf q\nu}$ |
| `4.4epw` | $\alpha^2F(\omega)$ and integrated coupling |
| `4.5epw` | Electron self-energy, linewidths, and quasiparticle lifetimes |

The `1scf`, `2ph`, and `3nscf` calculations were performed with Quantum ESPRESSO 7.5. Electron–phonon quantities were then evaluated and interpolated with EPW 6.0.

**First-principles analysis of electron–phonon scattering in fcc Pb using Quantum ESPRESSO 7.5 and EPW 6.0.**

---

## References

- [EPW School / Tutorial 01 (FCC Lead)](https://docs.epw-code.org/tutorials/tutorial_01/index.html)
- [EPW Input Variables](https://docs.epw-code.org/doc/Inputs.html)
- [EPW: Electron–phonon coupling using Wannier functions](https://arxiv.org/abs/1604.03525)
- G. D. Mahan, *Many-Particle Physics*, 3rd ed.
