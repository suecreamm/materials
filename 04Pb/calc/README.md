# Pb EPW Example

## Overview

This repository presents a compact first-principles EPW workflow for fcc Pb, connecting electron–phonon matrix elements to phonon linewidths, electron–phonon coupling strengths, electronic self-energies, and characteristic lifetimes. Pb is a conventional electron–phonon-mediated superconductor and therefore provides a useful benchmark system for analyzing these quantities within a consistent first-principles framework.

Rather than only demonstrating the computational workflow, this example illustrates how a microscopic scattering mechanism can be connected to energy dissipation, spectral broadening, and finite excitation lifetimes. Here, electron–phonon scattering is treated explicitly as one contribution to the damping, while other processes—such as anharmonic phonon–phonon scattering and defect-induced scattering—can provide additional decay channels that determine the observed lifetime.

When different scattering channels can be treated approximately independently, their contributions to the linewidth may be expressed schematically as

```math
\Gamma_{\mathbf q\nu}^{\mathrm{tot}}
\approx
\Gamma_{\mathbf q\nu}^{e\text{-}ph}
+
\Gamma_{\mathbf q\nu}^{ph\text{-}ph}
+
\Gamma_{\mathbf q\nu}^{\mathrm{defect}}
+\cdots .
```

More generally, individual scattering channels can be analyzed by including or omitting them separately to quantify their contributions to the resulting linewidth and lifetime. The present EPW calculation provides the electron–phonon contribution to this decomposition.

The same electron–phonon interaction can be viewed from two complementary perspectives:

**Phonon perspective**

```math
e\text{-}ph
\;\longrightarrow\;
\Gamma_{\mathbf q\nu}^{\mathrm{ph}}
\;\longrightarrow\;
\tau_{\mathbf q\nu}^{\mathrm{ph}}
```

Electron–phonon coupling gives a finite phonon linewidth, corresponding to a finite lifetime of the phonon excitation.

**Electronic perspective**

```math
e\text{-}ph
\;\longrightarrow\;
\mathrm{Im}\,\Sigma_{n\mathbf k}
\;\longrightarrow\;
\text{spectral broadening}
\;\longrightarrow\;
\tau_{n\mathbf k}^{\mathrm{el}}
```

From the electronic side, the same interaction contributes to the imaginary part of the electron self-energy, producing quasiparticle broadening and a finite electronic lifetime.

---

## Table of Contents

* [Workflow](#workflow)
* [Phonon Dispersion](#phonon-dispersion)
* [Phonon Linewidth: Electron–Phonon Contribution to Phonon Decay](#phonon-linewidth)
* [Eliashberg Spectral Function](#eliashberg-spectral-function)
* [Electronic Self-Energy and Quasiparticle Lifetime](#electronic-self-energy)
* [References](#references)

---

<a id="workflow"></a>

## Workflow

[![EPW workflow](EPW_workflow.png)](EPW_workflow.png)

*The workflow connects ground-state DFT and lattice-dynamical calculations to Wannier-interpolated electron–phonon quantities, allowing scattering properties to be resolved on dense momentum grids.*

---

<a id="phonon-dispersion"></a>

## Phonon Dispersion

<p align="center">
  <a href="phonon_dispersion.png">
    <img src="phonon_dispersion.png" alt="Phonon dispersion" width="600">
  </a>
</p>

Phonon dispersion shows the vibrational frequencies of a crystal across the Brillouin zone. It provides the phonon modes ((\mathbf q,\nu)) used in the subsequent electron–phonon coupling and linewidth calculations.

---

<a id="phonon-linewidth"></a>

## Phonon Linewidth: Electron–Phonon Contribution to Phonon Decay

The figure below shows the momentum- and mode-resolved variation of the phonon linewidth and electron–phonon coupling.

![Phonon linewidth and electron-phonon coupling](4epw/linewidth_lambda_T0000.075K.png)

Electron–phonon scattering gives each phonon mode a finite linewidth, (\Gamma_{\mathbf q\nu}), which is related to its decay rate. A larger linewidth means stronger electron–phonon damping and therefore a shorter phonon lifetime.

The width around each phonon branch visualizes the mode-dependent linewidth or coupling strength. Larger values identify the modes where electron–phonon interaction has a stronger effect.

---

<a id="eliashberg-spectral-function"></a>

## Eliashberg Spectral Function

The Eliashberg spectral function, $\alpha^2F(\omega)$, provides a frequency-resolved measure of the electron–phonon interaction by weighting the phonon spectrum with the corresponding electron–phonon coupling.

Its integral determines the dimensionless electron–phonon coupling constant,

```math
\lambda
=
2 \int_0^\infty
\frac{\alpha^2F(\omega)}{\omega}
\,d\omega .
```

For fcc Pb, the present calculation gives

```math
\lambda = 0.686 .
```

The spectral shape of $\alpha^2F(\omega)$ shows which phonon-frequency regions contribute most strongly to the total electron–phonon coupling, while $\lambda$ provides a compact measure of the overall interaction strength.

<p align="center">
  <a href="4epw/a2f.png">
    <img src="4epw/a2f.png" alt="Eliashberg spectral function" width="600">
  </a>
</p>

*The peaks in $\alpha^2F(\omega)$ reveal the phonon-energy ranges that contribute most strongly to the integrated electron–phonon coupling.*

---

<a id="electronic-self-energy"></a>

## Electronic Self-Energy and Quasiparticle Lifetime

The figure below shows the calculated electronic linewidth along the band structure and highlights regions in reciprocal space where electron–phonon scattering is strong. States close to the Fermi level are particularly relevant for low-energy electronic scattering processes.

<p align="center">
  <a href="4epw/elself_bands.png">
    <img src="4epw/elself_bands.png" alt="Electron self-energy and linewidth" width="600">
  </a>
</p>


When an electron interacts with phonons, its energy and lifetime are modified through the complex self-energy,

\Sigma_{n\mathbf k}
=
\mathrm{Re}\,\Sigma_{n\mathbf k}
+
i\,\mathrm{Im}\,\Sigma_{n\mathbf k}.

The real part shifts the quasiparticle energy, while the imaginary part produces spectral broadening. A larger (|\mathrm{Im}\Sigma|) therefore corresponds to a larger linewidth and a shorter quasiparticle lifetime.

\Gamma_{n\mathbf k}
\simeq
2|\mathrm{Im}\Sigma_{n\mathbf k}|,
\qquad
\tau_{n\mathbf k}
\simeq
\frac{\hbar}{\Gamma_{n\mathbf k}} .

<p align="center"> <a href="4epw/elself_bands.png"> <img src="4epw/elself_bands.png" alt="Electron self-energy and linewidth" width="600"> </a> </p>

Both marker color and size represent the electron linewidth. Yellow, larger markers indicate larger (\Gamma_{n\mathbf k}): stronger electron–phonon scattering, broader quasiparticle states, and shorter lifetimes.

For states within 50 meV of the Fermi level:

Mean linewidth: 16.198 meV
Maximum linewidth: 26.862 meV
Mean lifetime: 48.19 fs

`calc/4epw/99elself.py` post-processes `linewidth.elself.0.075K` and analyzes the electron self-energy along the selected k-path.

---

<a id="references"></a>

## References

* [EPW School / Tutorial 01 (FCC Lead)](https://docs.epw-code.org/tutorials/tutorial_01/index.html)
* [EPW Input Variables](https://docs.epw-code.org/doc/Inputs.html)
* [EPW: Electron–phonon coupling using Wannier functions](https://arxiv.org/abs/1604.03525)
* G. D. Mahan, *Many-Particle Physics*, 3rd ed., Ch. 8.
