# Pb EPW Example

## Overview

This repository presents a compact first-principles EPW workflow for fcc Pb, including phonon dispersion, electron–phonon coupling, Eliashberg spectral functions, and electron self-energy analysis. The example is intended to demonstrate my practical experience with first-principles calculations of lattice dynamics and electron–phonon interactions, with particular emphasis on linewidths, coupling strengths, and quasiparticle lifetimes.

## Table of Contents

- [Workflow](#workflow)
- [Phonon Dispersion](#phonon-dispersion)
- [Phonon Linewidth & Electron–Phonon Coupling](#phonon-linewidth--electronphonon-coupling)
- [Eliashberg Spectral Function](#eliashberg-spectral-function)
- [Electron Self-Energy](#electron-self-energy)
- [References](#references)

## Workflow

[![EPW workflow](EPW_workflow.png)](EPW_workflow.png)

## Phonon Dispersion

<p align="center">
  <a href="phonon_dispersion.png">
    <img src="phonon_dispersion.png" alt="Phonon dispersion" width="600">
  </a>
</p>

## Phonon Linewidth & Electron–Phonon Coupling

The phonon linewidth reflects the decay rate of a phonon mode due to electron–phonon scattering, while the mode-resolved coupling parameter, $\lambda_{\mathbf{q}\nu}$, quantifies the contribution of each phonon mode to the total electron–phonon interaction. Together, they provide how strongly different phonon branches couple to the electronic states across reciprocal space.

The figure below shows the momentum- and mode-resolved variation of these quantities across the Brillouin zone.

[![Phonon linewidth and electron-phonon coupling](4epw/linewidth_lambda_T0000.075K.png)](4epw/linewidth_lambda_T0000.075K.png)

## Eliashberg Spectral Function

The Eliashberg spectral function, $\alpha^2F(\omega)$, provides a frequency-resolved measure of the electron–phonon interaction by weighting the phonon density of states with the corresponding coupling matrix elements. Its integral determines the dimensionless electron–phonon coupling constant,

$$
\lambda = 2 \int_0^\infty \frac{\alpha^2F(\omega)}{\omega}\, d\omega
$$

For fcc Pb, the present calculation gives

$$
\lambda = 0.686
$$

The spectral shape of $\alpha^2F(\omega)$ shows which phonon modes contribute most strongly to the total electron–phonon coupling, while $\lambda$ provides a compact measure of the overall interaction strength.

<p align="center">
  <a href="4epw/a2f.png">
    <img src="4epw/a2f.png" alt="Eliashberg spectral function" width="600">
  </a>
</p>

## Electron Self-Energy

The electron self-energy provides a momentum- and band-resolved description of electron–phonon scattering. The imaginary part of the self-energy determines the quasiparticle linewidth,

$$
\Gamma_{n\mathbf{k}} = 2\mathrm{Im}\Sigma_{n\mathbf{k}}
$$

which is related to the quasiparticle scattering rate and, within the quasiparticle picture, as an inverse lifetime or scattering rate.

$$
\tau_{n\mathbf{k}} \approx \frac{\hbar}{\Gamma_{n\mathbf{k}}}
$$

The figure below shows the calculated linewidth along the electronic band structure and highlights regions in reciprocal space where electron–phonon scattering is strong. States close to the Fermi level are particularly relevant for low-energy electronic scattering processes.

<p align="center">
  <a href="4epw/elself_bands.png">
    <img src="4epw/elself_bands.png" alt="Electron self-energy" width="600">
  </a>
</p>

For states within 50 meV of the Fermi level:

- Mean linewidth: 16.198 meV
- Maximum linewidth: 26.862 meV
- Mean lifetime: 48.19 fs

These linewidths provide a simple estimate of the characteristic electron–phonon scattering timescale. The variation of $\Gamma_{n\mathbf{k}}$ across the selected states also illustrates that the scattering strength depends strongly on both band index and crystal momentum.

`calc/4epw/99elself.py` post-processes `linewidth.elself.0.075K` and analyzes the electron self-energy along the selected k-path.

## References

- [EPW School / Tutorial 01 (FCC Lead)](https://docs.epw-code.org/tutorials/tutorial_01/index.html)
- [EPW Input Variables](https://docs.epw-code.org/doc/Inputs.html)
- [EPW: Electron–phonon coupling using Wannier functions](https://arxiv.org/abs/1604.03525)
