# Pb EPW Example

## Overview

This repository presents a compact first-principles EPW workflow for fcc Pb, from phonon calculations to mode-resolved electron–phonon coupling, Eliashberg spectral functions, and electron self-energy analysis.

## Table of Contents

* [Workflow](#workflow)
* [Phonon Dispersion](#phonon-dispersion)
* [Phonon Linewidth & Electron–Phonon Coupling](#phonon-linewidth--electronphonon-coupling)
* [Eliashberg Spectral Function](#eliashberg-spectral-function)
* [Electron Self-Energy, Quasiparticle Linewidth, and Lifetime](#electron-self-energy-quasiparticle-linewidth-and-lifetime)
* [References](#references)

## Workflow

The calculation proceeds from DFT and DFPT to Wannier interpolation, mode-resolved electron–phonon analysis, $\alpha^2F(\omega)$, and electron self-energy calculations.

[![EPW workflow](EPW_workflow.png)](EPW_workflow.png)

## Phonon Dispersion

The phonon dispersion provides the spectrum used for the subsequent electron–phonon analysis along

\(\Gamma \rightarrow X \rightarrow W \rightarrow L \rightarrow \Gamma \rightarrow K.\)

<p align="center">
  <a href="phonon_dispersion.png">
    <img src="phonon_dispersion.png" alt="Phonon dispersion" width="650">
  </a>
</p>

## Phonon Linewidth & Electron–Phonon Coupling

The phonon linewidth $\gamma_{\mathbf{q}\nu}$ measures the electron–phonon contribution to phonon broadening, while $\lambda_{\mathbf{q}\nu}$ describes the mode-resolved coupling strength.

In the figure below, the marker size represents $\gamma_{\mathbf{q}\nu}$: larger markers indicate stronger electron–phonon damping.

<p align="center">
  <a href="4epw/linewidth_lambda_T0000.075K.png">
    <img src="4epw/linewidth_lambda_T0000.075K.png" alt="Phonon linewidth and electron-phonon coupling" width="650">
  </a>
</p>

## Eliashberg Spectral Function

The Eliashberg spectral function $\alpha^2F(\omega)$ shows how phonons at different frequencies contribute to the electron–phonon interaction.

The total coupling constant is

\(\lambda = 2 \int_0^\infty \frac{\alpha^2F(\omega)}{\omega}\,d\omega.\)

For the present calculation,

\(\lambda \approx 0.686.\)

<p align="center">
  <a href="4epw/a2f.png">
    <img src="4epw/a2f.png" alt="Eliashberg spectral function" width="650">
  </a>
</p>

## Electron Self-Energy, Quasiparticle Linewidth, and Lifetime

Electron–phonon scattering modifies the electronic states through the self-energy

\(\Sigma_{n\mathbf{k}}(\omega) = \mathrm{Re}\,\Sigma_{n\mathbf{k}}(\omega) + i\,\mathrm{Im}\,\Sigma_{n\mathbf{k}}(\omega).\)

The real part shifts the quasiparticle energy, while the imaginary part produces spectral broadening and a finite quasiparticle lifetime.

In the present post-processing, the plotted linewidth is defined as

\(\Gamma^{\mathrm{plot}}_{n\mathbf{k}} = 2\left|\mathrm{Im}\,\Sigma_{n\mathbf{k}}\right|,\)

with

\(\tau_{n\mathbf{k}} \approx \frac{\hbar}{\Gamma^{\mathrm{plot}}_{n\mathbf{k}}}.\)

Both marker color and size in the figure represent the calculated electron linewidth.

<p align="center">
  <a href="4epw/elself_bands.png">
    <img src="4epw/elself_bands.png" alt="Electron self-energy and quasiparticle linewidth" width="650">
  </a>
</p>

For sampled states along the selected high-symmetry path within $50,\mathrm{meV}$ of $E_F$:

* Mean linewidth: $16.198,\mathrm{meV}$
* Maximum linewidth: $26.862,\mathrm{meV}$
* Mean quasiparticle lifetime: $48.19,\mathrm{fs}$

`calc/4epw/99elself.py` post-processes `linewidth.elself.0.075K` and analyzes the electron self-energy along the selected k-path.

## References

* [EPW School / Tutorial 01 (FCC Lead)](https://docs.epw-code.org/tutorials/tutorial_01/index.html)
* [EPW Input Variables](https://docs.epw-code.org/doc/Inputs.html)
* [EPW: Electron–phonon coupling using Wannier functions](https://arxiv.org/abs/1604.03525)
* G. D. Mahan, *Many-Particle Physics*, 3rd ed.
