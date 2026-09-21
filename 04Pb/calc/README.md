# Pb EPW Example

## Overview

This repository presents a compact first-principles EPW workflow for fcc Pb, including phonon dispersion, phonon linewidths, mode-resolved electron–phonon coupling, Eliashberg spectral functions, and electron self-energy analysis. The example is intended to demonstrate my practical experience with first-principles calculations of lattice dynamics and electron–phonon interactions, with particular emphasis on linewidths, coupling strengths, and quasiparticle lifetimes.

## Table of Contents

* [Workflow](#workflow)
* [Phonon Dispersion](#phonon-dispersion)
* [Phonon Linewidth & Electron–Phonon Coupling](#phonon-linewidth--electronphonon-coupling)
* [Eliashberg Spectral Function](#eliashberg-spectral-function)
* [Electron Self-Energy](#electron-self-energy)
* [References](#references)

## Workflow

This example follows a compact EPW workflow for fcc Pb:

1. `1scf`: ground-state DFT calculation
2. `2phonon`: DFPT phonon calculation
3. `3nscf`: NSCF calculation on a dense k-mesh
4. `4.1epw`: Wannierization and coarse-grid electron–phonon matrix elements
5. `4.2epw`: interpolation check
6. `4.3epw`: phonon linewidth and mode-resolved coupling along the selected q-path
7. `4.4epw`: Eliashberg spectral function \(\alpha^2F(\omega)\) and total coupling \(\lambda\) on a uniform grid
8. `4.5epw`: electron self-energy, linewidth, and lifetime analysis along the selected k-path

[![EPW workflow](EPW_workflow.png)](EPW_workflow.png)

## Phonon Dispersion

The phonon dispersion shows the allowed lattice-vibration frequencies across the Brillouin zone. It provides the basic phonon spectrum used in the subsequent electron–phonon calculations and allows the behavior of different phonon branches to be compared along the selected high-symmetry path.

The calculated path is \(\Gamma \rightarrow X \rightarrow W \rightarrow L \rightarrow \Gamma \rightarrow K.\)

<p align="center">
  <a href="phonon_dispersion.png">
    <img src="phonon_dispersion.png" alt="Phonon dispersion" width="650">
  </a>
</p>

## Phonon Linewidth & Electron–Phonon Coupling

Electron–phonon scattering contributes to the linewidth of a phonon mode. For a phonon with wave vector \(\mathbf{q}\) and branch index \(\nu\), the corresponding linewidth is written as \(\gamma_{\mathbf{q}\nu}.\) A larger \(\gamma_{\mathbf{q}\nu}\) means that electron–phonon scattering contributes more strongly to the broadening and decay of that phonon mode.

The mode-resolved electron–phonon coupling strength is denoted by \(\lambda_{\mathbf{q}\nu}.\) While \(\gamma_{\mathbf{q}\nu}\) describes the electron–phonon contribution to phonon broadening, \(\lambda_{\mathbf{q}\nu}\) describes how strongly an individual phonon mode contributes to the total electron–phonon interaction. These two quantities are therefore related but not identical.

In the figure below, the marker size represents \(\gamma_{\mathbf{q}\nu}\), so larger markers indicate phonon modes with stronger electron–phonon damping.

<p align="center">
  <a href="4epw/linewidth_lambda_T0000.075K.png">
    <img src="4epw/linewidth_lambda_T0000.075K.png" alt="Phonon linewidth and electron-phonon coupling" width="650">
  </a>
</p>

## Eliashberg Spectral Function

The Eliashberg spectral function \(\alpha^2F(\omega)\) describes how phonons at different frequencies contribute to the electron–phonon interaction. The cumulative coupling strength is given by \(\lambda(\omega) = 2 \int_0^\omega \frac{\alpha^2F(\omega')}{\omega'}\,d\omega'.\)

The total electron–phonon coupling constant is obtained from \(\lambda = 2 \int_0^\infty \frac{\alpha^2F(\omega)}{\omega}\,d\omega.\)

For the present calculation, \(\lambda \approx 0.686.\)

The frequency dependence of \(\alpha^2F(\omega)\) shows which parts of the phonon spectrum contribute most strongly to the total electron–phonon coupling.

<p align="center">
  <a href="4epw/a2f.png">
    <img src="4epw/a2f.png" alt="Eliashberg spectral function" width="650">
  </a>
</p>

## Electron Self-Energy

Electron–phonon scattering modifies the electronic states through the electron self-energy, \(\Sigma_{n\mathbf{k}}(\omega) = \mathrm{Re}\,\Sigma_{n\mathbf{k}}(\omega) + i\,\mathrm{Im}\,\Sigma_{n\mathbf{k}}(\omega).\)

The real part, \(\mathrm{Re}\,\Sigma,\) shifts the quasiparticle energy, while the imaginary part, \(\mathrm{Im}\,\Sigma,\) gives the broadening caused by scattering.

In the present post-processing (`calc/4epw/99elself.py`), the plotted electron linewidth is defined as \(\Gamma^{\mathrm{plot}}_{n\mathbf{k}} = 2\left|\mathrm{Im}\,\Sigma_{n\mathbf{k}}\right|,\) and the corresponding lifetime estimate is \(\tau_{n\mathbf{k}} \approx \frac{\hbar}{\Gamma^{\mathrm{plot}}_{n\mathbf{k}}}.\)

In the figure below, the electronic bands show where the quasiparticle states occur in momentum and energy, while both the marker color and size represent the calculated electron linewidth. This makes it possible to identify the momentum and band regions where electron–phonon scattering is strongest.

<p align="center">
  <a href="4epw/elself_bands.png">
    <img src="4epw/elself_bands.png" alt="Electron self-energy" width="650">
  </a>
</p>

For the sampled electronic states along the selected high-symmetry path within \(50\,\mathrm{meV}\) of the Fermi level \(E_F\):

* Mean linewidth: \(16.198\,\mathrm{meV}\)
* Maximum linewidth: \(26.862\,\mathrm{meV}\)
* Mean quasiparticle lifetime: \(48.19\,\mathrm{fs}\)

These values provide a compact estimate of the characteristic electron–phonon scattering timescale for the selected low-energy states. The variation of \(\Gamma^{\mathrm{plot}}_{n\mathbf{k}}\) along the band structure also shows that the scattering strength depends strongly on both band index and crystal momentum.

`calc/4epw/99elself.py` post-processes `linewidth.elself.0.075K` and analyzes the electron self-energy along the selected k-path.

## References

* [EPW School / Tutorial 01 (FCC Lead)](https://docs.epw-code.org/tutorials/tutorial_01/index.html)
* [EPW Input Variables](https://docs.epw-code.org/doc/Inputs.html)
* [EPW: Electron–phonon coupling using Wannier functions](https://arxiv.org/abs/1604.03525)
