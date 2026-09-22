# Pb EPW Example

## Overview

This repository presents a compact first-principles EPW workflow for fcc Pb, including phonon dispersion, phonon linewidths, mode-resolved electron–phonon coupling, the Eliashberg spectral function, and electron self-energy analysis. The example demonstrates practical calculations of lattice dynamics and electron–phonon interactions, with emphasis on linewidths, coupling strengths, and quasiparticle lifetimes.

## Table of Contents

* [Workflow](#workflow)
* [Phonon Dispersion](#phonon-dispersion)
* [Phonon Linewidth & Electron–Phonon Coupling](#phonon-linewidth--electronphonon-coupling)
* [Eliashberg Spectral Function](#eliashberg-spectral-function)
* [Electron Self-Energy, Quasiparticle Linewidth, and Lifetime](#electron-self-energy-quasiparticle-linewidth-and-lifetime)
* [References](#references)

## Workflow

[![EPW workflow](EPW_workflow.png)](EPW_workflow.png)

## Phonon Dispersion

The phonon dispersion provides the phonon spectrum used for the subsequent electron–phonon analysis.

The selected high-symmetry path is

```math
\Gamma \rightarrow X \rightarrow W \rightarrow L \rightarrow \Gamma \rightarrow K .
```

<p align="center">
  <a href="phonon_dispersion.png">
    <img src="phonon_dispersion.png" alt="Phonon dispersion" width="650">
  </a>
</p>

## Phonon Linewidth & Electron–Phonon Coupling

Electron–phonon scattering contributes to the linewidth of a phonon mode. For wave vector $`\mathbf{q}`$ and branch index $`\nu`$, the corresponding linewidth is $`\gamma_{\mathbf{q}\nu}`$. A larger linewidth indicates a stronger electron–phonon contribution to phonon broadening and decay.

The mode-resolved coupling strength $`\lambda_{\mathbf{q}\nu}`$ describes the contribution of an individual phonon mode to the electron–phonon interaction. The two quantities describe different aspects of the same mode-resolved electron–phonon interaction: $`\gamma_{\mathbf{q}\nu}`$ characterizes phonon damping, while $`\lambda_{\mathbf{q}\nu}`$ gives the corresponding dimensionless coupling strength.

In the figure below, the phonon branches are surrounded by shaded envelopes used to visualize the mode-resolved quantities. The blue envelope is plotted as $`\omega \pm 6\gamma_{\mathbf{q}\nu}`$, so a broader blue region indicates a larger phonon linewidth. The red envelope is plotted as $`\omega \pm 2\lambda_{\mathbf{q}\nu}`$ as a visual encoding of the relative mode-resolved coupling strength. Because $`\lambda_{\mathbf{q}\nu}`$ is dimensionless, the red width should not be interpreted as a physical spectral linewidth.

<p align="center">
  <a href="4epw/linewidth_lambda_T0000.075K.png">
    <img src="4epw/linewidth_lambda_T0000.075K.png" alt="Phonon linewidth and electron-phonon coupling" width="650">
  </a>
</p>

## Eliashberg Spectral Function

The Eliashberg spectral function $`\alpha^2F(\omega)`$ shows how phonons at different frequencies contribute to the electron–phonon interaction.

The cumulative coupling strength is

```math
\lambda(\omega) = 2\int_0^\omega \frac{\alpha^2F(\omega')}{\omega'}\,d\omega' .
```

The total electron–phonon coupling constant is

```math
\lambda = 2\int_0^\infty \frac{\alpha^2F(\omega)}{\omega}\,d\omega .
```

For the present calculation,

```math
\lambda \approx 0.686 .
```

The frequency dependence of $`\alpha^2F(\omega)`$ identifies which parts of the phonon spectrum contribute most strongly to the total coupling.

<p align="center">
  <a href="4epw/a2f.png">
    <img src="4epw/a2f.png" alt="Eliashberg spectral function" width="650">
  </a>
</p>

## Electron Self-Energy, Quasiparticle Linewidth, and Lifetime

Electron–phonon scattering modifies the electronic states through the electron self-energy,

```math
\Sigma_{n\mathbf{k}}(\omega)
=
\mathrm{Re}\,\Sigma_{n\mathbf{k}}(\omega)
+
i\,\mathrm{Im}\,\Sigma_{n\mathbf{k}}(\omega) .
```

The real part shifts the quasiparticle energy, while the imaginary part produces spectral broadening associated with a finite quasiparticle lifetime.

In the present post-processing (`calc/4epw/99elself.py`), an on-shell electron-linewidth estimate is defined as

```math
\Gamma^{\mathrm{plot}}_{n\mathbf{k}}
=
2\left|\mathrm{Im}\,\Sigma_{n\mathbf{k}}\right| ,
```

with the corresponding lifetime estimate

```math
\tau_{n\mathbf{k}}
\approx
\frac{\hbar}{\Gamma^{\mathrm{plot}}_{n\mathbf{k}}} .
```

This estimate neglects the quasiparticle renormalization factor $`Z`$. A larger linewidth therefore corresponds to a shorter quasiparticle lifetime.

In the figure below, both marker color and size represent the calculated electron linewidth, allowing regions of stronger electron–phonon scattering to be identified along the electronic bands.

<p align="center">
  <a href="4epw/elself_bands.png">
    <img src="4epw/elself_bands.png" alt="Electron self-energy and quasiparticle linewidth" width="650">
  </a>
</p>

For the sampled electronic states along the selected high-symmetry path within $`50\,\mathrm{meV}`$ of the Fermi level $`E_F`$:

* Mean linewidth: $`16.198\,\mathrm{meV}`$
* Maximum linewidth: $`26.862\,\mathrm{meV}`$
* Mean quasiparticle lifetime: $`48.19\,\mathrm{fs}`$

These values provide an estimate of the characteristic electron–phonon scattering timescale for the selected low-energy states and show its dependence on band index and crystal momentum.

`calc/4epw/99elself.py` post-processes `linewidth.elself.0.075K` and analyzes the imaginary part of the electron self-energy along the selected k-path.

## References

* [EPW School / Tutorial 01 (FCC Lead)](https://docs.epw-code.org/tutorials/tutorial_01/index.html)
* [EPW Input Variables](https://docs.epw-code.org/doc/Inputs.html)
* [EPW: Electron–phonon coupling using Wannier functions](https://arxiv.org/abs/1604.03525)
* G. D. Mahan, *Many-Particle Physics*, 3rd ed.
