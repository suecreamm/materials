## Phonon Dispersion

The phonon dispersion shows the allowed lattice-vibration frequencies across the Brillouin zone. It provides the basic phonon spectrum used in the following electron–phonon calculations and allows the behavior of different phonon branches to be compared along the selected high-symmetry path.

The calculated path is

$$
\Gamma \rightarrow X \rightarrow W \rightarrow L \rightarrow \Gamma \rightarrow K .
$$

---

## Phonon Linewidth & Contribution to Phonon Decay

Electron–phonon scattering gives a finite linewidth to each phonon mode. For a phonon with wave vector $\mathbf{q}$ and branch index $\nu$, the linewidth is written as $\gamma_{\mathbf{q}\nu}$.

A larger $\gamma_{\mathbf{q}\nu}$ means that electron–phonon scattering contributes more strongly to the decay of that phonon mode.

The figure also shows the mode-resolved electron–phonon coupling strength $\lambda_{\mathbf{q}\nu}$. The **marker size represents $\lambda_{\mathbf{q}\nu}$**: larger markers correspond to phonon modes with stronger electron–phonon coupling.

Therefore, the plot provides two complementary pieces of information:

* $\gamma_{\mathbf{q}\nu}$ describes the electron–phonon contribution to phonon broadening.
* $\lambda_{\mathbf{q}\nu}$ describes the coupling strength of the individual phonon mode.

---

## Eliashberg Spectral Function

The Eliashberg spectral function $\alpha^2F(\omega)$ describes how phonons at different frequencies contribute to the electron–phonon interaction.

The cumulative electron–phonon coupling constant is

$$
\lambda(\omega)
=
2\int_0^\omega
\frac{\alpha^2F(\omega')}{\omega'}
\,d\omega' .
$$

The total coupling constant is obtained from

$$
\lambda
=
2\int_0^\infty
\frac{\alpha^2F(\omega)}{\omega}
\,d\omega .
$$

For this calculation,

$$
\lambda \approx 0.686 .
$$

The frequency dependence of $\alpha^2F(\omega)$ shows which parts of the phonon spectrum contribute most strongly to the total electron–phonon coupling.

---

## Electron Self-Energy and Quasiparticle Lifetime

Electron–phonon scattering modifies the electronic states through the electron self-energy

$$
\Sigma_{n\mathbf{k}}(\omega)
=
\operatorname{Re}\Sigma_{n\mathbf{k}}(\omega)
+
i\,\operatorname{Im}\Sigma_{n\mathbf{k}}(\omega) .
$$

The real part, $\operatorname{Re}\Sigma$, changes the quasiparticle energy, while the imaginary part, $\operatorname{Im}\Sigma$, gives the broadening caused by scattering.

The electron linewidth is related to the imaginary part of the self-energy by

$$
\Gamma_{n\mathbf{k}}
=
2\left|
\operatorname{Im}\Sigma_{n\mathbf{k}}
\right| .
$$

A finite linewidth means that the quasiparticle has a finite lifetime. The corresponding lifetime can be estimated as

$$
\tau_{n\mathbf{k}}
=
\frac{\hbar}
{\Gamma_{n\mathbf{k}}} .
$$

Therefore, stronger electron–phonon scattering produces a larger linewidth and a shorter quasiparticle lifetime.

In the figure, the electronic bands show where the quasiparticle states occur in momentum and energy, while the **color represents the calculated electron linewidth**. This makes it possible to identify the momentum and band regions where electron–phonon scattering is strongest.

For electronic states within $50~\mathrm{meV}$ of the Fermi level $E_F$:

* mean linewidth: $16.198~\mathrm{meV}$
* maximum linewidth: $26.862~\mathrm{meV}$
* mean quasiparticle lifetime: $48.19~\mathrm{fs}$
