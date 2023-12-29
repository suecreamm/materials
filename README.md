# Materials

## 1 Graphene

- VASP Input Files: POSCAR, INCAR, KPOINTS

### You can get graphene's ...
- Density of States
- Electronic Bandstructures (G-M-K-G)

## 2 Carbon Nanotubes

### Atomic Structures of CNT(VASP POSCAR & CIF)
#### Web CNT Generator
![Carbon Nanotubes Generator Demonstration](https://github.com/suecreamm/cnt_generator/raw/main/readmeImg.png)
<br><br><a href="https://suecreamm.github.io/cnt_generator/" target="_blank" rel="noopener noreferrer">🔮 Go to the Website</a><br>
With just a few clicks, you can generate the atomic structures of carbon nanotubes. No installation required!  

***

## Scripts

### B001_phononJobAuto.sh
Phonopy allows you to get the phonon dispersion.<br><br>
``` $phonopy -d --dim 2 2 2 ```
<br>This gets us a set of supercells.<br>

``` $bash B001_phononJobAuto.sh ``` <br>
Run this script to automatically assign cluster jobs for phonon calculations. My script may be more efficient when Phonopy provides more than 10 structures.
