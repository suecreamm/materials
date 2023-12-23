# Materials

## 1 Graphene

## 2 Carbon Nanotubes

### Carbon Nanotubes' Atomic Structures(VASP POSCAR & CIF)
#### Web CNT Generator
![Carbon Nanotubes Generator Demonstration](https://github.com/suecreamm/cnt_generator/raw/main/readmeImg.png)
<br><br><a href="https://suecreamm.github.io/cnt_generator/" target="_blank" rel="noopener noreferrer">🔮 Go to the Website</a><br>
This generates the atomic structures of carbon nanotubes with only a few clicks! No installation is needed!  

***

## Useful Scripts

### B001_phononJobAuto.sh
With Phonopy, you can get the phonon dispersion.<br><br>
``` $phonopy -d --dim 2 2 2 ```
<br>This gets us a set of supercells.<br>

``` $bash B001_phononJobAuto.sh ``` <br>
Just run this script! This would automatically assign your cluster jobs for the phonon calculations. My script may be efficient when Phonopy gives us more than 10 structures.
