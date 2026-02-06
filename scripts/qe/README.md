# QE / EPW Utility Scripts ⚙️🧪

Utilities for **Quantum ESPRESSO (QE)** and **EPW** workflows, organized by category.

---

## 🌐 Language / 언어
- 한국어
- English

==============================
한국어 🇰🇷
==============================


### 📚 Table of Contents

#### 🔗 EPW / PH Utilities
- link_ph_outputs.py

#### 📊 Electronic Structure (Bands / DOS)
- pdosPlot.py
- wannier_plot_check.py

#### 📈 Phonon Analysis
- plot_phonon_dispersion.py

#### 🧮 EPW Post-processing
- postprocess_epw.py

#### 🚀 Automation / SLURM
- postprocess.sh

------------------------------------------------------------
🔗 EPW / PH Utilities
------------------------------------------------------------

link_ph_outputs.py 🔗

Purpose  
EPW가 요구하는 phonon 파일 naming 규칙에 맞추기 위해  
PH 계산 결과(dyn, dvscf)에 대해 안전한 symbolic link를 생성합니다.

생성되는 이름:
- PREFIX.dyn_qN
- PREFIX.dvscf_qN
- PREFIX.dvscfN_1

Where to run  
- PREFIX.dyn* 파일들이 존재하는 디렉토리  
  (보통 ph.x 계산이 끝난 작업 디렉토리)

Usage
python link_ph_outputs.py PREFIX
python link_ph_outputs.py PREFIX --dvscf_dir ./tmp/_ph0

Notes
- 기존 regular file은 절대 덮어쓰지 않음
- self-link, 중복 link 자동 회피
- EPW 버전 차이에 따른 dvscf naming을 동시에 생성

------------------------------------------------------------
📊 Electronic Structure (Bands / DOS)
------------------------------------------------------------

pdosPlot.py 📊

Purpose  
QE projwfc.x 출력 PDOS 파일들을 자동 수집하여  
Fermi level 정렬 + spin 처리 포함 PDOS overlay plot을 생성합니다.

Usage
python pdosPlot.py

Output
<seedname>_pdos_overlay.png

Notes
- QE output에서 Fermi energy 자동 탐색
- spin down은 음수로 플롯

wannier_plot_check.py 🔍

Purpose  
DFT band와 Wannier-interpolated band를 비교하여  
Wannierization 품질을 시각적으로 검증합니다.

Usage
python wannier_plot_check.py scf_band.dat wannier_band.dat

------------------------------------------------------------
📈 Phonon Analysis
------------------------------------------------------------

plot_phonon_dispersion.py 📈

Purpose  
QE matdyn.x 결과로부터 phonon dispersion plot을 생성합니다.

Usage
python plot_phonon_dispersion.py --freq matdyn.freq.gp

------------------------------------------------------------
🧮 EPW Post-processing
------------------------------------------------------------

postprocess_epw.py 🧮

Purpose  
EPW 계산 결과(λ, α²F, Tc 등)를 후처리합니다.

Usage
python postprocess_epw.py

------------------------------------------------------------
🚀 Automation / SLURM
------------------------------------------------------------

postprocess.sh 🚀

Purpose  
SLURM 환경에서 post-processing 스크립트를 일괄 실행합니다.

Usage
sbatch postprocess.sh

==============================
English 🇬🇧
==============================

### 📚 Table of Contents

#### 🔗 EPW / PH Utilities
- link_ph_outputs.py

#### 📊 Electronic Structure (Bands / DOS)
- pdosPlot.py
- wannier_plot_check.py

#### 📈 Phonon Analysis
- plot_phonon_dispersion.py

#### 🧮 EPW Post-processing
- postprocess_epw.py

#### 🚀 Automation / SLURM
- postprocess.sh

------------------------------------------------------------
🔗 EPW / PH Utilities
------------------------------------------------------------

link_ph_outputs.py 🔗

Purpose  
Create safe symbolic links for PH outputs so EPW can locate phonon files.

Usage
python link_ph_outputs.py PREFIX
python link_ph_outputs.py PREFIX --dvscf_dir ./tmp/_ph0

------------------------------------------------------------
📊 Electronic Structure (Bands / DOS)
------------------------------------------------------------

pdosPlot.py 📊

Purpose  
Generate PDOS overlay plots with Fermi-level alignment and spin handling.

Usage
python pdosPlot.py

wannier_plot_check.py 🔍

Purpose  
Sanity check comparing DFT and Wannier-interpolated bands.

Usage
python wannier_plot_check.py scf_band.dat wannier_band.dat

------------------------------------------------------------
📈 Phonon Analysis
------------------------------------------------------------

plot_phonon_dispersion.py 📈

Purpose  
Generate headless phonon dispersion plots from QE outputs.

Usage
python plot_phonon_dispersion.py --freq matdyn.freq.gp

------------------------------------------------------------
🧮 EPW Post-processing
------------------------------------------------------------

postprocess_epw.py 🧮

Purpose  
Post-process EPW results such as λ, α²F, and Tc.

Usage
python postprocess_epw.py

------------------------------------------------------------
🚀 Automation / SLURM
------------------------------------------------------------

postprocess.sh 🚀

Purpose  
Run post-processing steps in batch mode via SLURM.

Usage
sbatch postprocess.sh
