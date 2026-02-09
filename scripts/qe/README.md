# QE / EPW Utility Scripts ⚙️🧪

Utilities for **Quantum ESPRESSO (QE)** and **EPW** workflows.  
Tested with **QE 7.5** and **EPW 6.0**.

---

## 🌐 Language / 언어
- 🇰🇷 [한국어](#korean)
- 🇬🇧 [English](#english)

---

<a id="korean"></a>
## 한국어 🇰🇷
---

### 📚 Table of Contents

#### 🧮 EPW
- [link_ph_outputs.py](#link_ph_outputspy-ko)
- [postprocess_epw.py](#postprocess_epwpy-ko)
- [wannier_plot_check.py](#wannier_plot_checkpy-ko)

#### 📊 Electronic Structure (Bands / DOS)
- [pdosPlot.py](#pdosplotpy-ko)

#### 📈 Phonon Analysis
- [plot_phonon_dispersion.py](#plot_phonon_dispersionpy-ko)

#### 🚀 Automation / SLURM
- [postprocess.sh](#postprocesssh-ko)

---

## 🧮 EPW

### link_ph_outputs.py 🔗 <a id="link_ph_outputspy-ko"></a>

**Purpose**  
EPW가 요구하는 phonon 파일 naming 규칙에 맞추기 위해  
PH 계산 결과(`dyn`, `dvscf`)에 대해 **안전한 symbolic link**를 생성합니다.

Quantum ESPRESSO 공식 `pp.py`가 환경 의존적으로 정상 동작하지 않는 문제를
회피하기 위한 **실용적인 대체 스크립트**입니다.

**Generated names**
- `PREFIX.dyn_qN`
- `PREFIX.dvscf_qN`
- `PREFIX.dvscfN_1`

**Where to run**
- `PREFIX.dyn*` 파일이 존재하는 디렉토리 (`ph.x` 실행 디렉토리)

**Usage**
```bash
python link_ph_outputs.py PREFIX
python link_ph_outputs.py PREFIX --dvscf_dir ./tmp/_ph0
```

**Notes**
- 기존 regular file 절대 overwrite 하지 않음
- self-link / 중복 link 자동 회피
- EPW 버전 차이 대응 dvscf naming 동시 생성

---

### postprocess_epw.py 🧮 <a id="postprocess_epwpy-ko"></a>

**Purpose**  
EPW 계산 결과(λ, α²F, Tc 등)를 자동 후처리합니다.

**Usage**
```bash
python postprocess_epw.py
```

**Notes**
- EPW output 파일이 현재 디렉토리에 있어야 함
- SLURM batch 이후 실행 용도

---

### wannier_plot_check.py 🔍 <a id="wannier_plot_checkpy-ko"></a>

**Purpose**  
DFT band와 Wannier-interpolated band를 비교하여  
Wannierization 품질을 시각적으로 검증합니다.

**Usage**
```bash
python wannier_plot_check.py scf_band.dat wannier_band.dat
```

**Notes**
- 정량 분석 목적 아님 (sanity check)

---

## 📊 Electronic Structure (Bands / DOS)

### pdosPlot.py 📊 <a id="pdosplotpy-ko"></a>

**Purpose**  
QE `projwfc.x` PDOS 출력 파일을 자동 수집하여  
Fermi level 정렬 + spin 처리 포함 PDOS overlay plot을 생성합니다.

**Usage**
```bash
python pdosPlot.py
```

**Output**
```
<seedname>_pdos_overlay.png
```

**Notes**
- QE output에서 Fermi energy 자동 탐색
- spin ↑ : positive / spin ↓ : negative

---

## 📈 Phonon Analysis

### plot_phonon_dispersion.py 📈 <a id="plot_phonon_dispersionpy-ko"></a>

**Purpose**  
QE `matdyn.x` 결과(`.freq.gp` 또는 `.freq`)로부터  
GUI 없이 phonon dispersion plot을 생성합니다.

**Usage**
```bash
python python plot_phonon_dispersion.py --freq TiS2_monolayer
python plot_phonon_dispersion.py --freq TiS2_monolayer --dos TiS2_monolayer.phdos --qpath qpath.in
python plot_phonon_dispersion.py --freq TiS2_monolayer --unit thz --emin -10 --emax 80
```

**Notes**
- `.freq.gp` 없으면 raw `.freq`도 파싱 가능

---

## 🚀 Automation / SLURM

### postprocess.sh 🚀 <a id="postprocesssh-ko"></a>

**Purpose**  
SLURM 환경에서 EPW 후처리 스크립트를 일괄 실행합니다.

**Usage**
```bash
sbatch postprocess.sh
```

**Notes**
- SLURM directive는 클러스터 환경에 맞게 수정 필요

---

<a id="english"></a>
## English 🇬🇧
---

### 📚 Table of Contents

#### 🧮 EPW
- [link_ph_outputs.py](#link_ph_outputspy-en)
- [postprocess_epw.py](#postprocess_epwpy-en)
- [wannier_plot_check.py](#wannier_plot_checkpy-en)

#### 📊 Electronic Structure (Bands / DOS)
- [pdosPlot.py](#pdosplotpy-en)

#### 📈 Phonon Analysis
- [plot_phonon_dispersion.py](#plot_phonon_dispersionpy-en)

#### 🚀 Automation / SLURM
- [postprocess.sh](#postprocesssh-en)

---

## 🧮 EPW

### link_ph_outputs.py 🔗 <a id="link_ph_outputspy-en"></a>

**Purpose**  
Creates **safe symbolic links** for PH calculation outputs (`dyn`, `dvscf`)  
to satisfy EPW-required phonon file naming conventions.

This script is a **practical replacement** for the official Quantum ESPRESSO
`pp.py`, which may behave inconsistently depending on the environment.

**Generated names**
- `PREFIX.dyn_qN`
- `PREFIX.dvscf_qN`
- `PREFIX.dvscfN_1`

**Where to run**
- Directory containing `PREFIX.dyn*` files (where `ph.x` finished)

**Usage**
```bash
python link_ph_outputs.py PREFIX
python link_ph_outputs.py PREFIX --dvscf_dir ./tmp/_ph0
```

**Notes**
- Never overwrites existing regular files
- Automatically avoids self-links and duplicate links
- Generates multiple dvscf naming variants for EPW version differences

---

### postprocess_epw.py 🧮 <a id="postprocess_epwpy-en"></a>

**Purpose**  
Post-processes EPW outputs such as λ, α²F, and Tc.

**Usage**
```bash
python postprocess_epw.py
```

**Notes**
- EPW output files must be in the current directory
- Intended for execution after SLURM jobs

---

### wannier_plot_check.py 🔍 <a id="wannier_plot_checkpy-en"></a>

**Purpose**  
Visually compares DFT and Wannier-interpolated band structures to validate
Wannierization quality.

**Usage**
```bash
python wannier_plot_check.py scf_band.dat wannier_band.dat
```

**Notes**
- Intended as a sanity check, not quantitative analysis

---

## 📊 Electronic Structure (Bands / DOS)

### pdosPlot.py 📊 <a id="pdosplotpy-en"></a>

**Purpose**  
Automatically collects QE `projwfc.x` PDOS outputs and generates
spin-resolved, Fermi-aligned PDOS overlay plots.

**Usage**
```bash
python pdosPlot.py
```

**Output**
```
<seedname>_pdos_overlay.png
```

**Notes**
- Automatically detects Fermi energy from QE output
- spin-up plotted positive, spin-down plotted negative

---

## 📈 Phonon Analysis

### plot_phonon_dispersion.py 📈 <a id="plot_phonon_dispersionpy-en"></a>

**Purpose**  
Generates phonon dispersion plots from QE `matdyn.x` outputs
(`.freq.gp` or raw `.freq`) in headless environments.

**Usage**
```bash
python plot_phonon_dispersion.py --freq matdyn.freq.gp
python plot_phonon_dispersion.py --freq matdyn.freq.gp --qpath qpath.in
python plot_phonon_dispersion.py --freq matdyn.freq.gp --unit mev --emin -10 --emax 80
```

**Notes**
- Uses Matplotlib Agg backend (no GUI required)
- Automatically marks Γ / M / K high-symmetry points
- Can parse raw `.freq` if `.freq.gp` is unavailable

---

## 🚀 Automation / SLURM

### postprocess.sh 🚀 <a id="postprocesssh-en"></a>

**Purpose**  
Wrapper script to batch-run EPW post-processing steps in SLURM environments.

**Usage**
```bash
sbatch postprocess.sh
```

**Notes**
- SLURM directives must be adapted to the target cluster
