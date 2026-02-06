# QE / EPW Utility Scripts ⚙️🧪

Utilities for **Quantum ESPRESSO (QE)** and **EPW** workflows.
Tested with QE version 7.5 and EPW version 6.0.

---

## 🌐 Language / 언어
- 🇰🇷 한국어
- 🇬🇧 English

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

이 스크립트는 **Quantum ESPRESSO 공식 패키지에 포함된 `pp.py`가
현재 사용 중인 계산 환경에서 정상적으로 동작하지 않아**,
이를 대체하기 위한 **실용적인 대안**으로 작성되었습니다.  
EPW 워크플로우에서 요구되는 파일 naming 및 디렉토리 구조를
보다 안정적으로 처리하는 데 목적이 있습니다.

생성되는 이름:
- `PREFIX.dyn_qN`
- `PREFIX.dvscf_qN`
- `PREFIX.dvscfN_1`

**Where to run**  
- `PREFIX.dyn*` 파일들이 존재하는 디렉토리  
  (보통 `ph.x` 계산이 끝난 작업 디렉토리)

**Basic usage**
```bash
python link_ph_outputs.py PREFIX
```

**Explicit dvscf directory**
```bash
python link_ph_outputs.py PREFIX --dvscf_dir ./tmp/_ph0
```

**Notes**
- ❌ 기존 regular file은 절대 덮어쓰지 않음
- 🔁 self-link, 중복 link 자동 회피
- ✅ EPW 버전 차이에 따른 dvscf naming을 동시에 생성
- 🔄 공식 `pp.py`가 환경 의존적으로 동작하지 않는 경우를 대비한 대체 스크립트

---

### postprocess_epw.py 🧮 <a id="postprocess_epwpy-ko"></a>

**Purpose**  
EPW 계산 결과(λ, α²F, Tc 등)를 후처리하기 위한 스크립트입니다.

**Typical usage**
```bash
python postprocess_epw.py
```

**Notes**
- 📂 EPW output 파일들이 현재 디렉토리에 있어야 함
- 🧵 SLURM 계산 이후 batch 후처리 단계용

---

### wannier_plot_check.py 🔍 <a id="wannier_plot_checkpy-ko"></a>

**Purpose**  
DFT band와 Wannier-interpolated band를 비교하여  
**Wannierization 품질을 시각적으로 검증**합니다.

**Typical usage**
```bash
python wannier_plot_check.py scf_band.dat wannier_band.dat
```

**Notes**
- 📐 DFT band 데이터와 Wannier band 데이터가 모두 필요
- ⚠️ 정량 분석이 아니라 sanity check 목적

---

## 📊 Electronic Structure (Bands / DOS)

### pdosPlot.py 📊 <a id="pdosplotpy-ko"></a>

**Purpose**  
QE `projwfc.x` 출력 PDOS 파일들을 자동 수집하여  
**Fermi level 정렬 + spin 처리 포함 PDOS overlay plot**을 생성합니다.

**Where to run**  
- `*.pdos_*` 파일들과 QE `*.out` 파일이 있는 디렉토리

**Basic usage**
```bash
python pdosPlot.py
```

**Output**
```text
<seedname>_pdos_overlay.png
```

**Notes**
- 🔎 QE output에서 `the Fermi energy is ... eV` 자동 탐색
- 🧲 spin-polarized 계산 시:
  - spin ↑ : positive
  - spin ↓ : plotted negative
- ⚙️ 별도 인자 없이 완전 자동 동작

---

## 📈 Phonon Analysis

### plot_phonon_dispersion.py 📈 <a id="plot_phonon_dispersionpy-ko"></a>

**Purpose**  
QE `matdyn.x` 결과(`.freq.gp` 또는 raw `.freq`)로부터  
GUI 없이(headless) **phonon dispersion plot**을 생성합니다.

**Where to run**  
- `.freq.gp` 또는 `.freq` 파일이 있는 디렉토리

**Basic usage (recommended)**
```bash
python plot_phonon_dispersion.py --freq matdyn.freq.gp
```

**With high-symmetry labels**
```bash
python plot_phonon_dispersion.py   --freq matdyn.freq.gp   --qpath qpath.in   --out phonon.png
```
**qpath.in**
```bash
$ more qpath.in
4
0.0  0.0  0.0    50  G
0.5  0.0  0.0    50  M
0.333333  0.333333  0.0  50  K
0.0  0.0  0.0    50  G
```
**Unit / range control**
```bash
python plot_phonon_dispersion.py   --freq matdyn.freq.gp   --unit mev   --emin -10 --emax 80
```


**Notes**
- 🖥️ Matplotlib `Agg` backend 사용 (GUI 불필요)
- 📍 Γ / M / K 등 high-symmetry point 자동 표시
- 🔄 `.freq.gp`가 없을 경우 raw `.freq`도 파싱 가능

---

## 🚀 Automation / SLURM

### postprocess.sh 🚀 <a id="postprocesssh-ko"></a>

**Purpose**  
SLURM 환경에서 EPW process 후 post-processing 스크립트들을  
일괄 실행하기 위한 wrapper 스크립트입니다.

**Usage**
```bash
sbatch postprocess.sh
```

**Notes**
- 🛠️ SLURM directive는 클러스터 환경에 맞게 수정 필요
- `postprocess_epw.py`와 job submit 스크립트(`7501q.sh`)와 같이 사용함.

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
Create **safe symbolic links** for PH outputs (`dyn`, `dvscf`) so that EPW can locate
phonon-related files using the naming conventions it expects.

This script was written because the **official `pp.py` script distributed with
the Quantum ESPRESSO package does not work reliably in the author’s computing
environment**.  
It is intended as a **practical alternative** that more robustly handles
file naming and directory layouts required by EPW workflows.

Generated filenames:
- `PREFIX.dyn_qN`
- `PREFIX.dvscf_qN`
- `PREFIX.dvscfN_1`

**Where to run**  
- In the directory containing `PREFIX.dyn*` files  
  (typically the working directory where `ph.x` has finished)

**Basic usage**
```bash
python link_ph_outputs.py PREFIX
```

**Explicit dvscf directory**
```bash
python link_ph_outputs.py PREFIX --dvscf_dir ./tmp/_ph0
```

**Notes**
- ❌ Never overwrites an existing regular file
- 🔁 Automatically avoids self-links and duplicate links
- ✅ Generates multiple dvscf naming variants for EPW compatibility
- 🔄 Intended as an alternative when the official `pp.py` does not function properly

---

### postprocess_epw.py 🧮 <a id="postprocess_epwpy-en"></a>

**Purpose**  
Post-process EPW results (e.g. λ, α²F, Tc).

**Typical usage**
```bash
python postprocess_epw.py
```

**Notes**
- 📂 EPW output files must be present in the current directory
- 🧵 Intended for batch post-processing after SLURM jobs finish

---

### wannier_plot_check.py 🔍 <a id="wannier_plot_checkpy-en"></a>

**Purpose**  
Visually compare DFT bands with Wannier-interpolated bands to  
**validate the quality of Wannierization**.

**Typical usage**
```bash
python wannier_plot_check.py scf_band.dat wannier_band.dat
```

**Notes**
- 📐 Requires both DFT band data and Wannier band data
- ⚠️ Intended as a sanity check, not for quantitative analysis

---

## 📊 Electronic Structure (Bands / DOS)

### pdosPlot.py 📊 <a id="pdosplotpy-en"></a>

**Purpose**  
Automatically collect QE `projwfc.x` PDOS outputs and generate a  
**PDOS overlay plot with Fermi-level alignment and proper spin handling**.

**Where to run**  
- In the directory containing `*.pdos_*` files and QE `*.out` outputs

**Basic usage**
```bash
python pdosPlot.py
```

---

## 📈 Phonon Analysis

### plot_phonon_dispersion.py 📈 <a id="plot_phonon_dispersionpy-en"></a>

**Purpose**  
Generate a headless phonon dispersion plot from QE `matdyn.x` results
(`.freq.gp` or raw `.freq`).

**Basic usage**
```bash
python plot_phonon_dispersion.py --freq matdyn.freq.gp
```

---

## 🚀 Automation / SLURM

### postprocess.sh 🚀 <a id="postprocesssh-en"></a>

**Purpose**  
SLURM wrapper script to run EPW-related post-processing steps
in a single batch job.

**Usage**
```bash
sbatch postprocess.sh
```
