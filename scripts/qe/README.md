# Scripts CLI Terminal Usage for QE / EPW Utility 
Utilities for **Quantum ESPRESSO (QE)** and **EPW** workflows:
- EPW symlink preparation
- PDOS overlay plotting
- Phonon dispersion plotting
- EPW post-processing
- Wannier band sanity checks
---

## Language / 언어
- [한국어](#한국어)
- [English](#english)

---

## 한국어

### Table of Contents
- [`link_ph_outputs.py`](#link_ph_outputspy-ko)
- [`pdosPlot.py`](#pdosplotpy-ko)
- [`plot_phonon_dispersion.py`](#plot_phonon_dispersionpy-ko)
- [`postprocess_epw.py`](#postprocess_epwpy-ko)
- [`wannier_plot_check.py`](#wannier_plot_checkpy-ko)
- [`postprocess.sh`](#postprocesssh-ko)

---

### `link_ph_outputs.py` <a id="link_ph_outputspy-ko"></a>

**Purpose**  
EPW가 요구하는 phonon 파일 naming 규칙에 맞추기 위해  
PH 계산 결과(`dyn`, `dvscf`)에 대해 **안전한 symbolic link**를 생성합니다.

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
- 기존 regular file은 절대 덮어쓰지 않음
- self-link, 중복 link 자동 회피
- EPW 버전 차이에 따른 dvscf naming을 동시에 생성

---

### `pdosPlot.py` <a id="pdosplotpy-ko"></a>

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
```
<seedname>_pdos_overlay.png
```

**Notes**
- QE output에서 `the Fermi energy is ... eV` 자동 탐색
- spin-polarized 계산 시:
  - spin ↑ : positive
  - spin ↓ : plotted negative
- 별도 인자 없이 완전 자동 동작

---

### `plot_phonon_dispersion.py` <a id="plot_phonon_dispersionpy-ko"></a>

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
python plot_phonon_dispersion.py \
  --freq matdyn.freq.gp \
  --qpath qpath.in \
  --out phonon.png
```

**Unit / range control**
```bash
python plot_phonon_dispersion.py \
  --freq matdyn.freq.gp \
  --unit mev \
  --emin -10 --emax 80
```

**Notes**
- Matplotlib `Agg` backend 사용 (GUI 불필요)
- Γ / M / K 등 high-symmetry point 자동 표시
- `.freq.gp`가 없을 경우 raw `.freq`도 파싱 가능

---

### `postprocess_epw.py` <a id="postprocess_epwpy-ko"></a>

**Purpose**  
EPW 계산 결과(λ, α²F, Tc 등)를 후처리하기 위한 스크립트입니다.

**Typical usage**
```bash
python postprocess_epw.py
```

**Notes**
- EPW output 파일들이 현재 디렉토리에 있어야 함
- SLURM 계산 이후 batch 후처리 단계용

---

### `wannier_plot_check.py` <a id="wannier_plot_checkpy-ko"></a>

**Purpose**  
DFT band와 Wannier-interpolated band를 비교하여  
**Wannierization 품질을 시각적으로 검증**합니다.

**Typical usage**
```bash
python wannier_plot_check.py
```

**Notes**
- DFT band 데이터와 Wannier band 데이터가 모두 필요
- 정량 분석용이 아니라 sanity check 목적

---

### `postprocess.sh` <a id="postprocesssh-ko"></a>

**Purpose**  
SLURM 환경에서 post-processing 스크립트들을  
일괄 실행하기 위한 wrapper 스크립트입니다.

**Usage**
```bash
sbatch postprocess.sh
```

**Notes**
- SLURM directive는 클러스터 환경에 맞게 수정 필요

---

## English

### Table of Contents
- [`link_ph_outputs.py`](#link_ph_outputspy-en)
- [`pdosPlot.py`](#pdosplotpy-en)
- [`plot_phonon_dispersion.py`](#plot_phonon_dispersionpy-en)
- [`postprocess_epw.py`](#postprocess_epwpy-en)
- [`wannier_plot_check.py`](#wannier_plot_checkpy-en)
- [`postprocess.sh`](#postprocesssh-en)

---

### `link_ph_outputs.py` <a id="link_ph_outputspy-en"></a>

**Purpose**  
Create **safe symbolic links** so EPW can find phonon outputs under
the naming conventions it expects:
- `PREFIX.dyn_qN`
- `PREFIX.dvscf_qN`
- `PREFIX.dvscfN_1`

**Where to run**  
- In the directory containing `PREFIX.dyn*` files  
  (typically the PH calculation directory)

**Basic usage**
```bash
python link_ph_outputs.py PREFIX
```

**Explicit dvscf directory**
```bash
python link_ph_outputs.py PREFIX --dvscf_dir ./tmp/_ph0
```

**Notes**
- Never overwrites an existing regular file
- Avoids self-referential and duplicate links
- Generates multiple dvscf naming variants for EPW compatibility

---

### `pdosPlot.py` <a id="pdosplotpy-en"></a>

**Purpose**  
Automatically scan QE `projwfc.x` PDOS outputs and generate a
**single overlay PDOS plot** with:
- Fermi-level detection and optional shifting
- Proper spin handling (spin-down plotted negative)

**Where to run**  
- In the directory containing `*.pdos_*` files and QE `*.out` outputs

**Basic usage**
```bash
python pdosPlot.py
```

**Output**
```
<seedname>_pdos_overlay.png
```

**Notes**
- Extracts EF from lines like `the Fermi energy is ... eV`
- Heuristically skips EF-shift if PDOS already appears EF-centered
- Fully automatic (no CLI arguments)

---

### `plot_phonon_dispersion.py` <a id="plot_phonon_dispersionpy-en"></a>

**Purpose**  
Headless phonon dispersion plotter for QE `matdyn.x` outputs:
- supports `.freq.gp` and raw `.freq`
- optional high-symmetry labels from QE-style q-path files
- PNG output only (Agg backend)

**Where to run**  
- In the directory containing `.freq.gp` or `.freq`

**Basic usage (recommended)**
```bash
python plot_phonon_dispersion.py --freq matdyn.freq.gp
```

**Notes**
- No GUI required
- Draws high-symmetry vertical lines and tick labels when `--qpath` is provided

---

### `postprocess_epw.py` <a id="postprocess_epwpy-en"></a>

**Purpose**  
Post-process EPW outputs (e.g. λ, α²F, Tc).

**Typical usage**
```bash
python postprocess_epw.py
```

**Notes**
- Requires EPW outputs to be present in the current directory
- Intended for batch post-processing after EPW runs

---

### `wannier_plot_check.py` <a id="wannier_plot_checkpy-en"></a>

**Purpose**  
Quick visual sanity check comparing DFT bands and Wannier-interpolated bands.

**Typical usage**
```bash
python wannier_plot_check.py
```

**Notes**
- Requires both DFT band data and Wannier band data

---

### `postprocess.sh` <a id="postprocesssh-en"></a>

**Purpose**  
SLURM wrapper script to run post-processing jobs in batch mode.

**Usage**
```bash
sbatch postprocess.sh
```

**Notes**
- Edit SLURM directives to match your cluster configuration
