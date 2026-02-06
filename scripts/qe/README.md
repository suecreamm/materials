# QE / EPW Utility Scripts ⚙️🧪

Utilities for **Quantum ESPRESSO (QE)** and **EPW** workflows.  
Tested with QE version 7.5 and EPW version 6.0.

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

이 스크립트는 **Quantum ESPRESSO 공식 패키지에 포함된 `pp.py`가
현재 사용 중인 계산 환경에서 정상적으로 동작하지 않아**,
이를 대체하기 위한 대안으로 작성되었습니다.  
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
to conform to the **phonon file naming conventions required by EPW**.

This script was written as a **practical replacement** for the official  
Quantum ESPRESSO `pp.py`, which does **not work reliably in the current environment**.  
Its purpose is to provide a **more stable and environment-independent solution**  
for EPW-required file naming and directory handling.

Generated names:
- `PREFIX.dyn_qN`
- `PREFIX.dvscf_qN`
- `PREFIX.dvscfN_1`

**Where to run**  
- A directory containing `PREFIX.dyn*` files  
  (typically the directory where `ph.x` finished)

**Basic usage**
```bash
python link_ph_outputs.py PREFIX
```

**Explicit dvscf directory**
```bash
python link_ph_outputs.py PREFIX --dvscf_dir ./tmp/_ph0
```

**Notes**
- ❌ Never overwrites existing regular files
- 🔁 Automatically avoids self-links and duplicate links
- ✅ Generates multiple dvscf naming variants to support EPW version differences
- 🔄 Designed as a fallback when the official `pp.py` behaves inconsistently
