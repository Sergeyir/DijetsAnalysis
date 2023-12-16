# Overview

This is a simple project that helps to study dijets in high-energy p+p collisions. It compares $d \sigma/dp_{T}$ for single jets and $d \sigma / d \Delta y$ for dijets for pdf set NNPDF31_lo_as_0118 in pythia8 generation with and without hadronization and radiation + fastjet3 reconstruction algorithm and basic computation. The realisation of comutation and generation is implemented with c++.

Using $d \sigma / dy_{1} dy_{2} dp^2_{T}$ from [Michelangelo L. Mangano "Introduction to QCD"](https://cds.cern.ch/record/454171/files/p53.pdf) (Eq.205) the following formula was obtained

```math
\frac{d \sigma}{dp_T dy_1 dy_2} = \frac{8 \pi p_T}{s} \sum_{ijkl} f_{ij}(x_{1}, \mu^2) f_{ij}(x_{2}, \mu^2) \frac{d \sigma_{ij \rightarrow kl}}{d \Omega}
```

The summation is performed by the simplest processes ($q$ is for quark, $g$ is for gluon):
- $q + q' \rightarrow q + q'$
- $q + q \rightarrow q + q$
- $q + \bar{q} \rightarrow q + \bar{q}$
- $q + \bar{q} \rightarrow q + '\bar{q}'$
- $q + \bar{q} \rightarrow g + g$
- $g+g \rightarrow q + \bar{q}$
- $g+q \rightarrow g + q$
- $g+g \rightarrow g + g$

You can find $d \sigma_{ij \rightarrow kl}/d \Omega$ in [PDG report](https://pdg.lbl.gov/2023/reviews/contents_sports.html) (Eq. 51.4 - 51.12)

The final formulas for the cross sections are
```math
\frac{d \sigma}{dp_T} = \int dy_1 \int dy_2 \frac{d \sigma}{dp_T dy_1 dy_2}
```

```math
\frac{d \sigma}{d \Delta y} = \int dp_T \int dy_1 \int dy_2 \frac{d \sigma}{dp_T dy_1 dy_2} \delta(\Delta y - |y_1 - y_2|)
```

# Requirements

- GNU GCC wiht C++17 or newer
- [ROOT](https://root.cern/) + compilation with python3
- [LHAPDF6](https://lhapdf.hepforge.org/) + pdf set NNPDF31_lo_as_0118 installed
- [PYTHIA8](https://pythia.org/) + compilation with LHAPDF6
- [FASTJET3](https://fastjet.fr/) 

# Installation

```sh
git clone https://github.com/Sergeyir/DijetsAnalysis
```

# Usage

First define paths `$PYTHIA`, `$FASTJET`, `$ROOT`, and `$LHAPDF` in your .bashrc file to the packages directories that they were compiled into. Also check `${LHAPDF6}/../src` - part of `$LHAPDF6_LIB` variable in Makefile.inc since the directories can mismatch.

There are 2 programs you can run: generate.cpp and analytic.cpp. First calculates the cross sections in pythia8 with fastjet3 algorithm and writes the output file. Second calculates the cross section by implementing monte-carlo integration and also writes the output file. You need to run them first to get the data you then can draw.

Frist make both programs
```sh
make analytic generate
```

Then you can launch pythia8+fastjet3 computation by typing
```sh
./generate.exe setup_name
```

For setup_name pass a name of a .cmnd file located in input directory without an extension. There are 2 files already presented in the input directory:
- setup1 - for calculation of cross sections with normal pythia8 setup.
- setup2 - for calculation of cross sections without hadronization and without initial and final state radiation
These files also contain important information for generation and contain parameters for energy, pdf set, pTHatMin, number of events to generate, $\Delta y_{max}$, and fastjet R parameter.

And you also can launch the analytic computation of the cross section by typing
```sh
./analytic.exe
```

After generating the data you can draw the result by running
```sh
python draw.py
```
