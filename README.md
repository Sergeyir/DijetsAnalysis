# Overview

This is a simple project that helps to study dijets in high-energy p+p collisions at LO. It compares $d \sigma/dp_{T}$ and $d \sigma / d \Delta y$ analytical computation with single parton in pythia8 and jets in pythia8 + fastjet3. Also the pdf set NNPDF31_lo_as_0118 is used. The realisation of comutation and generation is implemented with c++.

Using $d \sigma / dy_{1} dy_{2} dp^2_{T}$ from [PDG report](https://pdg.lbl.gov/2023/reviews/contents_sports.html) (Eq.51.42) the following formula was obtained

```math
\frac{d \sigma}{dp_T dy_1 dy_2} = \frac{8 \pi p_T}{\hat{s}} \sum_{ijkl} x_{1} f_{i}(x_{1}, \mu_F^2) x_{2} f_{j}(x_{2}, \mu_F^2) \frac{d \sigma_{ij \rightarrow kl}}{d \Omega}
```

The calculation in this project is implemented on the hard scale i.e. $\mu_F \approx p_T$

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

Kynematic variables can be calculated using formulas from [Michelangelo L. Mangano "Introduction to QCD"](https://cds.cern.ch/record/454171/files/p53.pdf) (Eq.206-211)

# Requirements

- GNU GCC wiht C++17 or newer
- [ROOT6](https://root.cern/) + compilation with python3
- [LHAPDF6](https://lhapdf.hepforge.org/) + pdf set NNPDF31_lo_as_0118 installed
- [PYTHIA8](https://pythia.org/) + compilation with LHAPDF6
- [FASTJET3](https://fastjet.fr/) 

# Installation

```sh
git clone https://github.com/Sergeyir/DijetsAnalysis
```

# Usage

First define paths `$PYTHIA`, `$FASTJET`, `$ROOT`, and `$LHAPDF` in your .bashrc file to the packages directories that they were compiled into. Also check `${LHAPDF6}/../src` - part of `$LHAPDF6_LIB` variable in Makefile.inc since the directories can mismatch.

There are 2 programs you can run: generate.cpp and analytic.cpp. First calculates the cross sections in pythia8 for partons and for jets. An the second one calculates the cross sections analyticaly by implementing monte-carlo integration. You need to run them first to get the data that you then can draw.

Frist make all programs
```sh
make analytic generate
```

To launch any program type while substituting 'name' by the name of the program

```sh
./'name'.exe
```

All input parameters are located in the top of .cpp files in struct called Par. Also there are input files that pass parameters for pythia generation in input directory.

After generating the data you can draw the result by running
```sh
python draw.py
```
Pictures will be drawn and saved in the output directory
