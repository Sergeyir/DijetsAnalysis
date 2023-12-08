# Install directory prefixes.
PREFIX_BIN=${PYTHIA}/bin
PREFIX_INCLUDE=${PYTHIA}/include
PREFIX_LIB=${PYTHIA}/lib
PREFIX_SHARE=${PYTHIA}/share/Pythia8

CXX=g++
CXX_COMMON=-Wpedantic -W -Wall -Werror
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

FASTJET3_USE=1
FASTJET3_CONFIG=${FASTJET}/src
FASTJET3_BIN=${FASTJET}/bin
FASTJET3_INCLUDE=-I${FASTJET}/include
FASTJET3_LIB=-Wl,-rpath,${FASTJET}/lib -lfastjettools -L${FASTJET}/lib -lfastjet -lm

ROOT_USE=1
ROOT_CONFIG=${ROOT}/install/bin/root-config
ROOT_BIN=${ROOT}/install/bin
ROOT_INCLUDE=${ROOT}/install/include
ROOT_LIB=-Wl,-rpath,${ROOT}/install/lib -L${ROOT}/install/lib -lCore

LHAPDF6_USE=1
LHAPDF6_CONFIG=${LHAPDF}bin/lhapdf-config
LHAPDF6_BIN=${LHAPDF}/bin
LHAPDF6_INCLUDE=-I${LHAPDF}/include
LHAPDF6_LIB=-Wl,-rpath,${LHAPDF}/lib -L${LHAPDF}/lib -lLHAPDF -L${LHAPDF}/../src
