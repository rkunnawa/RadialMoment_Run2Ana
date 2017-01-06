if [ $# -ne 1 ]
then
    echo "Usage: ./compileWithROOT.sh <inputFile>"
    exit 1
fi

path=$PWD

g++ $1 $(fastjet-config --cxxflags --libs) $(root-config --cflags --libs) -lRecursiveTools  -I ${CMSSW_BASE}/src -I ${CMSSW_RELEASE_BASE}/src -L ${CMSSW_BASE}/lib/${SCRAM_ARCH} -L ${CMSSW_RELEASE_BASE}/lib/${SCRAM_ARCH} -lCondFormatsJetMETObjects -Wall -Wextra -I $path -O2 -o "${1/%.C/}.exe"
