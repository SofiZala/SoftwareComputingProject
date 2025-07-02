#!/bin/bash

# Check to have 1 input number
if [ $# -ne 1 ]; then
    echo "Usage: ./run_pipeline.sh [0|1]"
    exit 1
fi

# Check to have 0 or 1
if [ "$1" != "0" ] && [ "$1" != "1" ]; then
    echo "Invalid argument: $1"
    echo "Usage: ./run_pipeline.sh [0|1]"
    exit 1
fi

#Save number in Brem
brem=$1

echo ">> Starting full pipeline with brem = $brem"

########################################
# 1. Setting
########################################
echo ">> Entering Setting/"
cd Setting || exit 1
echo ">> Running SetInput.cpp with brem = $brem"
root -l -b -q "SetInput.cpp($brem)"
cd .. || exit 1

########################################
# 2. Simulation
########################################
echo ">> Entering Simulation/"
cd Simulation || exit 1
echo ">> Running simulation.sh with brem = $brem"
./simulation.sh "$brem"
cd .. || exit 1

########################################
# 3. PreliminaryFit
########################################
echo ">> Entering PreliminaryFit/Brem$brem/"
cd PreliminaryFit/Brem$brem || exit 1
echo ">> Running fit.cpp"
root -l -b -q fit.cpp
cd ../.. || exit 1

########################################
# 4. FitToData
########################################
echo ">> Entering FitToData/Brem$brem/"
cd FitToData/Brem$brem || exit 1
echo ">> Running finalFit.cpp"
root -l -q finalFit.cpp
cd ../.. || exit 1

echo ">> Pipeline completed successfully."
