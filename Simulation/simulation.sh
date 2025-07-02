#!/bin/bash

mode=$1

# this code runs the RapidSim simulation for the D02Kpiee and D02Kpipienu or D02Kpipipi processes
# 50000 events are generated for each process

if [ "$mode" == "1" ]; then
  echo "Esecuzione in modalità 1: Kpiee + Kpipienu"
  
  cd D02Kpiee || exit
  /opt/RapidSim/bin/RapidSim.exe Kpiee 50000 1
  cd ..

  cd D02Kpipienu || exit
  /opt/RapidSim/bin/RapidSim.exe Kpipienu 50000 1
  cd ..

elif [ "$mode" == "0" ]; then
  echo "Esecuzione in modalità 0: Kpiee + Kpipipi"
  
  cd D02Kpiee || exit
  /opt/RapidSim/bin/RapidSim.exe Kpiee 50000 1
  cd ..

  cd D02Kpipipi || exit
  /opt/RapidSim/bin/RapidSim.exe Kpipipi 50000 1
  cd ..

fi
