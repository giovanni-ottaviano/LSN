# !/bin/bash

if [[ $# -ne 3  ]] || ! [[ $1 =~ [0-3] ]] || ! [[ $2 =~ ^[0-9]+$ ]] || ! [[ $3 =~ ^[0-9]+$ ]]; then
  echo -e "\nUsage [options 0, 1, 2, 3] [Number of blocks] [Number of steps]"
  echo -e "0 = start new MD simulation, 1 = restart from final configuration"
  echo -e "2 = new + equilibration, 3 = run after equilibration (same as 1 without moving config.final to config.0)\n"
  exit 1
fi

MDP=MolDyn_NVE
CLEAN=clean.sh
INPUT=input.dat
EQSTEPS=1000
STEPS=$3
BLOCKS=$2

#if you want tostart e new simulation, clean output files and initialize config
if [[ $1 -eq 0  ]]; then
    bash $CLEAN
    cp config.fcc config.0
    sed -e "s:restartdummy:0:" < input.dummy > $INPUT   #new start
    sed -i "s:stepsdummy:$STEPS:" $INPUT   #number of steps
    sed -i "s:blocksdummy:$BLOCKS:" $INPUT
    ./$MDP
elif [[ $1 -eq 1 ]]; then
    bash $CLEAN
    mv config.final config.0
    mv old.final old.0
    sed -e "s:restartdummy:1:" < input.dummy > $INPUT   #restart
    sed -i "s:stepsdummy:$STEPS:" $INPUT   #number of steps
    sed -i "s:blocksdummy:$BLOCKS:" $INPUT
    ./$MDP
elif [[ $1 -eq 2 ]]; then
    bash $CLEAN
    cp config.fcc config.0
    sed -e "s:restartdummy:0:" < input.dummy > $INPUT   #new start
    sed -i "s:stepsdummy:$EQSTEPS:" $INPUT   #steps for equilibration
    sed -i "s:blocksdummy:1:" $INPUT	#1 block for equilibration
    ./$MDP
    sed -e "s:restartdummy:1:" < input.dummy > $INPUT   #restart
    sed -i "s:stepsdummy:$EQSTEPS:" $INPUT   #steps for equilibration
    sed -i "s:blocksdummy:1:" $INPUT

#i < 5 -> liquid/gas, i < 2 -> solid
    for (( i = 0; i < 10; i++ )); do
        mv config.final config.0
        mv old.final old.0
        ./$MDP
	echo -e "\n----------------------"
	echo -e "Equilibration STEP $((i+1))"
	echo -e "DONE"
	echo -e "----------------------"
    done

    mv config.final config.0
    mv old.final old.0	#equilibration completed: ready to start MD simulation
else
    bash $CLEAN
    sed -e "s:restartdummy:1:" < input.dummy > $INPUT   #restart
    sed -i "s:stepsdummy:$STEPS:" $INPUT   #number of steps
    sed -i "s:blocksdummy:$BLOCKS:" $INPUT
    ./$MDP
fi

exit 0
