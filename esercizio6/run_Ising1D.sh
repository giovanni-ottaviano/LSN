#!/bin/bash

FLOAT='^[0-9]+([.][0-9]+)?$'

if [[ $# -ne 4 ]] || ! [[ $1 =~ [0-2] ]] || ! [[ $2 =~ $FLOAT ]] || [[ $3 = '^[0-9]+$' ]] || [[ $4 = '^[0-9]+$' ]]; then
  echo -e "\nUsage [options 0, 1, 2] [Temperature] [N_blocks] [Steps per block]"
  echo -e "0 = start new ISING simulation, 1 = restart from final configuration"
  echo -e "2 = equilibrate system and run simulation\n"
  exit 1
fi

ISINGEXE=Monte_Carlo_ISING_1D
INPUT=input.dat
DUMMY=input.dummy
TEMP=$2
BLK=$3
STEPS=$4
EQSTEPS=500
TEMPINCR=0.025

#if you want to start e new simulation, clean output files and initialize config
if [[ $1 -eq 0 ]]; then
    rm output.* seed.out config.final config.0
    sed -e "s/restartdummy/0/" \
	-e "s/tempdummy/$TEMP/" \
	-e "s/blocksdummy/$BLK/" \
	-e "s/stepsdummy/$STEPS/" < $DUMMY > $INPUT   #start new simulation, set temperature, number of blocks and steps 
    ./$ISINGEXE                                    #run ISING
elif [[ $1 -eq 1 ]]; then
    mv config.final config.0
    rm output.*
    sed -e "s/restartdummy/1/" \
	-e "s/tempdummy/$TEMP/" \
        -e "s/blocksdummy/$BLK/" \
        -e "s/stepsdummy/$STEPS/" < $DUMMY > $INPUT  #restart, set temperature, number of blocks and steps
    ./$ISINGEXE                                     #run ISING
elif [[ $1 -eq 2 ]]; then
    rm -f output.*	#remove outputs

    for (( i = 0; i < 61; i++ )); do
	rm -f seed.out config.final config.0

        NEWTEMP=$(awk -v T=$TEMP -v incr=$TEMPINCR -v step=$i 'BEGIN {printf T+incr*step; exit}')
	sed -e "s/restartdummy/0/" \
            -e "s/tempdummy/$NEWTEMP/" \
            -e "s/blocksdummy/1/" \
            -e "s/stepsdummy/500/" < $DUMMY > $INPUT  #prepare equilibration
	./$ISINGEXE                                 #run equilibration

	mv config.final config.0
	rm output.*.0
        sed -e "s/restartdummy/1/" \
            -e "s/tempdummy/$NEWTEMP/" \
            -e "s/blocksdummy/$BLK/" \
            -e "s/stepsdummy/$STEPS/" < $DUMMY > $INPUT  #prepare simulation
        ./$ISINGEXE                                #run ISING simulation

        mv output.ene.0 output.ene.$((i+1))        #change output names
        mv output.heat.0 output.heat.$((i+1))
        mv output.mag.0 output.mag.$((i+1))
        mv output.chi.0 output.chi.$((i+1))

	echo -e "\nSimulation at temperature = $NEWTEMP ... DONE\n\n"
    done
fi

exit 0
