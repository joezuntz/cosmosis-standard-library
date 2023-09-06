#!/usr/bin/env bash

set -e

# Bash is insane. If you don't do this next bit then it
# decides that you wanted to pass all the arguments from this
# script on to cosmosis-configure. Madness.
function DoSource() { source cosmosis-configure ; } ; DoSource


# for demo 11 first run
export HALOFIT=takahashi

if [ "$1" == "" ]; then
    demos=$(seq 1 19)
else
    demos=$1
fi

echo "Running demos: " $demos

for i in $demos
do
    if [ "$i" == "3" ]; then
        args="-p grid.nsample_dimension=10"
    elif [ "$i" == "4" ]; then
        args="-v cosmological_parameters.n_s=0.962 cosmological_parameters.h0=0.680 -p maxlike.tolerance=1.0"
    elif [ "$i" == "5" ]; then
        args="-p emcee.samples=30 emcee.nsteps=10 emcee.walkers=20"
    elif [ "$i" == "9" ]; then
        args="-p multinest.live_points=50"
    elif [ "$i" == "10" ]; then
        args="-p grid.nsample_dimension=10"
    elif [ "$i" == "11" ]; then
        args="-p grid.nsample_dimension=10"
    elif [ "$i" == "13" ]; then
        args="-p snake.threshold=3"
    elif [ "$i" == "14" ]; then
        args="-p zeus.samples=15 zeus.nsteps=5"
    elif [ "$i" == "17" ]; then
        args="-v cosmological_parameters.n_s=0.96 cosmological_parameters.omega_b=0.0468 cosmological_parameters.h0=0.6881 -p test.fatal_errors=T"
    elif [ "$i" == "19" ]; then
        args="-p grid.nsample_dimension=20"
    else
        args=
    fi
    if [ "$args" == "" ]; then
        args="-p runtime.verbosity=debug"
    else
        args="$args  runtime.verbosity=debug"
    fi
    echo Running demo $i:  cosmosis demos/demo${i}.ini $args
    time cosmosis demos/demo${i}.ini $args
    echo Done demo $i
    echo ""
    echo ""
    echo "Postprocessing demo $i" 
    if [ "$i" == "10" ]; then
        # Run the second part of demo 10
        export HALOFIT=mead2020
        echo "Running demo 10 (part 2):  cosmosis demos/demo${i}.ini $args"
        time cosmosis demos/demo${i}.ini $args

        # Postprocess demo10
        cosmosis-postprocess output/demo10_mead2020.txt output/demo10_takahashi.txt -o output/plots/demo10
    else
        cosmosis-postprocess demos/demo${i}.ini  -o output/plots/demo${i}
    fi
    echo ""
    echo ""
    echo ""
    echo ""
done
