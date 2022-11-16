#!/usr/bin/env bash

set -e

mkdir output

# Check the downloader works
./examples/get-planck-data.sh


# Run some examples
cosmosis examples/planck.ini | tee output/planck.log
grep 'Likelihood =  -6288.61' output/planck.log


cosmosis examples/planck_class.ini | tee output/class.log
# This may need updating as we modify the class interface
# The settings are not optimized
# grep 'Likelihood =  -5968.93' output/class.log
# This seems to give different results on different systems


cosmosis examples/wmap.ini | tee output/wmap.log

cosmosis examples/pantheon.ini -p emcee.samples=20
cosmosis-postprocess examples/pantheon.ini -o output/pantheon
test -f output/pantheon/cosmological_parameters--omega_m.png

cosmosis examples/des-y1.ini | tee output/des-y1.log
grep 'Likelihood =  5237.37' output/des-y1.log

cosmosis examples/des-y1.ini -p 2pt_shear.file=./shear/cl_to_corr/cl_to_corr.py 2pt_shear.corr_type=xi | tee output/des-y1.log
grep 'Likelihood =  5237.3' output/des-y1.log

cosmosis examples/des-y3.ini -p pk_to_cl_gg.save_kernels=T pk_to_cl.save_kernels=T | tee output/des-y3.log
grep 'Likelihood =  6043.23' output/des-y3.log

cosmosis examples/des-y3-mira-titan.ini | tee output/des-y3-mt.log
grep 'Likelihood =  6048.0' output/des-y3-mt.log

cosmosis examples/des-y3.ini -v halo_model_parameters.logT_AGN=8.2 -p camb.halofit_version=mead2020_feedback | tee output/des-y3.log
grep 'Likelihood =  6049.94' output/des-y3.log

cosmosis examples/euclid-emulator.ini
          test -f output/euclid-emulator/matter_power_nl/p_k.txt

cosmosis examples/w_model.ini






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
        args="-v cosmological_parameters.n_s=0.96 cosmological_parameters.omega_b=0.0468 cosmological_parameters.h0=0.6881"
    elif [ "$i" == "18" ]; then
        echo "Skipping demo 18"
        continue
    elif [ "$i" == "19" ]; then
        args="-p grid.nsample_dimension=20"
    else
        args=
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
