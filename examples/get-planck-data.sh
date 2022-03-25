#!/usr/bin/env bash
set -x

pushd likelihood/planck2018
URL="http://pla.esac.esa.int/pla-sl/data-action?COSMOLOGY.COSMOLOGY_OID=151902"
FILENAME="COM_Likelihood_Data-baseline_R3.00.tar.gz"
wget -O ${FILENAME} ${URL}
tar -zxvf ${FILENAME}
popd
