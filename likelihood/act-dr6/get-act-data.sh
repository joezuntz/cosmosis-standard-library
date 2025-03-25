#!/usr/bin/env bash


if [ -d data/ACTDR6MFLike/v1.0 ]
then
    echo ACT DR6 Lensing data already downloaded
else
    wget https://portal.nersc.gov/project/act/dr6_data/dr6_data.tar.gz
    tar -zxvf dr6_data.tar.gz
    mkdir -p data/ACTDR6MFLike/v1.0
    mv dr6_data.fits data/ACTDR6MFLike/v1.0
    rm dr6_data.tar.gz
    popd
fi