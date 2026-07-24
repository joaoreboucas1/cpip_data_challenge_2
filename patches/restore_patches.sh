#!/bin/bash
set -xe

if [[ -z ${ROOTDIR} ]] then
    echo "ERROR: environment variable ROOTDIR not defined"
    exit 1
fi

echo "Restoring original files in ${ROOTDIR}/projects/roman_real/likelihood/"

mv ${ROOTDIR}/projects/roman_real/likelihood/_cosmolike_prototype_base.py.bak ${ROOTDIR}/projects/roman_real/likelihood/_cosmolike_prototype_base.py

mv ${ROOTDIR}/projects/roman_real/likelihood/cosmic_shear.yaml.bak ${ROOTDIR}/projects/roman_real/likelihood/cosmic_shear.yaml