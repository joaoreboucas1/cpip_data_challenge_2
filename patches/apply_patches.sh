#!/bin/bash
set -xe

if [[ -z ${ROOTDIR} ]] then
    echo "ERROR: environment variable ROOTDIR not defined"
    exit 1
fi

echo "Applying patches to ${ROOTDIR}/projects/roman_real/likelihood/"

mv ${ROOTDIR}/projects/roman_real/likelihood/_cosmolike_prototype_base.py ${ROOTDIR}/projects/roman_real/likelihood/_cosmolike_prototype_base.py.bak
cp ${ROOTDIR}/projects/cpip_data_challenge_2/patches/_cosmolike_prototype_base.py ${ROOTDIR}/projects/roman_real/likelihood/

mv ${ROOTDIR}/projects/roman_real/likelihood/cosmic_shear.yaml ${ROOTDIR}/projects/roman_real/likelihood/cosmic_shear.yaml.bak
cp ${ROOTDIR}/projects/cpip_data_challenge_2/patches/cosmic_shear.yaml ${ROOTDIR}/projects/roman_real/likelihood/