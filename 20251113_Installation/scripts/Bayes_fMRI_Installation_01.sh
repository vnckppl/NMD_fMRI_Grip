#!/usr/bin/env bash
# * Install prerequisites for compiling BayesfMRI in R

# ** For 'sf'
sudo apt-get install -y \
     libudunits2-dev \
     libgdal-dev \
     libgeos-dev \
     libproj-dev

# ** For 'excursions'
sudo apt-get install -y libgsl-dev
