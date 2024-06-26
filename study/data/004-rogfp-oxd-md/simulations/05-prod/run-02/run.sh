#!/usr/bin/env bash

echo 'Starting 08_prod_npt'
date
pmemd.cuda -O -i $SLURM_SUBMIT_DIR/inputs/08_prod_npt.in -o $SLURM_SUBMIT_DIR/outputs/08_prod_npt.out -c $SLURM_SUBMIT_DIR/../../04-relax/run-02/outputs/07_relax_npt.rst -p $SLURM_SUBMIT_DIR/inputs/mol.prmtop -r $SLURM_SUBMIT_DIR/outputs/08_prod_npt.rst -x $SLURM_SUBMIT_DIR/outputs/08_prod_npt.nc -ref $SLURM_SUBMIT_DIR/inputs/mol.inpcrd -inf $SLURM_SUBMIT_DIR/outputs/08_prod_npt.mdinfo

echo 'Starting 09_prod_npt'
date
pmemd.cuda -O -i $SLURM_SUBMIT_DIR/inputs/09_prod_npt.in -o $SLURM_SUBMIT_DIR/outputs/09_prod_npt.out -c $SLURM_SUBMIT_DIR/outputs/08_prod_npt.rst -p $SLURM_SUBMIT_DIR/../../02-prep/mol.prmtop -r $SLURM_SUBMIT_DIR/outputs/09_prod_npt.rst -x $SLURM_SUBMIT_DIR/outputs/09_prod_npt.nc -ref $SLURM_SUBMIT_DIR/../../02-prep/mol.inpcrd -inf $SLURM_SUBMIT_DIR/outputs/09_prod_npt.mdinfo

echo 'Starting 10_prod_npt'
date
pmemd.cuda -O -i $SLURM_SUBMIT_DIR/inputs/10_prod_npt.in -o $SLURM_SUBMIT_DIR/outputs/10_prod_npt.out -c $SLURM_SUBMIT_DIR/outputs/09_prod_npt.rst -p $SLURM_SUBMIT_DIR/../../02-prep/mol.prmtop -r $SLURM_SUBMIT_DIR/outputs/10_prod_npt.rst -x $SLURM_SUBMIT_DIR/outputs/10_prod_npt.nc -ref $SLURM_SUBMIT_DIR/../../02-prep/mol.inpcrd -inf $SLURM_SUBMIT_DIR/outputs/10_prod_npt.mdinfo
