
echo 'Starting 05_relax_nvt_r'
date
pmemd.cuda -O -i $SLURM_SUBMIT_DIR/inputs/05_relax_nvt_r.in -o $SLURM_SUBMIT_DIR/outputs/05_relax_nvt_r.out -c $SLURM_SUBMIT_DIR/../03-min/outputs/04_min.rst -p $SLURM_SUBMIT_DIR/inputs/mol.prmtop -r $SLURM_SUBMIT_DIR/outputs/05_relax_nvt_r.rst -x $SLURM_SUBMIT_DIR/outputs/05_relax_nvt_r.nc -ref $SLURM_SUBMIT_DIR/inputs/mol.inpcrd -inf $SLURM_SUBMIT_DIR/outputs/05_relax_nvt_r.mdinfo

echo 'Starting 06_relax_npt_r'
date
pmemd.cuda -O -i $SLURM_SUBMIT_DIR/inputs/06_relax_npt_r.in -o $SLURM_SUBMIT_DIR/outputs/06_relax_npt_r.out -c $SLURM_SUBMIT_DIR/outputs/05_relax_nvt_r.rst -p $SLURM_SUBMIT_DIR/inputs/mol.prmtop -r $SLURM_SUBMIT_DIR/outputs/06_relax_npt_r.rst -x $SLURM_SUBMIT_DIR/outputs/06_relax_npt_r.nc -ref $SLURM_SUBMIT_DIR/inputs/mol.inpcrd -inf $SLURM_SUBMIT_DIR/outputs/06_relax_npt_r.mdinfo

echo 'Starting 07_relax_npt'
date
pmemd.cuda -O -i $SLURM_SUBMIT_DIR/inputs/07_relax_npt.in -o $SLURM_SUBMIT_DIR/outputs/07_relax_npt.out -c $SLURM_SUBMIT_DIR/outputs/06_relax_npt_r.rst -p $SLURM_SUBMIT_DIR/inputs/mol.prmtop -r $SLURM_SUBMIT_DIR/outputs/07_relax_npt.rst -x $SLURM_SUBMIT_DIR/outputs/07_relax_npt.nc -ref $SLURM_SUBMIT_DIR/inputs/mol.inpcrd -inf $SLURM_SUBMIT_DIR/outputs/07_relax_npt.mdinfo
