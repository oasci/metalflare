echo 'Starting 01_min'
date
pmemd.cuda -O -i $SLURM_SUBMIT_DIR/01_min.in -o $SLURM_SUBMIT_DIR/01_min.out -c $SLURM_SUBMIT_DIR/mol.inpcrd -p $SLURM_SUBMIT_DIR/mol.prmtop -r $SLURM_SUBMIT_DIR/01_min.rst -x $SLURM_SUBMIT_DIR/01_min.nc -ref $SLURM_SUBMIT_DIR/mol.inpcrd -inf $SLURM_SUBMIT_DIR/01_min.mdinfo

echo 'Starting 02_min'
date
pmemd.cuda -O -i $SLURM_SUBMIT_DIR/02_min.in -o $SLURM_SUBMIT_DIR/02_min.out -c $SLURM_SUBMIT_DIR/01_min.rst -p $SLURM_SUBMIT_DIR/mol.prmtop -r $SLURM_SUBMIT_DIR/02_min.rst -x $SLURM_SUBMIT_DIR/02_min.nc -ref $SLURM_SUBMIT_DIR/mol.inpcrd -inf $SLURM_SUBMIT_DIR/02_min.mdinfo

echo 'Starting 03_min'
date
pmemd.cuda -O -i $SLURM_SUBMIT_DIR/03_min.in -o $SLURM_SUBMIT_DIR/03_min.out -c $SLURM_SUBMIT_DIR/02_min.rst -p $SLURM_SUBMIT_DIR/mol.prmtop -r $SLURM_SUBMIT_DIR/03_min.rst -x $SLURM_SUBMIT_DIR/03_min.nc -ref $SLURM_SUBMIT_DIR/mol.inpcrd -inf $SLURM_SUBMIT_DIR/03_min.mdinfo

echo 'Starting 04_min'
date
pmemd.cuda -O -i $SLURM_SUBMIT_DIR/04_min.in -o $SLURM_SUBMIT_DIR/04_min.out -c $SLURM_SUBMIT_DIR/03_min.rst -p $SLURM_SUBMIT_DIR/mol.prmtop -r $SLURM_SUBMIT_DIR/04_min.rst -x $SLURM_SUBMIT_DIR/04_min.nc -ref $SLURM_SUBMIT_DIR/mol.inpcrd -inf $SLURM_SUBMIT_DIR/04_min.mdinfo
