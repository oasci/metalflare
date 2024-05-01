# Jobs

Directory to put temporary jobs.

## Slurm

Here is a basic template for a Python job on slurm.

```bash
#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=6-00:00:00
#SBATCH --cluster=invest
#SBATCH --partition=jdurrant
#SBATCH --account=jdurrant
#SBATCH --job-name=001-sdf

# Load environment
module purge
source /ihome/jdurrant/amm503/miniconda3/bin/activate
conda activate metalflare-dev
```
