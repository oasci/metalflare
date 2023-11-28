from metalflare.simulation.amber.run import AmberRunPrep


class TestAmberRunPrep:
    def test_amber_run_command_prep_01_min(self, amber_simulation_standard_context):
        simulation_context = AmberRunPrep.prepare_context(
            amber_simulation_standard_context
        )
        run_command = AmberRunPrep.get_stage_run_command(simulation_context.get())
        assert len(run_command) == 4
        assert run_command[0] == "echo 'Starting 01_min'"
        assert run_command[1] == "date"
        assert run_command[-1] == ""

        amber_command = run_command[2]
        assert "None" not in amber_command
        assert "mpirun -np 12 pmemd.MPI -O" == amber_command[:26]
        amber_command_split = amber_command[26:].split(" -")[1:]
        assert "i " == amber_command_split[0][:2]
        assert "metalflare/tests/01_min.in" in amber_command_split[0]
        assert "o " == amber_command_split[1][:2]
        assert "metalflare/tests/01_min.out" in amber_command_split[1]
        assert "c " == amber_command_split[2][:2]
        assert "metalflare/tests/mol.inpcrd" in amber_command_split[2]
        assert "p " == amber_command_split[3][:2]
        assert "metalflare/tests/mol.prmtop" in amber_command_split[3]
        assert "r " == amber_command_split[4][:2]
        assert "metalflare/tests/01_min.rst" in amber_command_split[4]
        assert "x " == amber_command_split[5][:2]
        assert "metalflare/tests/01_min.nc" in amber_command_split[5]
        assert "ref " == amber_command_split[6][:4]
        assert "metalflare/tests/mol.inpcrd" in amber_command_split[6]
        assert "inf " == amber_command_split[7][:4]
        assert "metalflare/tests/01_min.mdinfo" in amber_command_split[7]

    def test_amber_prepare_stage(self, amber_simulation_standard_context):
        simulation_context = AmberRunPrep.prepare_context(
            amber_simulation_standard_context
        )
        stage_input_lines, _ = AmberRunPrep.prepare_stage(
            simulation_context.get(), write=False
        )
        stage_input_lines_ref = [
            "01_min",
            "&cntrl",
            "    irest=1,",
            "    ntx=5,",
            "    ig=-1,",
            "    dt=0.002,",
            "    nstlim=500000,",
            "    nscm=500,",
            "    ntr=1,",
            "    restraint_wt=0.5,",
            "    restraintmask=\"!(:WAT) & (@C,CA,N,O,O5',P,O3',C3',C4',C5')\",",
            "    ntb=2,",
            "    ntf=2,",
            "    ntc=2,",
            "    cut=10.0,",
            "    ntt=3,",
            "    temp0=300.0,",
            "    gamma_ln=5.0,",
            "    ntp=1,",
            "    barostat=2,",
            "    pres0=1.01325,",
            "    mcbarint=100,",
            "    comp=44.6,",
            "    taup=1.0,",
            "    ntxo=2,",
            "    ntwr=5000,",
            "    ntpr=500,",
            "    ntwx=5000,",
            "    ioutfm=1,",
            "    iwrap=1,",
            "&end",
        ]
        assert stage_input_lines == stage_input_lines_ref
