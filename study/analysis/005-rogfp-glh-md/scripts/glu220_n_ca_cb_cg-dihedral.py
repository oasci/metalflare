study/analysis/005-rogfp-glh-md/scripts/cys202_n_ca_cb_sg-dihedral.py

    u = mda.Universe(topology_path, trajectory_paths)
    n_frames = len(u.trajectory)

    atoms = u.select_atoms("(resid 220 and name N CA CB CG)")

    atoms_npy_path = os.path.join(data_dir, "glu220_n_ca_cb_cg-dihedral.npy")
    atoms_dihedral_array = np.full((n_frames,), np.nan, dtype=np.float64)

    for i, ts in enumerate(u.trajectory):
        coords = atoms.positions
        dihedral_angle = mda.lib.distances.calc_dihedrals(*coords)
        atoms_dihedral_array[i] = dihedral_angle

    np.save(atoms_npy_path, atoms_dihedral_array)

    print(atoms_dihedral_array)


if __name__ == "__main__":
    main()
