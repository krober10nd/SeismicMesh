import SeismicMesh


def example_write():
    # Name of SEG-Y file containg velocity model.
    fname = "velocity_models/vel_z6.25m_x12.5m_exact.segy"
    vp = SeismicMesh.read_segy(fname)
    bbox = (-12e3, 0, 0, 67e3)

    # Construct mesh sizing object from velocity model
    ef = SeismicMesh.MeshSizeFunction(
        bbox=bbox,
        velocity_grid=vp,
        domain_ext=2e3,
        dt=0.001,
        freq=5,
        wl=5,
        hmax=1e3,
        hmin=50.0,
        grade=0.05,
    )

    ef = ef.build()

    # Write it to an hdf5 file for later use.
    ef.write_velocity_model("BP2004")

    # Save your mesh size function options.
    ef.save_mesh_size_function_options("BP2004_SizingFunction")


if __name__ == "__main__":

    example_write()
