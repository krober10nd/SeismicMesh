from MeshSizeFunction import MeshSizeFunction

fname = 'velocity_models/vel_z6.25m_x12.5m_exact.segy'

ef = MeshSizeFunction(
    bbox=(-12e3, 0, 0, 67e3),
    segy=fname,
    hmin=10)

fh = ef.build()

