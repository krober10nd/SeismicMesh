import numpy as np
import distmesh as dm
import matplotlib.pyplot as plt

from EdgeFunction import EdgeFunction

fname = 'vel_z6.25m_x12.5m_exact.segy'

ef = EdgeFunction(
            bbox=(-12e3,0,0,67e3),
            segy=fname,
            wl=5,freq=5,
            hmin=10)

