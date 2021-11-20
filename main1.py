#!/usr/bin/env python3
# ---------------------
import numpy as np
from numba import njit
import matplotlib.pyplot as plt
from netCDF4 import Dataset

NTIME = 200

NX = 200
NY = 200

dx = 1.
dy = 1.
dt = 0.001

f = 1.e-4
g = 9.81


OUTPUTTIME = 1


def init_output():
    ncout = Dataset("out.nc", "w", "NETCDF4")
    ncout.createDimension("x", NX)
    ncout.createDimension("y", NY)
    ncout.createDimension("time", None)
    ncout.createVariable("eta", "float64", ("time", "x", "y"))
    ncout.createVariable("u", "float64", ("time", "x", "y"))
    ncout.createVariable("v", "float64", ("time", "x", "y"))
    ncout.close()


itime_write = 0


def write_output(eta, u, v):
    global itime_write
    with Dataset("out.nc", "a") as f:
        for variable, varname in [(eta, "eta"), (u, "u"), (v, "v")]:
            varout = f.variables[varname]
            varout[itime_write, :, :] = variable
    itime_write += 1


def advect(itime, eta_now, eta_new, u_now, v_now, u_new, v_new, H):
    # h
    def u_at_h(iu, ju): return (u_now[iu, ju] + u_now[iu+1, ju])/2.
    def v_at_h(iv, jv): return (v_now[iv, jv] + v_now[iv, jv+1])/2.
    def v_at_u(iv, jv): return (
        v_now[iv, jv] + v_now[iv, j+1] + v_now[iv-1, j+1] + v_now[iv-1, jv])/4.

    def h_at_u(iu, ju): return (eta_now[iu, ju] + eta_now[iu-1, ju])/2.
    def h_at_v(iv, jv): return (eta_now[iv, jv] + eta_now[iv, jv-1])/2.
    def u_at_v(iu, ju): return (
        u_now[iu, ju] + u_now[iu+1, ju] + u_now[iu, ju-1] + u_now[iu+1, ju-1])/2.

    for i in range(2, NX-2):
        for j in range(2, NY-2):

            dudx = (u_at_h(i+1, j) - u_at_h(i-1, j))/(2*dx)  # on h grid
            dvdy = (v_at_h(i, j+1) - v_at_h(i, j-1))/(2*dy)  # on h grid

            eta_new[i, j] = eta_now[i, j] - H[i, j]*(dudx + dvdy)*dt

    # u
    for i in range(2, NX-2):
        for j in range(2, NY-2):
            dhdx = (h_at_u(i+1, j) - h_at_u(i-1, j))/(2*dx)  # dh/dx at u
            u_new[i, j] = u_now[i, j] + (f*v_at_u(i, j) - g*dhdx)

    # v
    for i in range(2, NX-2):
        for j in range(2, NY-2):
            dhdy = (h_at_v(i, j+1) - h_at_v(i, j-1))/(2*dy)  # dh/dy at v
            v_new[i, j] = v_now[i, j] - (f*u_at_v(i, j) + g*dhdy)

    return eta_new, u_new, v_new


def main():

    H = np.ones((NX, NY))*100.

    # Initial conditions
    eta_now = np.zeros_like(H)
    eta_new = np.zeros_like(H)
    eta_now[NX//2-5: NX//2+5, NY//2-5:NY//2+5] = 0.3

    h_total = H + eta_now

    u_now = np.zeros((NX, NY))
    v_now = np.zeros((NX, NY))
    u_new = np.zeros((NX, NY))
    v_new = np.zeros((NX, NY))

    # Output
    init_output()

    # Integration
    itime = 0
    while (itime < NTIME):
        eta_new, u_new, v_new = advect(
            itime, eta_now, eta_new, u_now, v_now, u_new, v_new, H)

        if (itime % OUTPUTTIME == 0):
            write_output(eta_new, u_new, v_new)

        eta_now = eta_new
        u_now = u_new
        v_now = v_new

        itime += 1


if __name__ == "__main__":
    main()
