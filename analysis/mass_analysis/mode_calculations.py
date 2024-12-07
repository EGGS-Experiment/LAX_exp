import numpy as np
from scipy.optimize import fsolve

"""Fundamental Constants"""
qe = 1.60217662e-19
e0 = 8.8541878188e-12

"""Units"""
Hz = 1
kHz = 1e3
MHz = 1e6
m = 1
mm = 1e-3
um = 1e-6
nm = 1e-9
amu = 1.6605390e-27

"""Trap and Ion Parameters"""
Omega_RF = 2 * np.pi * 19.03908 * MHz
r0 = 550 * um
mCa = 40 * amu
z0 = 2.2 * mm


def calc_a_params(wsec_list, Omega_RF):
    wx, wy, wz = np.flip(np.sort(wsec_list))

    ar = 2 * (wx ** 2 - wy ** 2) / Omega_RF ** 2
    az = (2 * wz ** 2) / Omega_RF ** 2
    return np.array([ar - az, -ar - az, az])


def calc_q_params(wsec_list, a_params_list, Omega_RF):
    wx, wy, wz = np.flip(np.sort(wsec_list))
    ax, ay, az = a_params_list
    ar = (ax - ay) / 2

    qx = np.sqrt(2 * ((4 * wx ** 2) / (Omega_RF ** 2) + (az - ar)))
    return np.array([qx, -qx, 0])


def calc_mass_indep_params(a_params_list, q_params_list):
    ax, ay, az = a_params_list
    qx, qy, qz = q_params_list
    qxp = (mCa / amu) * qx
    qyp = (mCa / amu) * qy
    qzp = (mCa / amu) * qz
    axp = (mCa / amu) * ax
    ayp = (mCa / amu) * ay
    azp = (mCa / amu) * az

    return np.array([qxp, qyp, qzp, axp, ayp, azp])


def calc_new_sec_freqs(new_mass, q_params_list, a_params_list, Omega_RF):
    qx, qy, qz = np.array(q_params_list) / new_mass
    ax, ay, az = np.array(a_params_list) / new_mass

    return np.array([1 / 2 * Omega_RF * np.sqrt(ax + 1 / 2 * qx ** 2), 1 / 2 * Omega_RF * np.sqrt(ay + 1 / 2 * qy ** 2),
                     1 / 2 * Omega_RF * np.sqrt(2 * az)])


def calc_equilibrium_pos(z, mass_arr, omega_z):
    eqn_list = mass_arr * omega_z ** 2 * z
    for i in range(len(z)):

        if z[:i].size != 0:
            eqn_list[i] += np.sum(1 / ((z[i] - z[:i]) ** 2))
        if z[i + 1:].size != 0:
            eqn_list[i] -= np.sum(1 / ((z[i] - z[i + 1:]) ** 2))

    return eqn_list


def calc_normal_modes(mass_arr, sec_freqs_arr, z_arr):
    axes = sec_freqs_arr.shape[1]

    num_ions = len(mass_arr)
    normal_modes = np.zeros((num_ions, axes))
    normal_vecs = np.zeros((num_ions,3,axes))

    B = np.array(np.zeros((num_ions, num_ions, 3)))
    B += sec_freqs_arr ** 2 * np.repeat(np.eye(num_ions), 3, axis=1).reshape((num_ions, num_ions, 3))

    for axis in range(axes):

        A = B[:, :, axis]

        for n in range(num_ions):
            for k in range(num_ions):

                if n == k and axis == 2:
                    A[n, k] += 2 / mass_arr[n] * np.sum(1 / (np.abs(z_arr[n] - np.delete(z_arr, n)) ** 3))

                elif n != k and axis == 2:
                    A[n, k] -= 2 / np.sqrt(mass_arr[n] * mass_arr[k]) * 1 / (np.abs(z_arr[n] - z_arr[k]) ** 3)

                elif n == k and axis != 2:
                    A[n, k] -= 1 / mass_arr[n] * np.sum(1 / (np.abs(z_arr[n] - np.delete(z_arr, n)) ** 3))

                elif n != k and axis != 2:
                    A[n, k] += 1 / np.sqrt(mass_arr[n] * mass_arr[k]) * 1 / (np.abs(z_arr[n] - z_arr[k]) ** 3)

        eigenvalues, eigenvectors = np.linalg.eig(A)

        normal_modes[:, axis] = np.sqrt(eigenvalues)
        # normal_vecs[:,:, axis] = eigenvectors

        for i in range(len(eigenvalues)):


            print(np.sqrt(eigenvalues[i])/(2 * np.pi * MHz), eigenvectors[:,i])

    return np.transpose(normal_modes)


if __name__ == "__main__":
    wsec_list = 2 * np.pi * np.array([1.302, 1.5894, 0.7024]) * MHz
    wsec_list = np.flip(np.sort(wsec_list))
    a_params_list = calc_a_params(wsec_list, Omega_RF)
    q_params_list = calc_q_params(wsec_list, a_params_list, Omega_RF)
    mass_indep_params_list = calc_mass_indep_params(a_params_list, q_params_list)
    new_mass = 36
    new_sec_freqs = calc_new_sec_freqs(new_mass, mass_indep_params_list[:3], mass_indep_params_list[3:], Omega_RF)

    wzs = np.array([wsec_list[-1], new_sec_freqs[-1]])

    z = np.zeros(len(wzs))
    for i in range(len(wzs)):
        z[i] = i * um

    mass_arr = np.array([mCa / amu, new_mass]) * amu
    sec_freqs_arr = np.array([wsec_list, new_sec_freqs])

    pos_arr = fsolve(calc_equilibrium_pos, np.flip(np.sort(z)), args=(np.array([40, new_mass]) * amu, wzs))
    normal_modes = calc_normal_modes(mass_arr, sec_freqs_arr, pos_arr)
    # print(normal_modes / (2 * np.pi * MHz))
