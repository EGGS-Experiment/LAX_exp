import matplotlib.pyplot as plt
import numpy as np
from artiq.experiment import *

from LAX_exp.analysis import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXExperiment
from LAX_exp.system.subsequences import (InitializeQubit, Readout, RescueIon,
                                         SidebandCoolContinuous, SidebandReadout)
from copy import deepcopy
from sipyco import pyon
from LAX_exp.extensions.physics_constants import *
from LAX_exp.extensions.conversions import *
import h5py


class Test(EnvExperiment):
    """
    Calibration: Data

    Examines the resonance for the EGGS feedthrough and calculates the appropriate amplitude scalings.
    """

    def build(self):
        self.scheduler = self.get_device("scheduler")

    def prepare(self):
        pass

    def run(self):
        pass

    def analyze(self):

        file_path = r'\\eric.physics.ucla.edu\\groups\\motion\\Data\\2024-05\\2024-05-02\\000058719-None.h5'

        with (h5py.File(file_path) as f):
            dataset = np.array(f['datasets']['results'])
            reps = f['arguments'].attrs['repetitions']

        sub_reps = 1

        # set up target voltages
        target_V_dipole = 0.015
        target_V_sideband = 8
        # values from experiment
        time_readout_s = list(dataset[:,5])[0] * 1e-6  # assume single readout time!!!
        time_readout_s = 127e-6
        secular_freq_hz = list(dataset[:,3])[0]  ### assume one secular frequncy used
        secular_freq_mhz = secular_freq_hz/ 1e6
        secular_freq_ang = 2 * np.pi * secular_freq_hz
        x0 = np.sqrt(hbar / (2 * mCa * secular_freq_ang))

        # grab data for each configuration
        sorting_col_num = 2

        rsb_qvsa_dataset = dataset[[bool(val) for val in dataset[:, 6]], :]
        bsb_qvsa_dataset = dataset[[bool(val) for val in dataset[:, 7]], :]
        squeezed_dataset = dataset[[bool(val) for val in dataset[:, 8]], :]

        # extract phonon number
        # tuple with (Ratios, ave_rsb, avd_bsb, std_rsb, std_bsb, scan_freqs)
        rsb_qvsa = extract_ratios(rsb_qvsa_dataset, sorting_col_num,
                                  1, 0, 75, sub_reps)


        rsb_qvsa_phonons = convert_ratios_to_coherent_phonons(rsb_qvsa[0])

        plt.figure(100)
        plt.plot(rsb_qvsa_phonons)

        bsb_qvsa = extract_ratios(bsb_qvsa_dataset, sorting_col_num,
                                  1, 0, 75, sub_reps)

        bsb_qvsa_phonons = convert_ratios_to_coherent_phonons(bsb_qvsa[0])

        squeezed = extract_ratios(squeezed_dataset, sorting_col_num,
                                  1, 0, reps, sub_reps)

        squeezed_phonons = convert_ratios_to_squeezed_phonons(squeezed[0])

        print(rsb_qvsa_phonons)
        print(bsb_qvsa_phonons)
        print(squeezed_phonons)

        assert len(bsb_qvsa_phonons) == len(rsb_qvsa_phonons) == len(squeezed_phonons), "Length Mismath"
        dipole_scaling_coeffs = dict()
        quadrupole_scaling_coeffs = dict()

        Vrs = np.zeros(len(rsb_qvsa_phonons))
        Vbs = np.zeros(len(rsb_qvsa_phonons))
        Vds = np.zeros(len(rsb_qvsa_phonons))
        carrier_freqs_mhz = np.zeros(len(rsb_qvsa_phonons))

        # try converting this from for loop to numpy array operations
        for idx, r_phonon in enumerate(rsb_qvsa_phonons):
            b_phonon = bsb_qvsa_phonons[idx]
            s_phonon = squeezed_phonons[idx]

            carrier_freq_mhz = rsb_qvsa[5][idx]
            carrier_freq_hz = rsb_qvsa[5][idx] * 1e6
            carrier_freq_ang = 2 * np.pi * carrier_freq_hz

            numerator_r = np.sqrt(
                (carrier_freq_ang - secular_freq_ang) * (carrier_freq_ang + secular_freq_ang)) * hbar * np.sqrt(
                np.arcsinh(np.sqrt(s_phonon))) * (r_phonon) ** (1 / 4) * r0 ** 2
            denominator_r = qe * np.sqrt(time_readout_s * secular_freq_ang) * (b_phonon) ** (1 / 4) * x0 ** 2
            Vr = numerator_r / denominator_r

            numerator_b = np.sqrt(
                (carrier_freq_ang - secular_freq_ang) * (carrier_freq_ang + secular_freq_ang)) * hbar * np.sqrt(
                np.arcsinh(np.sqrt(s_phonon))) * (b_phonon) ** (1 / 4) * r0 ** 2
            denominator_b = qe * np.sqrt(time_readout_s * secular_freq_ang) * (r_phonon) ** (1 / 4) * x0 ** 2
            Vb = numerator_b / denominator_b

            numerator_d = 2 * np.sqrt(
                (carrier_freq_ang - secular_freq_ang) * (carrier_freq_ang + secular_freq_ang) * (hbar ** 2) * np.sqrt(
                    r_phonon) * np.sqrt(b_phonon) * r0 ** 2)
            denominator_d = np.sqrt(qe ** 2 * time_readout_s * secular_freq_ang * x0 ** 2) * np.sqrt(
                np.arcsinh(np.sqrt(s_phonon)))
            Vd = numerator_d / denominator_d

            Vbs[idx] = Vb
            Vrs[idx] = Vr
            Vds[idx] = Vd
            carrier_freqs_mhz[idx] = carrier_freq_mhz


        f, (ax, ax2) = plt.subplots(1, 2, sharey=True, facecolor='w')
        ax.plot(carrier_freqs_mhz-secular_freq_mhz, Vrs, 'r', label="Red Sideband Voltages")
        ax2.plot(carrier_freqs_mhz+secular_freq_mhz, Vbs, 'b', label = "Blue Sideband Voltages")
        ax.set_xlim(np.min(carrier_freqs_mhz-secular_freq_mhz), np.max(carrier_freqs_mhz-secular_freq_mhz))
        ax2.set_xlim(np.min(carrier_freqs_mhz + secular_freq_mhz), np.max(carrier_freqs_mhz + secular_freq_mhz))
        ax.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        f.legend()
        plt.figure(2)
        plt.plot(carrier_freqs_mhz, Vds, 'black', label="Carrier Voltages")
        plt.legend()
        plt.show()


        #     if carrier_freq_mhz in dipole_scaling_coeffs.keys():
        #         dipole_scaling_coeffs[carrier_freq_mhz] = np.mean((target_V_dipole / Vd),
        #                                                           dipole_scaling_coeffs[carrier_freq_mhz])
        #     else:
        #         dipole_scaling_coeffs[carrier_freq_mhz] = target_V_dipole / Vd
        #
        #     if (carrier_freq_mhz + secular_freq_mhz) in quadrupole_scaling_coeffs.keys():
        #         quadrupole_scaling_coeffs[carrier_freq_mhz + secular_freq_mhz] = np.mean((target_V_sideband / Vb),
        #                                                                                  quadrupole_scaling_coeffs[
        #                                                                                      carrier_freq_mhz + secular_freq_mhz])
        #     else:
        #         quadrupole_scaling_coeffs[carrier_freq_mhz + secular_freq_mhz] = target_V_sideband / Vb
        #
        #     if (carrier_freq_mhz - secular_freq_mhz) in quadrupole_scaling_coeffs.keys():
        #         quadrupole_scaling_coeffs[carrier_freq_mhz - secular_freq_mhz] = np.mean((target_V_sideband / Vr),
        #                                                                                  quadrupole_scaling_coeffs[
        #                                                                                      carrier_freq_mhz - secular_freq_mhz])
        #         quadrupole_scaling_coeffs[carrier_freq_mhz - secular_freq_mhz] = target_V_sideband / Vr
        #     else:
        #         quadrupole_scaling_coeffs[carrier_freq_mhz - secular_freq_mhz] = target_V_sideband / Vr
        #
        # for key in quadrupole_scaling_coeffs.keys():
        #     if key in calibrations_eggs_scaling_coeffs_quadrupole.keys():
        #         calibrations_eggs_scaling_coeffs_quadrupole[key] *= quadrupole_scaling_coeffs[key]
        #     else:
        #         calibrations_eggs_scaling_coeffs_quadrupole[key] = quadrupole_scaling_coeffs
        #
        # for key in dipole_scaling_coeffs.keys():
        #     if key in calibrations_eggs_scaling_coeffs_dipole.keys():
        #         calibrations_eggs_scaling_coeffs_dipole[key] *= dipole_scaling_coeffs[key]
        #     else:
        #         calibrations_eggs_scaling_coeffs_dipole[key] = dipole_scaling_coeffs
