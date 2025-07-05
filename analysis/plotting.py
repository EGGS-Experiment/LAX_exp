"""
LAX.analysis.plotting

idk josh wrote this probably
"""
__all__ = [
    'plot_sidebands', 'plot_sb_probs', 'plot_phonons'
]

import numpy as np
import matplotlib.pyplot as plt
from LAX_exp.analysis import *


def plot_sidebands(rsb_freqs_MHz: np.array, bsb_freqs_MHz:np.array, probs_rsb: np.array,
                 probs_bsb: np.array, std_rsb: np.array, std_bsb: np.array, fit_rsb = None, fit_bsb = None):
    """
    Plot sideband spectra via qubit linescan.
    Args:
        rsb_freqs_MHz: rsb frequencies in MHz
        bsb_freqs_MHz: bsb frequencies in MHz
        probs_rsb: probabilitiies of rsb excitation
        probs_bsb: probabilitiies of bsb excitation
        std_rsb: standard deviation from averaging rsb excitations
        std_bsb: standard deviation from averaging bsb excitations
        fit_rsb: fit generated for rsb excitations over a frequncy span
        fit_bsb: fit generated for bsb excitations over a frequncy span
    """
    fig, axs = plt.subplots(2)

    axs[0].errorbar(rsb_freqs_MHz, probs_rsb, std_rsb, fmt='o', color="red")
    if fit_rsb is not None:
        axs[0].plot(np.linspace(np.min(rsb_freqs_MHz), np.max(rsb_freqs_MHz), len(fit_rsb)), fit_rsb, color="black")
    axs[0].set_ylim([0, 1])
    axs[0].set_xlim([np.min(rsb_freqs_MHz), np.max(rsb_freqs_MHz)])
    axs[0].ticklabel_format(useOffset=False)
    axs[0].set_ylabel('D State')

    axs[1].errorbar(bsb_freqs_MHz, probs_bsb, std_bsb, fmt='o', color="blue")
    if fit_bsb is not None:
        axs[1].plot(np.linspace(np.min(bsb_freqs_MHz), np.max(bsb_freqs_MHz), len(fit_bsb)), fit_bsb, color="black")
    axs[1].set_ylim([0, 1])
    axs[1].set_xlim([np.min(bsb_freqs_MHz), np.max(bsb_freqs_MHz)])
    axs[1].ticklabel_format(useOffset=False)
    axs[1].set_xlabel('EGGS Frequency (MHz)')
    axs[1].set_ylabel('D State')

def plot_sb_probs(scanning_freq_MHz: np.array,  probs_rsb: np.array,
                 probs_bsb: np.array, std_rsb: np.array, std_bsb: np.array):
    """
    Plot probability of sideband excitation for a frequency span
    Args:
        scanning_freq_MHz: frequency span
        probs_rsb: probabilitiies of rsb excitation
        probs_bsb: probabilitiies of bsb excitation
        std_rsb: standard deviation from averaging rsb excitations
        std_bsb: standard deviation from averaging bsb excitations
    """
    fig, axs = plt.subplots(2)

    axs[0].errorbar(scanning_freq_MHz, probs_rsb, std_rsb, fmt='o', color="red")
    axs[0].set_ylim([0, 1])
    axs[0].set_xlim([np.min(scanning_freq_MHz), np.max(scanning_freq_MHz)])
    axs[0].ticklabel_format(useOffset=False)
    axs[0].set_ylabel('D State')

    axs[1].errorbar(scanning_freq_MHz, probs_bsb, std_bsb, fmt='o', color="blue")
    axs[1].set_ylim([0, 1])
    axs[1].set_xlim([np.min(scanning_freq_MHz), np.max(scanning_freq_MHz)])
    axs[1].ticklabel_format(useOffset=False)
    axs[1].set_xlabel('EGGS Frequency (MHz)')
    axs[1].set_ylabel('D State')

def plot_phonons(scanning_freq_MHz: np.array, phonons: np.array, fit_data = None):
    """
    Plot phonons for a frequency span
    Args:
        scanning_freq_MHz: frequency span
        phonons: phonons
        fit_data: curve generated for phonons over a frequnecy span
    """
    fig, axs = plt.subplots(1)

    axs.errorbar(scanning_freq_MHz, phonons, fmt="o-", color="black")
    if fit_data is not None:
        axs.plot(np.linspace(np.min(scanning_freq_MHz), np.max(scanning_freq_MHz), len(fit_data)), fit_data, color="green")
    axs.set_ylim([0, 1])
    axs.set_xlim([np.min(scanning_freq_MHz), np.max(scanning_freq_MHz)])
    axs.ticklabel_format(useOffset=False)
    axs.set_xlabel('EGGS Frequency (MHz)')
    axs.set_ylabel('Phonons')
