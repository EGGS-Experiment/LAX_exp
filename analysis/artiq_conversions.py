from numpy import pi, int32, int64

__all__ = []

# conversions
def ftw_to_frequency_mhz(ftw):
    return ftw * 4.656612880000000e-07/2

def ftw_to_frequency_khz(ftw):
    return ftw* 4.656612880000000e-07/2 * 1e3

def ftw_to_frequency(ftw):
    return ftw* 4.656612880000000e-07/2 * 1e6

def mu_to_us(time_mu):
    return time_mu/1e3

def mu_to_ms(time_mu):
    return time_mu/1e6

def mu_to_s(time_mu):
    return time_mu/1e9

def pow_to_turns(pow):
    return pow*1/(2**16. - 1.)

def pow_to_rads(pow):
    return 2*pi*pow*1/(2**16. - 1.)

def pow_to_turns(pow):
    return pow * 1 / (2 ** 16. - 1.)
# helper functions

def adc_mu_to_volt(data, gain=0, corrected_fs=True):
    """Convert ADC data in machine units to volts.

    :param data: 16-bit signed ADC word
    :param gain: PGIA gain setting (0: 1, ..., 3: 1000)
    :param corrected_fs: use corrected ADC FS reference.
        Should be ``True`` for Sampler revisions after v2.1. ``False`` for v2.1 and earlier.
    :return: Voltage in volts

    from https://m-labs.hk/artiq/manual-beta/_modules/artiq/coredevice/sampler.html
    """
    if gain == 0:
        volt_per_lsb = 20.48 / (1 << 16) if corrected_fs else 20. / (1 << 16)
    elif gain == 1:
        volt_per_lsb = 2.048 / (1 << 16) if corrected_fs else 2. / (1 << 16)
    elif gain == 2:
        volt_per_lsb = .2048 / (1 << 16) if corrected_fs else .2 / (1 << 16)
    elif gain == 3:
        volt_per_lsb = 0.02048 / (1 << 16) if corrected_fs else .02 / (1 << 16)
    else:
        raise ValueError("invalid gain")
    return data * volt_per_lsb

__all__.extend(['ftw_to_frequency_mhz', 'ftw_to_frequency_khz', 'mu_to_us',
                'mu_to_ms', 'pow_to_turns', 'pow_to_rads', 'pow_to_turns',
                'adc_mu_to_volt'])
