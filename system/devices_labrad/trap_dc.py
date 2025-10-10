from artiq.experiment import *
from LAX_exp.base import LAXDevice

import labrad
from os import environ
from EGGS_labrad.config.dc_config import dc_config

from time import sleep


class TrapDC(LAXDevice):
    """
    High-level API functions for configuring trap voltages via LabRAD interface.
    """
    name = "trap_dc"

    def prepare_device(self):
        # establish labrad client connections
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.trap_dc = self.cxn.dc_server

        # get relevant DC channels
        try:
            self.EAST_ENDCAP_CHANNEL =  dc_config.channeldict['E Endcap']['num']
            self.WEST_ENDCAP_CHANNEL =  dc_config.channeldict['W Endcap']['num']
            self.V_SHIM_CHANNEL =       dc_config.channeldict['V Shim']['num']
            self.H_SHIM_CHANNEL =       dc_config.channeldict['H Shim']['num']
            self.ARAMP_1_CHANNEL =      dc_config.channeldict['A-Ramp 1']['num']
            self.ARAMP_2_CHANNEL =      dc_config.channeldict['A-Ramp 2']['num']
        except Exception as e:
            print("Failed to get all trap DC channels: {}".format(e))

    @rpc
    def initialize_device(self) -> TNone:
        """
        Set up the AMO8 box to prevent device communication from being interrupted.
        """
        self.trap_dc.polling(False)
        self.trap_dc.alarm(False)
        self.trap_dc.serial_write('remote.w 1\r\n')
        self.trap_dc.serial_read('\n')


    """
    Generic functions
    """
    @rpc
    def voltage_fast(self, channel_num: TInt32, voltage: TFloat) -> TNone:
        """
        Generic fast voltage update.
        Warning: for speed, does not emit any signals to any clients.
        :param channel_num: the AMO8 DAC channel to update
        :param voltage: the voltage to set on the channel
        """
        self.trap_dc.voltage_fast(channel_num, voltage)

    @rpc
    def voltage_get(self, channel_num: TInt32) -> TFloat:
        """
        Generic voltage read.
        :param channel_num: the AMO8 DAC channel to read
        :return: the channel's voltage
        """
        return self.trap_dc.voltage(channel_num)


    """
    Toggle Channels
    """
    @rpc
    def east_endcap_toggle(self, status: TBool) -> TNone:
        """
        Toggle East Endcap on/off.
        :param status: channel output status.
        """
        self.trap_dc.toggle(self.EAST_ENDCAP_CHANNEL, status)

    @rpc
    def west_endcap_toggle(self, status: TBool) -> TNone:
        """
        Toggle West Endcap on/off.
        :param status: channel output status.
        """
        self.trap_dc.toggle(self.WEST_ENDCAP_CHANNEL, status)

    @rpc
    def v_shim_toggle(self, status: TBool) -> TNone:
        """
        Toggle V Shim on/off.
        :param status: channel output status.
        """
        self.trap_dc.toggle(self.V_SHIM_CHANNEL, status)

    @rpc
    def h_shim_toggle(self, status: TBool) -> TNone:
        """
        Toggle H Shim channel on/off.
        :param status: channel output status.
        """
        self.trap_dc.toggle(self.H_SHIM_CHANNEL, status)

    @rpc
    def aramp_toggle(self, status: TBool) -> TNone:
        """
        Toggle A-Ramp on/off.
        :param status: channel output status.
        """
        self.trap_dc.toggle(self.ARAMP_2_CHANNEL, status)


    """
    Set Voltages
    """
    @rpc
    def set_east_endcap_voltage(self, voltage: TFloat) -> TNone:
        """
        Set east endcap voltage.
        :param voltage: voltage to set for channel.
        """
        if not 0 <= voltage <= 500:
            raise ValueError("voltage must be between 0V and 500V")
        self.trap_dc.voltage_fast(self.EAST_ENDCAP_CHANNEL, voltage)

    @rpc
    def set_west_endcap_voltage(self, voltage: TFloat) -> TNone:
        """
        Set west endcap voltage.
        :param voltage: voltage to set for channel.
        """
        if not 0 <= voltage <= 500:
            raise ValueError("voltage must be between 0V and 500V")
        self.trap_dc.voltage_fast(self.WEST_ENDCAP_CHANNEL, voltage)

    @rpc
    def set_h_shim_voltage(self, voltage: TFloat) -> TNone:
        """
        Set H Shim voltage.
        :param voltage: voltage to set for channel.
        """
        if not 0 <= voltage <= 200:
            raise Exception("voltage must be between 0V and 200V")
        self.trap_dc.voltage_fast(self.H_SHIM_CHANNEL, voltage)

    @rpc
    def set_v_shim_voltage(self, voltage: TFloat) -> TNone:
        """
        Set V Shim voltage.
        :param voltage: voltage to set for channel.
        """
        if not 0 <= voltage <= 200:
            raise Exception("voltage must be between 0V and 200V")
        self.trap_dc.voltage_fast(self.V_SHIM_CHANNEL, voltage)

    @rpc
    def set_aramp_voltage(self, voltage: TFloat) -> TNone:
        """
        Set A-Ramp voltage.
        :param voltage: voltage to set for channel.
        """
        if not 0 <= voltage <= 100:
            raise Exception("voltage must be between 0V and 100V")
        self.trap_dc.voltage_fast(self.ARAMP_2_CHANNEL, voltage)


    """
    Ramp Voltages
    """
    @rpc
    def ramp_east_endcap(self, voltage: TFloat, rate: TFloat=100) -> TNone:
        """
        Ramp east endcap.
        :param voltage: ending voltage (in volts)
        :param rate: ramp rate (in volts/second).
        """
        if not 0 <= voltage <= 400:
            raise Exception("voltage must be between 0V and 400V")
        self.trap_dc.ramp(self.EAST_ENDCAP_CHANNEL, voltage, rate)

    @rpc
    def ramp_west_endcap(self, voltage: TFloat, rate: TFloat=100) -> TNone:
        """
        Ramp west endcap.
        :param voltage: ending voltage (in volts)
        :param rate: ramp rate (in volts/second).
        """
        if not 0 <= voltage <= 400:
            raise Exception("voltage must be between 0V and 300V")
        self.trap_dc.ramp(self.WEST_ENDCAP_CHANNEL, voltage, rate)

    @rpc
    def ramp_both_endcaps(self, voltages: TList(TFloat), rates: TList(TFloat)) -> TNone:
        """
        Simultaneously ramp both endcaps.
        :param voltages: list of voltages (east endcap voltage in index 0, west endcap voltage in index 1)
        :param rates: list of ramp rates (east endcap ramp rate in index 0, west endcap ramp rate in index 1)
        """
        # sanitize input
        if not 0 <= max(voltages) <= 400:
            raise Exception("voltage must be between 0V and 400V")
        if (len(voltages) != 2) or (len(rates) != 2):
            raise Exception("Arguments must be lists of [e_end_param, w_end_param].")

        try:
            self.trap_dc.ramp_multiple([self.EAST_ENDCAP_CHANNEL, self.WEST_ENDCAP_CHANNEL], voltages, rates)
        except Exception as e:
            print(repr(e))

        sleep(3)
