from artiq.experiment import *
from LAX_exp.base import LAXDevice

import time
import labrad
from os import environ
from EGGS_labrad.config.dc_config import dc_config


class TrapDC(LAXDevice):
    """
    High-level API functions for configuring trap voltages
    """
    name = "trap_dc"

    # todo: make it get from dc_config somehow (preferably registry)
    EAST_ENDCAP_CHANNEL =   dc_config.channeldict['E Endcap']['num']
    WEST_ENDCAP_CHANNEL =   dc_config.channeldict['W Endcap']['num']
    V_SHIM_CHANNEL =        dc_config.channeldict['V Shim']['num']
    H_SHIM_CHANNEL =        dc_config.channeldict['H Shim']['num']
    ARAMP_1_CHANNEL =       dc_config.channeldict['A-Ramp 1']['num']
    ARAMP_2_CHANNEL =       dc_config.channeldict['A-Ramp 2']['num']

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'], port=7682, tls_mode='off', username='', password='lab')
        self.trap_dc = self.cxn.dc_server


    """Toggle Channels"""
    @rpc
    def east_endcap_toggle(self, status: TBool) -> TNone:
        """
        Toggle East Endcap on/off.
        """
        self.trap_dc.toggle(self.EAST_ENDCAP_CHANNEL, status)

    @rpc
    def west_endcap_toggle(self, status: TBool) -> TNone:
        """
        Toggle West Endcap on/off.
        """
        self.trap_dc.toggle(self.WEST_ENDCAP_CHANNEL, status)

    @rpc
    def v_shim_toggle(self, status: TBool) -> TNone:
        """
        Toggle V Shim on/off.
        """
        self.trap_dc.toggle(self.V_SHIM_CHANNEL, status)

    @rpc
    def h_shim_toggle(self, status: TBool) -> TNone:
        """
        Toggle H Shim channel on/off.
        """
        self.trap_dc.toggle(self.H_SHIM_CHANNEL, status)

    @rpc
    def aramp_toggle(self, status: TBool) -> TNone:
        """
        Toggle A-Ramp on/off.
        """
        self.trap_dc.toggle(self.ARAMP_2_CHANNEL, status)


    """Set Voltages"""

    @rpc
    def set_east_endcap_voltage(self, voltage: TFloat) -> TNone:
        """
        Set east endcap voltage.
        Args:
            voltage: voltage to set channel to
        Returns:
            None
        """
        if not 0 <= voltage <= 400:
            raise Exception("voltage must be between 0V and 400V")
        self.trap_dc.voltage_fast(self.EAST_ENDCAP_CHANNEL, voltage)

    @rpc
    def set_west_endcap_voltage(self, voltage: TFloat) -> TNone:
        """
        Set west endcap voltage.
        Args:
            voltage: voltage to set channel to
        Returns:
            None
        """
        if not 0 <= voltage <= 400:
            raise Exception("voltage must be between 0V and 400V")
        self.trap_dc.voltage_fast(self.WEST_ENDCAP_CHANNEL, voltage)

    @rpc
    def set_h_shim_voltage(self, voltage: TFloat) -> TNone:
        """
        Set H Shim voltage.
        Args:
            voltage: voltage to set channel to
        """
        if not 0 <= voltage <= 150:
            raise Exception("voltage must be between 0V and 150V")
        self.trap_dc.voltage_fast(self.H_SHIM_CHANNEL, voltage)

    @rpc
    def set_v_shim_voltage(self, voltage: TFloat) -> TNone:
        """
        Set V Shim voltage.
        Args:
            voltage: voltage to set channel to
        """
        if not 0 <= voltage <= 150:
            raise Exception("voltage must be between 0V and 150V")
        self.trap_dc.voltage_fast(self.V_SHIM_CHANNEL, voltage)

    @rpc
    def set_aramp_voltage(self, voltage: TFloat) -> TNone:
        """
        Set A-Ramp voltage.
        Args:
            voltage: voltage to set channel to
        """
        if not 0 <= voltage <= 100:
            raise Exception("voltage must be between 0V and 100V")
        self.trap_dc.voltage_fast(self.ARAMP_2_CHANNEL, voltage)


    """Ramp Voltages"""
    @rpc
    def ramp_east_endcap(self, voltage: TFloat, rate: TFloat=100) -> TNone:
        """
        Ramp east endcap.
        Args:
            voltage: ending voltage
            rate: ramp rate
        """
        if not 0 <= voltage <= 400:
            raise Exception("voltage must be between 0V and 400V")
        self.trap_dc.ramp(self.EAST_ENDCAP_CHANNEL, voltage, rate)

    @rpc
    def ramp_west_endcap(self, voltage: TFloat, rate: TFloat=100) -> TNone:
        """
        Ramp west endcap.
        Args:
            voltage: ending voltage
            rate: ramp rate
        """
        if not 0 <= voltage <= 400:
            raise Exception("voltage must be between 0V and 300V")
        self.trap_dc.ramp(self.WEST_ENDCAP_CHANNEL, voltage, rate)

    @rpc
    def ramp_both_endcaps(self, voltages: TList(TFloat), rates: TList(TFloat)) -> TNone:
        """
        Simultaneously ramp both endcaps.
        Args:
            voltages    (TList(TFloat)): list of voltages (east endcap voltage in index 0, west endcap voltage in index 1)
            rates       (TList(TFloat)): list of ramp rates (east endcap ramp rate in index 0, west endcap ramp rate in index 1)
        Returns:
            None
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

        time.sleep(3)
