import time
from artiq.experiment import *

from LAX_exp.base import LAXDevice

from os import environ
import labrad
from EGGS_labrad.config.dc_config import dc_config


class TrapDC(LAXDevice):
    """
    High-level api functions for configuring trap voltages
    """

    name = "trap_dc"

    EAST_ENDCAP_CHANNEL = 28
    WEST_ENDCAP_CHANNEL = 27
    V_SHIM_CHANNEL = 20
    H_SHIM_CHANNEL = 25
    ARAMP_1_CHANNEL = 24
    ARAMP_2_CHANNEL = 23

    def prepare_device(self):
        self.cxn = labrad.connect(environ['LABRADHOST'],
                                  port=7682, tls_mode='off',
                                  username='', password='lab')
        self.trap_dc = self.cxn.dc_server

    """Toggle Channels ON"""

    @rpc
    def east_endcap_on(self) -> TNone:
        """
        Turn on east endcap
        """
        self.trap_dc.toggle(self.EAST_ENDCAP_CHANNEL, 1)

    @rpc
    def west_endcap_on(self) -> TNone:
        """
        Turn on west endcap
        """
        self.trap_dc.toggle(self.WEST_ENDCAP_CHANNEL, 1)

    @rpc
    def v_shim_on(self) -> TNone:
        """
        Turn on v shim channel
        """
        self.trap_dc.toggle(self.V_SHIM_CHANNEL, 1)

    @rpc
    def h_shim_on(self) -> TNone:
        """
        Turn on h shim channel
        """
        self.trap_dc.toggle(self.H_SHIM_CHANNEL, 1)

    @rpc
    def aramp2_on(self) -> TNone:
        """
        Turn on a ramp2 channel
        """
        self.trap_dc.toggle(self.ARAMP_2_CHANNEL, 1)

    """Toggle Channels OFF"""

    @rpc
    def east_endcap_off(self) -> TNone:
        """
        Turn off east endcap
        """
        self.trap_dc.toggle(self.EAST_ENDCAP_CHANNEL, 0)

    @rpc
    def west_endcap_off(self) -> TNone:
        """
        Turn off west endcap
        """
        self.trap_dc.toggle(self.WEST_ENDCAP_CHANNEL, 0)

    @rpc
    def v_shim_off(self) -> TNone:
        """
        Turn off v shim channel
        """
        self.trap_dc.toggle(self.V_SHIM_CHANNEL, 0)

    @rpc
    def h_shim_off(self) -> TNone:
        """
        Turn off h shim channel
        """
        self.trap_dc.toggle(self.H_SHIM_CHANNEL, 0)

    @rpc
    def aramp2_off(self) -> TNone:
        """
        Turn off a ramp2 channel
        """
        self.trap_dc.toggle(self.ARAMP_2_CHANNEL, 0)

    """Set Voltages"""

    @rpc
    def set_east_endcap_voltage(self, voltage):
        """
        Set east endcap voltage

        Args:
            voltage: voltage to set channel to

        Returns:
            None
        """
        if not 0 <= voltage <= 300:
            raise Exception("voltage must be between 0V and 300V")
        self.trap_dc.voltage_fast(self.EAST_ENDCAP_CHANNEL, voltage)

    @rpc
    def set_west_endcap_voltage(self, voltage):
        """
        Set west endcap voltage
        """
        if not 0 <= voltage <= 300:
            raise Exception("voltage must be between 0V and 300V")
        self.trap_dc.voltage_fast(self.WEST_ENDCAP_CHANNEL, voltage)

    @rpc
    def set_h_shim_voltage(self, voltage):
        """
        Set h shim voltage

        Args:
            voltage: voltage to set channel to

        Returns:
            None
        """
        if not 0 <= voltage <= 150:
            raise Exception("voltage must be between 0V and 150V")
        self.trap_dc.voltage_fast(self.H_SHIM_CHANNEL, voltage)

    @rpc
    def set_v_shim_voltage(self, voltage):
        """
        Set v shim voltage

        Args:
            voltage: voltage to set channel to

        Returns:
            None
        """
        if not 0 <= voltage <= 150:
            raise Exception("voltage must be between 0V and 150V")
        self.trap_dc.voltage_fast(self.V_SHIM_CHANNEL, voltage)

    @rpc
    def set_aramp2_voltage(self, voltage):
        """
        Set a ramp2 voltage

        Args:
            voltage: voltage to set channel to

        Returns:
            None
        """

        if not 0 <= voltage <= 100:
            raise Exception("voltage must be between 0V and 100V")
        self.trap_dc.voltage_fast(self.ARAMP_2_CHANNEL, voltage)

    """Ramp Voltages"""

    @rpc
    def ramp_east_endcap(self, voltage, rate=100):
        """
        Ramp east endcap

        Args:
            voltage: ending voltage
            rate: ramp rate

        Returns:
            None
        """

        if not 0 <= voltage <= 300:
            raise Exception("voltage must be between 0V and 300V")
        self.trap_dc.ramp(self.EAST_ENDCAP_CHANNEL, voltage, rate)

    @rpc
    def ramp_west_endcap(self, voltage, rate=100):
        """
        Ramp west endcap

        Args:
            voltage: ending voltage
            rate: ramp rate

        Returns:
            None
        """

        if not 0 <= voltage <= 300:
            raise Exception("voltage must be between 0V and 300V")
        self.trap_dc.ramp(self.WEST_ENDCAP_CHANNEL, voltage, rate)

    @rpc
    def ramp_both_endcaps(self, voltages, rates):

        """
        Simultaneously ramp both endcaps

        Args:
            voltages: list of ending voltages (east endcap voltage in index 0, west endcap voltage in index 1)
            rates: list of ramp rates (east endcap ramp rate in index 0, west endcap ramp rate in index 1)

        Returns:
            None
        """

        if not 0 <= max(voltages) <= 300:
            raise Exception("voltage must be between 0V and 300V")

        try:
            self.trap_dc.ramp_multiple([self.EAST_ENDCAP_CHANNEL, self.WEST_ENDCAP_CHANNEL], voltages, rates)
        except Exception as e:
            print(repr(e))

        time.sleep(3)

    """CLEAR ALL VOLTAGES"""

    @rpc
    def clear_all_voltages(self):
        """
        Clear all voltages
        """

        self.trap_dc.clear()
