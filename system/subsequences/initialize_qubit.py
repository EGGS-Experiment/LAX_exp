from artiq.experiment import *
from LAX_exp.extensions import *
from LAX_exp.base import LAXSubsequence

class InitializeQubit(LAXSubsequence):
    """
    Subsequence: Initialize Qubit

    Initialize the ion in the S-1/2 mj=-1/2 state, then cool to the doppler limit.
    """
    name = 'initialize_qubit'
    kernel_invariants = {
        "time_spinpol_mu",
        "time_repump_qubit_mu",
        "time_doppler_cooling_mu"
    }

    def build_subsequence(self):
        self.setattr_device('probe')
        self.setattr_device('pump')
        self.setattr_device('repump_cooling')
        self.setattr_device('repump_qubit')
        self.setattr_device('pmt')
        self.setattr_device('aperture')

    def prepare_subsequence(self):
        self.time_spinpol_mu =          self.get_parameter('time_spinpol_us', group='timing', override=True,
                                                           conversion_function=seconds_to_mu, units=us)
        self.time_repump_qubit_mu =     self.get_parameter('time_repump_qubit_us', group='timing', override=True,
                                                           conversion_function=seconds_to_mu, units=us)
        self.time_doppler_cooling_mu =  self.get_parameter('time_doppler_cooling_us', group='timing', override=True,
                                                           conversion_function=seconds_to_mu, units=us)\

        self.aperture_open_time_mu = self.core.seconds_to_mu(5)
        self.aperture_close_time_mu = self.core.seconds_to_mu(1)
        self.doppler_cooling_counts_threshold = 0.5

    @kernel(flags={"fast-math"})
    def run(self, detect_collision: TBool = False) -> TNone:
        """
        Quench ion from D-5/2, doppler, then run spinpol.
        Typical stateprep sequence sans SBC.
        """
        # ensure 397nm spinpol beam is off before starting
        self.probe.off()

        # set cooling waveform
        self.pump.cooling()

        # enable cooling repump (866nm)
        # 2025/03/27: evidently we don't turn 866nm off - should we?
        self.repump_cooling.on()

        # repump pulse
        self.repump_qubit.on()
        delay_mu(self.time_repump_qubit_mu)
        # 2025/03/27: why bother turning off BEFORE doppler?
        # should leave on instead b/c 397 has 393 component, which can scatter into D-5/2
        # self.repump_qubit.off()

        # doppler cooling
        self.pump.on()
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.off()

        # spin polarization
        self.probe.on()
        delay_mu(self.time_spinpol_mu)
        self.probe.off()
        # 2025/03/27: ensure 854nm off
        self.repump_qubit.off()

    @kernel(flags={"fast-math"})
    def initialize_with_collison_check(self):
        """
        Quench ion from D-5/2, doppler, then run spinpol.
        Typical stateprep sequence sans SBC.
        """
        # ensure 397nm spinpol beam is off before starting
        self.probe.off()

        # set cooling waveform
        self.pump.cooling()

        # enable cooling repump (866nm)
        # 2025/03/27: evidently we don't turn 866nm off - should we?
        self.repump_cooling.on()

        # repump pulse
        self.repump_qubit.on()
        delay_mu(self.time_repump_qubit_mu)
        # 2025/03/27: why bother turning off BEFORE doppler?
        # should leave on instead b/c 397 has 393 component, which can scatter into D-5/2
        # self.repump_qubit.off()

        counts = self.doppler_cool_and_count()

        for attempt in range(5):
            if counts < self.doppler_cooling_counts_threshold:
                print('Ion chain decrystallized - opening aperture.')
                self.pump.readout()
                self.core.wait_until_mu(now_mu())
                self.aperture.pulse_aperture_open(5)
                self.core.break_realtime()
                delay_mu(self.aperture_close_time_mu)

                counts = self.doppler_cool_and_count()

                if counts >= self.doppler_cooling_counts_threshold:
                    print('Ion chain crystallized - resuming experiment')

            if counts >= self.doppler_cooling_counts_threshold:
                break
        if counts < self.doppler_cooling_counts_threshold:
            print("Could not recrystallize ion chain")
            self.cancel_all_experiments()
            raise TerminationRequested

        # spin polarization
        self.probe.on()
        delay_mu(self.time_spinpol_mu)
        self.probe.off()
        # 2025/03/27: ensure 854nm off
        self.repump_qubit.off()

    @kernel(flags={"fast-math"})
    def doppler_cool_and_count(self) -> TInt32:
        # doppler cooling
        self.pump.cooling()
        self.pump.on()
        self.repump_cooling.on()
        self.pmt.count(self.time_doppler_cooling_mu)
        delay_mu(self.time_doppler_cooling_mu)
        self.pump.off()
        counts = self.pmt.fetch_count()

        return counts

    @kernel(flags={"fast-math"})
    def quench(self) -> TNone:
        """
        Quick quench function to simplify quenching for other exps.
        Note: 866nm is also on.
        """
        # set target profile
        self.pump.readout()

        # run quench
        # note: no 8ns delays here b/c done judiciously
        self.repump_qubit.on()
        self.repump_cooling.on()
        delay_mu(self.time_repump_qubit_mu)
        self.repump_cooling.off()
        self.repump_qubit.off()

    @kernel(flags={"fast-math"})
    def spin_polarize(self) -> TNone:
        """
        Quick spinpol function to simplify spinpol for other exps.
        """
        # set target profile
        self.pump.readout()

        # run spinpol
        # note: yes 8ns delays here b/c spinpol & 866 have no ext sw
        self.probe.on()
        delay_mu(8)
        self.repump_cooling.on()
        delay_mu(self.time_spinpol_mu)
        self.probe.off()
        delay_mu(8)
        self.repump_cooling.off()
        delay_mu(8)

    @kernel(flags={"fast-math"})
    def slack_rescue(self) -> TNone:
        """
        Leave relevant beams on in rescue mode.
        This function is intended to be called at the end of an experimental shot
            (or whenever we have long, nondeterministic delays), such that any
            extra time/slack can be used to rescue/cool the ion.
        """
        # set rescue profile
        self.pump.rescue()
        # turn on rescuing beams
        self.pump.on()
        self.repump_qubit.on()
        self.repump_cooling.on()
        delay_mu(8)

