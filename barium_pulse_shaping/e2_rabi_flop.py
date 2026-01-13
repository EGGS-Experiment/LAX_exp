from artiq.experiment import *
from artiq.coredevice import ad9910

from numpy import int64

import barium.lib.scripts.artiq_sequences.subsequences.doppler_cooling_2 as doppler_cooling
import barium.lib.scripts.artiq_sequences.subsequences.e2_laser_sub as e2laser 
import barium.lib.scripts.artiq_sequences.subsequences.state_detection as state_detection
import barium.lib.scripts.artiq_sequences.subsequences.deshelving_sub as deshelving
import barium.lib.scripts.artiq_sequences.subsequences.state_prep as state_prep
import barium.lib.scripts.artiq_sequences.subsequences.sideband_cooling as sideband_cooling
import barium.lib.scripts.artiq_sequences.subsequences.raman_sideband_cooling as raman_sideband_cooling
import barium.lib.scripts.artiq_sequences.subsequences.state_prep_138 as state_prep_138
import barium.lib.scripts.artiq_sequences.subsequences.microwaves as microwaves

### modified ####
from barium.lib.scripts.artiq_sequences.subsequences.pulse_shaping import dds_pulse_shape
### modified ####


class e2_laser_sweep(EnvExperiment): 
    def build(self):
        self.setattr_device("core")
        self.dc = doppler_cooling.DopplerCooling(self)
        self.ss = e2laser.E2Laser(self)
        self.sd = state_detection.StateDetection(self)
        self.ds = deshelving.Deshelving(self)
        self.sp = state_prep.StatePrep(self)
        self.sbc = sideband_cooling.SidebandCooling(self)
        self.sbc_r = raman_sideband_cooling.RamanSidebandCooling(self)
        self.sp_138 = state_prep_138.StatePrep138(self)
        self.microwaves = microwaves.Microwaves(self)

        ### modified ####
        # prepare values etc for pulseshape subsequence
        self.dds_to_pulse_shape = self.get_device("urukul1_ch1")
        self.dds_freq_to_set_in_ftw = self.dds_to_pulse_shape.frequency_to_ftw(125.3 * MHz)
        self.pulseshape_profile = 3

        # instantiate pulseshape subsequence
        self.pulseshape = dds_pulse_shape.DDSPulseShape(
            self,
            dds_target=self.dds_to_pulse_shape,

            # the DDS profile to be used (must be in [0,7])
            # note: this profile should NOT be used by other sequences etc
            #   (unless you really have to)
            ram_profile=self.pulseshape_profile,

            # the RAM start address - the AD9910 has [0, 1023] possible addresses
            ram_addr_start=0,

            # number of samples to use for the waveform
            # note: num_samples + ram_addr_start must be less than 1023
            #   (i.e. the waveform must fit in the available addresses)
            num_samples=500,

            # max amplitude the waveform should reach (must be in (0, 100), exclusive)
            # note: don't use 0 or 100, since there's some silly behavior in those cases
            # note: this can't be adjusted later on - if you really need to change the amplitude,
            #   you will need to use the attenuator
            ampl_max_pct=50,

            # the pulse shape we want to use - available ones are in pulse_shaping.PulseShapes.py
            pulse_shape="blackman",
        )
        ### modified ####

        self.dc_counts = []
        self.sd_counts = []
        self.col_time = 100


    ### modified ####
    @kernel(flags={"fast-math"})
    def run_single_pulse(self):
        """
        Sample operation for e.g. a single pulse (e.g. rabi flopping, linescan).
        """
        self.core.reset()

        ### modified ####
        # initialize hardware for pulseshape sequence
        self.pulseshape.sequence_initialize()
        ### modified ####

        for i in range(self.cycles):
            self.core.break_realtime()

            ### modified ####
            # this is the time we're going to use for our pulse
            time_to_scan_desired_mu = int64(5000) # 5us

            # configure the DDS we're going to pulse shape with (need to set frequency and attenuation)
            # note: we use set_ftw instead of set or set_mu since, when in RAM mode, the DDS ignores
            #   the frequency selected by set/set_mu, using instead the value in set_ftw.
            # todo: note about phase accumulator when not in RAM mode - will accumulate phase
            # note: this also doesn't HAVE to be set every loop - if set_ftw doesn't have to be changed every
            #   cycle, and it not changed by anyone else, we could instead set it immediately
            #   after the sequence has initialized (i.e. after self.pulseshape.initialize())
            self.dds_to_pulse_shape.set_ftw(self.dds_freq_to_set_in_ftw)
            self.dds.set_att(10. * dB)

            # note: if we want the phase to be continuous/coherent, we need to set the frequency
            #   in the single-tone profile to be the same as that in the FTW register
            #   this is b/c whenever the AD9910 is not in RAM mode, it uses the single-tone profiles, and will
            #   thus accumulate phase BETWEEN pulses according to the frequency in the single-tone profile.
            self.dds.set_mu(self.dds_freq_to_set_in_ftw, asf=0x2000,
                            phase_mode=ad9910.PHASE_MODE_CONTINUOUS,
                            profile=self.pulseshape_profile)

            # we want to calculate the ACTUAL pulse time, which can be different from the desired pulse
            #   time b/c the DDS has limited time step size
            # we also need to get the RAM timestep size, which we'll need for pulseshape.run()
            time_ram_step, time_to_scan_actual_mu = self.pulseshape.calculate_time(time_to_scan_desired_mu)
            ### modified ####

            """
            INITIALIZE
            """
            self.dc.run(self.dc_time, self.dc_cool_power,self.dc_repump_power) # doppler cool
            self.sp.run(self.sp_time, self.sp_sp_power,self.sp_repump_power,self.sp_use_133) # state prep
            if self.MM_dur > 1.0: # microwave stuff
                self.microwaves.run(self.MM_dur, self.MM_Freq, self.MM_Amp, self.MM_Att, 1.0, 1, 0, self.phases)


            """
            SBC?
            """
            if self.use_sbc == 1:
                self.sbc.run(self.sbc_dur, self.sbc_deshelve_time, self.sbc_sbc_cycles,
                             self.sbc_sp_cycles, self.sbc_sp_time, self.sbc_deshelve_power)
            if self.use_raman_sbc_r == 1:
                self.sbc_r.run(self.sbc_r_dur, self.sbc_r_rf_freq_scan, self.sbc_r_rf_amp_scan,
                               self.sbc_r_rf_freq_static, self.sbc_r_rf_amp_static, self.sbc_r_sp_power, self.sbc_r_repump_power,
                               self.sbc_r_pulsed, self.sbc_r_pulse_dur, self.sbc_r_op_dur, self.sbc_r_cycles)

            """
            FINAL RUN? idk
            """
            # state prep - 1762
            if self.sp_use_1762 == 1: self.sp_138.run(self.sp_shelve_dur, self.sp_deshelve_dur, self.sp_cycles)
            if self.wait_time > 1.0: delay(self.wait_time*us)

            ### modified ####
            # i'm going to assume that this is the spectroscopy function (i.e. the one that runs e.g. rabiflopping
            #   or the linescan)
            # self.ss.run(self.pulse_duration, self.sh_deshelve_power)

            # for a single pulse:
            self.pulseshape.run(time_ram_step, time_to_scan_desired_mu, phas_pow=0x0)
            ### modified ####


            # state detect
            self.sd.run(self.sd_time, self.sd_cool_power,self.sd_repump_power, False, self.sd_cooling_freq,
                        self.sd_detect_freq,self.sd_cooling_amp,self.sd_detect_amp,self.sd_change_freq)

            # deshelve/reset
            self.ds.run(self.deshelve_time, self.deshelve_power, self.deshelve_cool_power, self.deshelve_repump_power)

            # retrieve counts
            self.dc_counts[i] = self.dc.fetch_count()
            self.sd_counts[i] = self.sd.fetch_count()

        ### modified ####
        # clean up pulseshape sequence to leave hardware in a safe state
        self.pulseshape.sequence_cleanup()
        ### modified ####

    ### modified ####
    @kernel(flags={"fast-math"})
    def run_pulse_train(self):
        """
        Sample operation for a phase-coherent pulse train.
        """
        self.core.reset()

        ### modified ####
        # initialize hardware for pulseshape sequence
        self.pulseshape.sequence_initialize()
        ### modified ####

        for i in range(self.cycles):
            self.core.break_realtime()

            ### modified ####
            # this is the time we're going to use for our pulses
            time_pulse_train_single_shot_mu = int64(5000)  # 5us

            # configure the DDS we're going to pulse shape with (need to set frequency and attenuation)
            # note: we use set_ftw instead of set or set_mu since, when in RAM mode, the DDS ignores
            #   the frequency selected by set/set_mu, using instead the value in set_ftw.
            # todo: note about phase accumulator when not in RAM mode - will accumulate phase
            # note: this also doesn't HAVE to be set every loop - if set_ftw doesn't have to be changed every
            #   cycle, and it not changed by anyone else, we could instead set it immediately
            #   after the sequence has initialized (i.e. after self.pulseshape.initialize())
            self.dds_to_pulse_shape.set_ftw(self.dds_freq_to_set_in_ftw)
            self.dds.set_att(10. * dB)
            ### modified ####

            """
            INITIALIZE
            """
            self.dc.run(self.dc_time, self.dc_cool_power, self.dc_repump_power)  # doppler cool
            self.sp.run(self.sp_time, self.sp_sp_power, self.sp_repump_power, self.sp_use_133)  # state prep
            if self.MM_dur > 1.0:  # microwave stuff
                self.microwaves.run(self.MM_dur, self.MM_Freq, self.MM_Amp, self.MM_Att, 1.0, 1, 0, self.phases)

            """
            SBC?
            """
            if self.use_sbc == 1:
                self.sbc.run(self.sbc_dur, self.sbc_deshelve_time, self.sbc_sbc_cycles,
                             self.sbc_sp_cycles, self.sbc_sp_time, self.sbc_deshelve_power)
            if self.use_raman_sbc_r == 1:
                self.sbc_r.run(self.sbc_r_dur, self.sbc_r_rf_freq_scan, self.sbc_r_rf_amp_scan,
                               self.sbc_r_rf_freq_static, self.sbc_r_rf_amp_static, self.sbc_r_sp_power,
                               self.sbc_r_repump_power,
                               self.sbc_r_pulsed, self.sbc_r_pulse_dur, self.sbc_r_op_dur, self.sbc_r_cycles)

            """
            FINAL RUN? idk
            """
            # state prep - 1762
            if self.sp_use_1762 == 1: self.sp_138.run(self.sp_shelve_dur, self.sp_deshelve_dur, self.sp_cycles)
            if self.wait_time > 1.0: delay(self.wait_time * us)

            ### modified ####
            # if possible, we want to configure as close to the pulses as we can
            # note: this will cost ~10-ish us
            time_to_scan_actual_mu = self.pulseshape.configure_train(time_pulse_train_single_shot_mu)

            # run a bunch of pulses with different phases (for e.g. dynamical decoupling sequences)
            for phas_pow in (0x0, 0x4000, 0x8000, 0x4000, 0x1234, 0x4321, 0x888):
                self.pulseshape.configure_train(phas_pow)
            ### modified ####

            # state detect
            self.sd.run(self.sd_time, self.sd_cool_power, self.sd_repump_power, False, self.sd_cooling_freq,
                        self.sd_detect_freq, self.sd_cooling_amp, self.sd_detect_amp, self.sd_change_freq)

            # deshelve/reset
            self.ds.run(self.deshelve_time, self.deshelve_power, self.deshelve_cool_power, self.deshelve_repump_power)

            # retrieve counts
            self.dc_counts[i] = self.dc.fetch_count()
            self.sd_counts[i] = self.sd.fetch_count()

        ### modified ####
        # clean up pulseshape sequence to leave hardware in a safe state
        self.pulseshape.sequence_cleanup()
        ### modified ####
        
    def get_counts(self):
        self.total_counts[0] = self.dc_counts
        self.total_counts[1] = self.sd_counts
        return self.total_counts

    def set_vals(self, keys, vals):
        self.p = dict(zip(keys, vals))
        self.cycles = int(self.p['E2RabiFlop.Sequences_Per_Point'])
        self.wait_time = int(self.p['E2RabiFlop.Wait_Time'])
        self.dc_time = float(self.p['DopplerCooling.duration'])
        self.dc_cool_power = float(self.p['DopplerCooling.CoolingPower'])
        self.dc_repump_power = float(self.p['DopplerCooling.RepumpPower'])
        
        self.sp_time = float(self.p['StatePrep.duration'])
        self.sp_sp_power = float(self.p['StatePrep.StatePrepPower'])
        self.sp_repump_power = float(self.p['StatePrep.RepumpPower'])
        self.sp_use_133 = float(self.p['StatePrep.Use_133'])

        self.sp_use_1762 = float(self.p['StatePrep138.Use_1762'])
        self.sp_shelve_dur = float(self.p['StatePrep138.Shelve_Duration'])
        self.sp_deshelve_dur = float(self.p['StatePrep138.Deshelve_Duration'])
        self.sp_cycles = float(self.p['StatePrep138.Prep_Cycles'])
        
        self.sd_time = float(self.p['StateDetection.duration'])
        self.sd_cool_power = float(self.p['StateDetection.CoolingPower'])
        self.sd_repump_power = float(self.p['StateDetection.RepumpPower'])
        self.sd_cooling_freq = float(self.p['StateDetection.CoolingFrequency'])
        self.sd_cooling_amp = float(self.p['StateDetection.CoolingAmplitude'])
        self.sd_detect_freq = float(self.p['StateDetection.DetectionFrequency'])
        self.sd_detect_amp = float(self.p['StateDetection.DetectionAmplitude'])
        self.sd_change_freq = float(self.p['StateDetection.ChangeFrequency'])
        
        self.deshelve_time = float(self.p['DeshelvingSub.duration'])
        self.deshelve_power = float(self.p['DeshelvingSub.DeshelvingPower'])
        self.deshelve_cool_power = float(self.p['DeshelvingSub.CoolingPower'])
        self.deshelve_repump_power = float(self.p['DeshelvingSub.RepumpPower'])
        self.pulse_duration = float(self.p['E2LaserSub.Pulse_Duration'])
        self.sh_deshelve_power = float(self.p['E2LaserSub.DeshelvingPower'])
        self.sbc_dur = float(self.p['SidebandCooling.duration'])
        self.sbc_sp_cycles = float(self.p['SidebandCooling.StatePrepCycles'])
        self.sbc_sp_time = float(self.p['SidebandCooling.StatePrepTime'])
        self.sbc_deshelve_power = float(self.p['SidebandCooling.DeshelvingPower'])
        self.sbc_deshelve_time = float(self.p['SidebandCooling.DeshelvingTime'])
        self.use_sbc = float(self.p['SidebandCooling.UseSidebandCooling'])
        self.sbc_sbc_cycles = float(self.p['SidebandCooling.SidebandCycles'])

        self.sbc_r_dur = float(self.p['RamanSidebandCooling.Duration'])
        self.sbc_r_rf_amp_scan = float(self.p['RamanSidebandCooling.RF_Scan_Amp'])
        self.sbc_r_rf_amp_static = float(self.p['RamanSidebandCooling.RF_Static_Amp'])
        self.sbc_r_rf_freq_scan = float(self.p['RamanSidebandCooling.Scan_AOM_Frequency'])
        self.sbc_r_rf_freq_static = float(self.p['RamanSidebandCooling.Static_AOM_Frequency'])
        self.sbc_r_sp_power = float(self.p['RamanSidebandCooling.StatePrepPower'])
        self.sbc_r_repump_power = float(self.p['RamanSidebandCooling.RepumpPower'])
        self.use_raman_sbc_r = float(self.p['RamanSidebandCooling.UseSidebandCooling'])
        self.sbc_r_pulsed = float(self.p['RamanSidebandCooling.Pulsed'])
        self.sbc_r_pulse_dur = float(self.p['RamanSidebandCooling.Pulse_Duration'])
        self.sbc_r_op_dur = float(self.p['RamanSidebandCooling.OP_Duration'])
        self.sbc_r_cycles = int(self.p['RamanSidebandCooling.Pulse_Cycles'])

        self.MM_dur = float(self.p['MicrowavesSub.Pulse_Duration'])
        self.MM_Freq = float(self.p['MicrowavesSub.Frequency'])
        self.MM_Amp = float(self.p['MicrowavesSub.RF_Amp'])
        self.MM_Att = int(self.p['MicrowavesSub.RF_Att'])
        self.phases = []


        self.pulse_duration = float(self.p['time'])
        self.dc_counts = [0]*self.cycles
        self.sd_counts = [0]*self.cycles
        self.total_counts = [[0]*self.cycles]*2


