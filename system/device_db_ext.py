# Custom devices using the LAXDevice class

device_db_ext = {



    # LABRAD Devices
    "aperture": {
        "type": "local",
        "module": "LAX_exp.system.devices.aperture",
        "class": "Aperture"
    },

    "camera": {
        "type": "local",
        "module": "LAX_exp.system.devices.andor_camera",
        "class": "AndorCamera"
    },

    # BEAMS
    "probe": {
        "type": "local",
        "module": "LAX_exp.system.devices.beam_397_probe",
        "class": "Beam397Probe"
    },
    "pump": {
        "type": "local",
        "module": "LAX_exp.system.devices.beam_397_pump",
        "class": "Beam397Pump"
    },
    "repump_cooling": {
        "type": "local",
        "module": "LAX_exp.system.devices.beam_866",
        "class": "Beam866"
    },
    "repump_qubit": {
        "type": "local",
        "module": "LAX_exp.system.devices.beam_854",
        "class": "Beam854"
    },
    "qubit": {
        "type": "local",
        "module": "LAX_exp.system.devices.beam_729",
        "class": "Beam729"
    },


    # DDS
    "dds_dipole": {
        "type": "local",
        "module": "LAX_exp.system.devices.dds_dipole",
        "class": "DDSDipole"
    },
    "dds_parametric": {
        "type": "local",
        "module": "LAX_exp.system.devices.dds_parametric",
        "class": "DDSParametric"
    },
    "phaser_eggs": {
        "type": "local",
        "module": "LAX_exp.system.devices.phaser_eggs",
        "class": "PhaserEGGS"
    },

    # PMT
    "pmt": {
        "type": "local",
        "module": "LAX_exp.system.devices.pmt",
        "class": "PMTCounter"
    },

    # TTL
    "trigger_line": {
        "type": "local",
        "module": "LAX_exp.system.devices.trigger_line",
        "class": "TriggerLine"
    },
    "trigger_rf": {
        "type": "local",
        "module": "LAX_exp.system.devices.trigger_rf",
        "class": "TriggerRF"
    }
}
