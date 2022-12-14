# Custom devices using the LAXDevice class

device_db_ext = {

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
    "tickle": {
        "type": "local",
        "module": "LAX_exp.system.devices.beam_tickle",
        "class": "BeamTickle"
    },

    # PMT
    "pmt": {
        "type": "local",
        "module": "LAX_exp.system.devices.pmt",
        "class": "PMTCounter"
    },

    # TTL
    "linetrigger": {
        "type": "local",
        "module": "LAX_exp.system.devices.trigger_linetrigger",
        "class": "Linetrigger"
    },
    "rf_sync": {
        "type": "local",
        "module": "LAX_exp.system.devices.trigger_rf_modulation",
        "class": "RFSync"
    }
}
