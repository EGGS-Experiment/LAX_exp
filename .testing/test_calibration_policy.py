from LAX_exp.calibration import CalibrationPolicy, CalibrationState


def test_calibrates_when_missing():
    policy = CalibrationPolicy(calibrate_if_missing=True)
    state = CalibrationState()

    decision = policy.should_calibrate(state, now=0.0)

    assert decision.should_calibrate
    assert decision.reason == "missing"


def test_does_not_calibrate_when_missing_is_disabled():
    policy = CalibrationPolicy(calibrate_if_missing=False)
    state = CalibrationState()

    decision = policy.should_calibrate(state, now=0.0)

    assert not decision.should_calibrate


def test_calibrates_after_n_experiments():
    policy = CalibrationPolicy(every_n_experiments=3)
    state = CalibrationState(
        last_calibration_time=0.0,
        runs_since_calibration=3,
    )

    decision = policy.should_calibrate(state, now=100.0)

    assert decision.should_calibrate
    assert decision.reason == "run_count"


def test_calibrates_after_time_interval():
    policy = CalibrationPolicy(max_age_s=10.0)
    state = CalibrationState(
        last_calibration_time=0.0,
        runs_since_calibration=1,
    )

    decision = policy.should_calibrate(state, now=10.0)

    assert decision.should_calibrate
    assert decision.reason == "time_interval"


def test_calibrates_after_error():
    policy = CalibrationPolicy(calibrate_after_error=True)
    state = CalibrationState(
        last_calibration_time=0.0,
        runs_since_calibration=1,
        previous_run_failed=True,
    )

    decision = policy.should_calibrate(state, now=1.0)

    assert decision.should_calibrate
    assert decision.reason == "after_error"
