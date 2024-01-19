"""State noise compensation (SNC) test module."""

import numpy as np

from iapetus.propagation.dynamics.nonlinear.astro import snc

qr = 1.2
qi = 2.3
qc = 3.4
qcd = 4.5
R = np.eye(3)


def test_qric():
    Q_fcn = snc.QRic(qr_mps2=qr, qi_mps2=qi, qc_mps2=qc)

    dt = 60
    Q = Q_fcn(dt)
    assert Q[0, 0] == qr**2 * dt**4 / 3
    assert Q[0, 3] == qr**2 * dt**3 / 2
    assert Q[1, 1] == qi**2 * dt**4 / 3
    assert Q[1, 4] == qi**2 * dt**3 / 2
    assert Q[2, 2] == qc**2 * dt**4 / 3
    assert Q[2, 5] == qc**2 * dt**3 / 2

    assert Q[3, 0] == qr**2 * dt**3 / 2
    assert Q[3, 3] == qr**2 * dt**2
    assert Q[4, 1] == qi**2 * dt**3 / 2
    assert Q[4, 4] == qi**2 * dt**2
    assert Q[5, 2] == qc**2 * dt**3 / 2
    assert Q[5, 5] == qc**2 * dt**2

    assert Q[0, 1] == 0
    assert Q[0, 2] == 0
    assert Q[0, 4] == 0
    assert Q[0, 5] == 0
    assert Q[1, 0] == 0
    assert Q[1, 2] == 0
    assert Q[1, 3] == 0
    assert Q[1, 5] == 0
    assert Q[2, 0] == 0
    assert Q[2, 1] == 0
    assert Q[2, 3] == 0
    assert Q[2, 4] == 0
    assert Q[3, 1] == 0
    assert Q[3, 2] == 0
    assert Q[3, 4] == 0
    assert Q[3, 5] == 0
    assert Q[4, 0] == 0
    assert Q[4, 2] == 0
    assert Q[4, 3] == 0
    assert Q[4, 5] == 0
    assert Q[5, 0] == 0
    assert Q[5, 1] == 0
    assert Q[5, 3] == 0
    assert Q[5, 4] == 0


def test_q_inertial_ric():
    Q_fcn = snc.QInertialRic(qr_mps2=qr, qi_mps2=qi, qc_mps2=qc)

    q_inertial = Q_fcn.q_inertial(R)
    assert q_inertial[0] == qr**2
    assert q_inertial[1] == qi**2
    assert q_inertial[2] == qc**2

    dt = 60
    Q = Q_fcn(dt, R)

    assert Q[0, 0] == qr**2 * dt**4 / 3
    assert Q[0, 3] == qr**2 * dt**3 / 2
    assert Q[1, 1] == qi**2 * dt**4 / 3
    assert Q[1, 4] == qi**2 * dt**3 / 2
    assert Q[2, 2] == qc**2 * dt**4 / 3
    assert Q[2, 5] == qc**2 * dt**3 / 2

    assert Q[3, 0] == qr**2 * dt**3 / 2
    assert Q[3, 3] == qr**2 * dt**2
    assert Q[4, 1] == qi**2 * dt**3 / 2
    assert Q[4, 4] == qi**2 * dt**2
    assert Q[5, 2] == qc**2 * dt**3 / 2
    assert Q[5, 5] == qc**2 * dt**2

    assert Q[0, 1] == 0
    assert Q[0, 2] == 0
    assert Q[0, 4] == 0
    assert Q[0, 5] == 0
    assert Q[1, 0] == 0
    assert Q[1, 2] == 0
    assert Q[1, 3] == 0
    assert Q[1, 5] == 0
    assert Q[2, 0] == 0
    assert Q[2, 1] == 0
    assert Q[2, 3] == 0
    assert Q[2, 4] == 0
    assert Q[3, 1] == 0
    assert Q[3, 2] == 0
    assert Q[3, 4] == 0
    assert Q[3, 5] == 0
    assert Q[4, 0] == 0
    assert Q[4, 2] == 0
    assert Q[4, 3] == 0
    assert Q[4, 5] == 0
    assert Q[5, 0] == 0
    assert Q[5, 1] == 0
    assert Q[5, 3] == 0
    assert Q[5, 4] == 0


def test_q_inertial_cd():
    Q_fcn = snc.QInertialCd(qr_mps2=qr, qi_mps2=qi, qc_mps2=qc, qcd=qcd)

    dt = 60
    Q = Q_fcn(dt, R)

    assert Q[6, 6] == qcd**2

    assert Q[0, 6] == 0
    assert Q[1, 6] == 0
    assert Q[2, 6] == 0
    assert Q[3, 6] == 0
    assert Q[4, 6] == 0
    assert Q[5, 6] == 0

    assert Q[6, 0] == 0
    assert Q[6, 1] == 0
    assert Q[6, 2] == 0
    assert Q[6, 3] == 0
    assert Q[6, 4] == 0
    assert Q[6, 5] == 0

    assert Q[0, 0] == qr**2 * dt**4 / 3
    assert Q[0, 3] == qr**2 * dt**3 / 2
    assert Q[1, 1] == qi**2 * dt**4 / 3
    assert Q[1, 4] == qi**2 * dt**3 / 2
    assert Q[2, 2] == qc**2 * dt**4 / 3
    assert Q[2, 5] == qc**2 * dt**3 / 2

    assert Q[3, 0] == qr**2 * dt**3 / 2
    assert Q[3, 3] == qr**2 * dt**2
    assert Q[4, 1] == qi**2 * dt**3 / 2
    assert Q[4, 4] == qi**2 * dt**2
    assert Q[5, 2] == qc**2 * dt**3 / 2
    assert Q[5, 5] == qc**2 * dt**2

    assert Q[0, 1] == 0
    assert Q[0, 2] == 0
    assert Q[0, 4] == 0
    assert Q[0, 5] == 0
    assert Q[1, 0] == 0
    assert Q[1, 2] == 0
    assert Q[1, 3] == 0
    assert Q[1, 5] == 0
    assert Q[2, 0] == 0
    assert Q[2, 1] == 0
    assert Q[2, 3] == 0
    assert Q[2, 4] == 0
    assert Q[3, 1] == 0
    assert Q[3, 2] == 0
    assert Q[3, 4] == 0
    assert Q[3, 5] == 0
    assert Q[4, 0] == 0
    assert Q[4, 2] == 0
    assert Q[4, 3] == 0
    assert Q[4, 5] == 0
    assert Q[5, 0] == 0
    assert Q[5, 1] == 0
    assert Q[5, 3] == 0
    assert Q[5, 4] == 0
