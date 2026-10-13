from qevsdftk.parse_qe import _parse_clock, parse_qe_wall_time


def test_parse_clock_formats():
    assert abs(_parse_clock("0.94s") - 0.94) < 1e-9
    assert abs(_parse_clock("1m 5.20s") - 65.2) < 1e-9
    assert abs(_parse_clock("1h 2m 3.00s") - 3723.0) < 1e-9


def test_parse_qe_wall_time(tmp_path):
    f = tmp_path / "x.out"
    f.write_text("     PWSCF        :      0.83s CPU      1m 2.50s WALL\n")
    assert abs(parse_qe_wall_time(f) - 62.5) < 1e-9
