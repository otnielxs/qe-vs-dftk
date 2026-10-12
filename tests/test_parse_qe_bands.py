from pathlib import Path

import numpy as np

from qevsdftk.parse_qe import parse_qe_bands

FIXTURE = Path(__file__).parent / "fixtures" / "si_bands.out"


def test_parse_qe_bands_shape_and_order():
    bands = np.array(parse_qe_bands(FIXTURE, 8))
    assert bands.shape == (101, 8)                 # 5 segments x 20 points + 1
    assert np.all(np.diff(bands, axis=1) >= -1e-9)  # bands sorted by energy
