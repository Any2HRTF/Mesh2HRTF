import pytest
import os
import numpy as np
import mesh2scattering as m2s
import numpy.testing as npt


def test_import():
    from mesh2scattering import output
    assert output


def test_read_ram_estimates():

    estimates = m2s.output.read_ram_estimates(os.path.join(
        os.path.dirname(__file__), "resources", "SHTF", "NumCalc", "source_1"))

    assert isinstance(estimates, np.ndarray)
    assert estimates.shape == (60, 3)
    npt.assert_allclose([1.00000e+00, 1.00000e+02, 4.16414e-02], estimates[0])
    npt.assert_allclose([6.00000e+01, 6.00000e+03, 7.22010e-02], estimates[-1])


def test_read_ram_estimates_assertions():
    """test assertions for read_ram_estimates"""

    with pytest.raises(ValueError, match="does not contain a Memory.txt"):
        m2s.output.read_ram_estimates(os.getcwd())
