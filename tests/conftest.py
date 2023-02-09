import mesh2scattering as m2s
import pytest
import numpy as np

@pytest.fixture
def source_coords_10deg():
    source_azimuth_deg = np.arange(0, 95, 10)
    source_colatitude_deg = np.arange(10, 85, 10)
    source_radius = 10

    return m2s.input.create_source_positions(
        source_azimuth_deg, source_colatitude_deg, source_radius)




