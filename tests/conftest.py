import pytest
import numpy as np
from mesh2scattering import input


@pytest.fixture
def source_coords_10deg():
    source_azimuth_deg = np.arange(0, 95, 10)
    source_colatitude_deg = np.arange(10, 85, 10)
    source_radius = 10

    return input.create_source_positions(
        source_azimuth_deg, source_colatitude_deg, source_radius)
