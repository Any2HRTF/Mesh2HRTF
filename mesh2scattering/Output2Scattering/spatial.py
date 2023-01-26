import numpy as np
import pyfar as pf


def angles2coords(
        azimuth, colatitude,
        radius: float = 1., unit='rad') -> pf.Coordinates:
    """
    ``data.cshape`` fits the cshape of ```coords``. Data get shifed throght
    the ``coords`` Object around azimuth by ``shift_azimuth``.
    """
    azimuth = np.array(azimuth)
    colatitude = np.array(colatitude)
    if unit == 'deg':
        azimuth = azimuth * np.pi / 180.
        colatitude = colatitude * np.pi / 180.
    elif unit != 'rad':
        raise TypeError("Unknown Unit")
    phi, theta = np.meshgrid(azimuth, colatitude, indexing='ij')
    return pf.Coordinates(
        phi, theta, np.ones(phi.shape)*radius, 'sph')


def shift_data_coords(
        data: pf.FrequencyData, coords: pf.Coordinates, shift_azimuth: float
        ) -> pf.FrequencyData:
    """
    ``data.cshape`` fits the cshape of ```coords``. Data get shifed throght
    the ``coords`` Object around azimuth by ``shift_azimuth``.
    """
    # test input
    if not isinstance(data, pf.FrequencyData):
        raise TypeError(
            f'Data should be of type FrequencyData not {type(data)}')
    if not isinstance(coords, pf.Coordinates):
        raise TypeError(
            f'coords should be of type Coordinates not {type(coords)}')
    if not isinstance(shift_azimuth, (float, int)):
        raise TypeError(
            f'shift_azimuth should be of type float not {type(shift_azimuth)}')

    if shift_azimuth == 0:
        return data.copy()
    coords_ref = coords.copy()
    coords_cp = coords.copy()
    sph = coords_cp.get_sph(unit='deg')
    # shift azimuth by shift_azimuth in deg
    azimuth = np.remainder(sph[..., 0] + shift_azimuth, 360)
    coords_cp.set_sph(azimuth, sph[..., 1], sph[..., 2], unit='deg')
    xyz = coords_ref.get_cart()
    data_mask, _ = coords_cp.find_nearest_k(xyz[..., 0], xyz[..., 1], xyz[..., 2])
    data_mask = data_mask.flatten()
    freq = data.freq.copy()
    freq_new = freq[:, data_mask, ...]
    data_out = pf.FrequencyData(
        freq_new, data.frequencies)
    return data_out


def reshape_to_az_by_el(
        data: pf.FrequencyData, coords_in: pf.Coordinates,
        coords_out: pf.Coordinates, cdim: int = 0) -> (pf.FrequencyData):
    if cdim > 0:
        data.freq = np.moveaxis(data.freq, cdim, 0)
    freq_shape = list(coords_out.cshape)
    if len(data.cshape) > 1:
        for dim in data.cshape[1:]:
            freq_shape.append(dim)
    freq_shape.append(data.n_bins)
    freq = np.zeros(freq_shape, dtype=complex)
    data_in = data.freq
    xyz = coords_out.get_cart()
    index, _ = coords_in.find_nearest_k(xyz[..., 0], xyz[..., 1], xyz[..., 2])
    for iaz in range(coords_out.cshape[0]):
        res_data = data_in[index[iaz, :], ...]
        freq[iaz, ...] = res_data
    if cdim > 0:
        freq = np.moveaxis(freq, 0, cdim+1)
        freq = np.moveaxis(freq, 0, cdim+1)
    data_out = pf.FrequencyData(freq, data.frequencies)
    return data_out


def apply_symmetry_circular(
        data: pf.FrequencyData, coords_mic: pf.Coordinates,
        coords_inc: pf.Coordinates, coords_inc_out: pf.Coordinates):
    """apply semmetry for cicular symmetrical surfaces.

    Parameters
    ----------
    data : pf.FrequencyData
        data which is rotated, cshape need to be (#theta_coords_inc,
        #coords_mic)
    coords_mic : pf.Coordinates
        Coordinate object from the receiver positions of the current reference
        plate of cshape (#theta_coords_inc)
    coords_inc : pf.Coordinates
        Coordinate object from the source positions of the reference of cshape
        (#coords_inc_reference.
    coords_inc_out : pf.Coordinates
        Coordinate object from the source positions of the sample of cshape
        (#coords_inc_sample).

    Returns
    -------
    pf.FrequencyData
        _description_
    """
    shape = [coords_inc_out.cshape[0], data.cshape[1], len(data.frequencies)]
    freq = np.empty(shape, dtype=complex)
    thetas = coords_inc.get_sph(unit='deg')[..., 1]
    for ii in range(coords_inc_out.cshape[0]):
        az = coords_inc_out[ii].get_sph(unit='deg')[0, 0]
        theta = coords_inc_out[ii].get_sph(unit='deg')[0, 1]
        data_in = data[np.abs(thetas-theta) < 1e-5, :]
        freq[ii, ...] = shift_data_coords(
            data_in, coords_mic, float(az)).freq.copy()
    data_out = pf.FrequencyData(freq, data.frequencies)
    return data_out


def apply_symmetry_mirror(data, coords_mic, incident_coords, mirror_axe=None):
    shape = list(data.cshape)
    shape.append(data.n_bins)
    shape[mirror_axe] = 2*shape[mirror_axe]-1
    index_min = int(shape[mirror_axe]/2)
    index_max = int(shape[mirror_axe])
    freq = np.empty(shape, dtype=complex)
    freq[:] = np.nan
    freq = np.moveaxis(freq, mirror_axe, -1)
    freq_in = np.moveaxis(data.freq, mirror_axe, -1)
    azimuths = incident_coords.get_sph()[:, 1, 0]
    max_aimuth = np.max(azimuths)
    elevations = incident_coords.get_sph()[0, :, 1]
    radius = np.median(incident_coords.get_sph()[:, :, 2])
    azimuths_new = []
    max_index = -1
    for iaz in range(index_max):
        if iaz > index_min:
            idx = index_max-iaz-1
            az = (max_aimuth-azimuths[idx]) * 2
            if azimuths[idx] + az > 2 * np.pi:
                max_index = iaz
                break
            data_in = shift_data_coords(data, coords_mic, az/np.pi*180).freq
            azimuths_new.append(azimuths[idx] + az)
        else:
            data_in = data.freq
            idx = iaz
            azimuths_new.append(azimuths[iaz])
        freq_in = np.moveaxis(data_in, mirror_axe, -1)
        freq[..., iaz] = freq_in[..., idx]
    if max_index > 0:
        freq = freq[..., :max_index]
    freq = np.moveaxis(freq, -1, mirror_axe)
    data_out = pf.FrequencyData(freq, data.frequencies)
    new_inc_coords = angles2coords(np.array(azimuths_new), elevations, radius)
    return data_out, new_inc_coords
