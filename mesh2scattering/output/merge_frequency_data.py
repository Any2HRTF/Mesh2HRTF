"""
Python tools for Mesh2HRTF including functions to generate SOFA files
containing the HRTF/HRIR data, merge SOFA files containing data for the left
and right ear and generate evaluation grids.
"""
import numpy as np
import pyfar as pf


def merge_frequency_data(data_list):
    data_out = data_list[0].copy()
    frequencies = data_out.frequencies.copy()
    for idx in range(1, len(data_list)):
        data = data_list[idx]
        assert data_out.cshape == data.cshape
        frequencies = np.append(frequencies, data.frequencies)
        frequencies = np.array([i for i in set(frequencies)])
        frequencies = np.sort(frequencies)

        data_new = []
        for f in frequencies:
            if any(data_out.frequencies == f):
                freq_index = np.where(data_out.frequencies == f)
                data_new.append(data_out.freq[..., freq_index[0][0]])
            elif any(data.frequencies == f):
                freq_index = np.where(data.frequencies == f)
                data_new.append(data.freq[..., freq_index[0][0]])

        data_new = np.moveaxis(np.array(data_new), 0, -1)
        data_out = pf.FrequencyData(data_new, frequencies)
    return data_out
