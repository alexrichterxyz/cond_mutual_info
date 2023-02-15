#!/usr/bin/env python3
import ctypes
import pandas as pd
import numpy as np

info = ctypes.cdll.LoadLibrary('./build/info.so')


info.make_dblvec_vec.restype = ctypes.c_void_p
info.make_dblvec_vec.argtypes = None

info.delete_dblvec_vec.restype = None
info.delete_dblvec_vec.argtypes = [ctypes.c_void_p]

info.attach_dblvec.restype = None
info.attach_dblvec.argtypes = [
    ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t]

info.conditional_mi.restype = None
info.conditional_mi.argtypes = [ctypes.c_void_p, ctypes.c_void_p,
                                ctypes.c_void_p, ctypes.c_size_t, ctypes.c_void_p, ctypes.c_void_p]


def conditional_mi(data=None, x=None, y=None, z=None, p_samples=100):
    # design matrices
    x_array = None
    y_array = None
    z_array = None

    def parse_df_fields(fields):
        selection = data[fields].values.T

        if selection.ndim == 1:
            selection = np.expand_dims(selection, 0)

        assert selection.dtype == float, 'Data must be float'
        return selection

    def parse_iterable(d):
        if isinstance(d, pd.DataFrame):
            return parse_iterable(d.values)
        elif isinstance(d, np.ndarray):
            assert d.dtype == float, 'Data must be float'
            if d.ndim == 1:
                return np.expand_dims(d, 0).T
            elif d.ndim == 2:
                return d.T
            else:
                raise 'Data must be 1 or 2D'
        elif isinstance(d, list):
            return parse_iterable(np.array(d))
        else:
            assert 'Invalid type'

    if isinstance(data, pd.DataFrame):
        x_array = parse_df_fields(x)
        y_array = parse_df_fields(y)
        z_array = parse_df_fields(z)
    elif data is None:
        x_array = parse_iterable(x)
        y_array = parse_iterable(y)
        z_array = parse_iterable(z)

    assert x_array.shape[1] == y_array.shape[1] and \
        y_array.shape[1] == z_array.shape[1], 'Data must have the same length'


    def array_to_dblvec(array):
        dblvec = info.make_dblvec_vec()
        size = array.shape[1]
        c_double_arr = ctypes.c_double * size
        for a in array:
            info.attach_dblvec(dblvec, (c_double_arr)(*a.tolist()), size)
        return dblvec
    
    x_vec = array_to_dblvec(x_array)
    y_vec = array_to_dblvec(y_array)
    z_vec = array_to_dblvec(z_array)

    p_val = ctypes.c_double()
    cmi_val = ctypes.c_double()
    p_samples = ctypes.c_size_t(p_samples)

    info.conditional_mi(x_vec, y_vec, z_vec, p_samples, ctypes.byref(cmi_val), ctypes.byref(p_val))

    info.delete_dblvec_vec(x_vec)
    info.delete_dblvec_vec(y_vec)
    info.delete_dblvec_vec(z_vec)

    return cmi_val.value, p_val.value