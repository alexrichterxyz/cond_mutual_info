#!/usr/bin/env python3
import ctypes
import pandas as pd
import numpy as np
import math
import pathlib

info = ctypes.cdll.LoadLibrary(pathlib.Path(__file__).parent.parent / 'build/info.so')


info.make_dblvecvec.restype = ctypes.c_void_p
info.make_dblvecvec.argtypes = None

info.delete_dblvecvec.restype = None
info.delete_dblvecvec.argtypes = [ctypes.c_void_p]

info.attach_dblvec.restype = None
info.attach_dblvec.argtypes = [
	ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t]

info.cond_mutual_info.restype = None
info.cond_mutual_info.argtypes = [ctypes.c_void_p, ctypes.c_void_p,
								  ctypes.c_void_p, ctypes.c_size_t, ctypes.c_double, ctypes.c_void_p, ctypes.c_void_p]


def cond_mutual_info(data=None, x=None, y=None, z=None, p_samples=10000, base=math.e):
	# design matrices
	x_array = None
	y_array = None
	z_array = None

	def parse_df_fields(fields):
		selection = data[fields].values.T

		if selection.ndim == 1:
			selection = np.expand_dims(selection, 0)

		if selection.dtype == int:
			selection = selection.astype(float)

		assert selection.dtype == float, 'Data must be float'
		return selection

	def parse_iterable(d):
		if isinstance(d, pd.DataFrame):
			return parse_iterable(d.values)
		elif isinstance(d, np.ndarray):
			if d.dtype == int:
				d = d.astype(float)

			assert d.dtype == float, 'Data must be float'

			if d.ndim == 1:
				return np.expand_dims(d, 0).T
			elif d.ndim == 2:
				return d.T
			else:
				raise 'Data must be 1 or 2D'
		elif isinstance(d, list):
			d = np.array(d)

			if d.ndim == 1:
				d = np.expand_dims(d, 0)
				
			return parse_iterable(d.T)
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

	print(x_array.shape)
	print(y_array.shape)
	print(z_array.shape)
	assert x_array.shape[1] == y_array.shape[1] and \
		y_array.shape[1] == z_array.shape[1], 'Data must have the same length'

	def array_to_dblvec(array):
		dblvec = info.make_dblvecvec()
		size = array.shape[1]
		c_double_arr = ctypes.c_double * size
		for a in array:
			info.attach_dblvec(dblvec, (c_double_arr)(*a.tolist()), size)
		return dblvec

	x_vec = array_to_dblvec(x_array)
	y_vec = array_to_dblvec(y_array)
	z_vec = array_to_dblvec(z_array)

	c_p_val = ctypes.c_double()
	c_cmi = ctypes.c_double()
	c_p_samples = ctypes.c_size_t(p_samples)
	c_base = ctypes.c_double(base)

	info.cond_mutual_info(x_vec, y_vec, z_vec, c_p_samples,
						  c_base, ctypes.byref(c_cmi), ctypes.byref(c_p_val))

	info.delete_dblvecvec(x_vec)
	info.delete_dblvecvec(y_vec)
	info.delete_dblvecvec(z_vec)

	return c_cmi.value, c_p_val.value
