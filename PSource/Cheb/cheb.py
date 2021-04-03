"""
Python wrapper functions for taking Chebyshev pseudospectral derivatives
over the interval [lower, upper].
"""
import os, ctypes
import numpy as np

_p_cheb = os.getcwd()
_l_cheb = ctypes.CDLL(_p_cheb+'/lib_cheb.so')

_cheb_initialized = False
#=============================================================================
_cheb_init          = _l_cheb.init 
_cheb_init.restype  = None
_cheb_init.argtypes = [
      ctypes.c_size_t,
      ctypes.c_double, 
      ctypes.c_double
      ]

def init(n:int, lower:float, upper:float) -> None:
   """
   Initialize internal FFTW pointers, etc.
   """
   global _cheb_initialized
   _cheb_initialized = True
   _cheb_init(
         ctypes.c_size_t(n), 
         ctypes.c_double(lower), 
         ctypes.c_double(upper)
         )
#=============================================================================
_cheb_cleanup          = _l_cheb.cleanup
_cheb_cleanup.restype  = None
_cheb_cleanup.argtypes = []

def cleanup() -> None:
   """
   Free internal FFTW pointers, etc.
   """
   global _cheb_initialized
   assert(_cheb_initialized)
   _cheb_initialized = False
   _cheb_cleanup()
#=============================================================================
_cheb_n          = _l_cheb.n
_cheb_n.restype  = ctypes.c_size_t
_cheb_n.argtypes = []

def n() -> int:
   """
   Number of Chebyshev collocation points.
   """
   assert(_cheb_initialized)
   return _cheb_n()
#=============================================================================
_cheb_lower          = _l_cheb.lower
_cheb_lower.restype  = ctypes.c_double
_cheb_lower.argtypes = []

def lower() -> float:
   """
   Lower boundary of domain in real space.
   """
   assert(_cheb_initialized)
   return _cheb_lower()
#=============================================================================
_cheb_upper          = _l_cheb.upper
_cheb_upper.restype  = ctypes.c_double
_cheb_upper.argtypes = []

def upper() -> float:
   """
   Upper boundary of domain in real space.
   """
   assert(_cheb_initialized)
   return _cheb_upper()
#=============================================================================
_cheb_pt          = _l_cheb.pt
_cheb_pt.restype  = ctypes.c_double
_cheb_pt.argtypes = [
      ctypes.c_size_t
      ]

def pt(i:int) -> float:
   """
   Location of ith Chebyshev point in real space.
   """
   assert(_cheb_initialized)
   assert(i>=0 and i<n())
   return _cheb_pt(ctypes.c_size_t(i))
#=============================================================================
_cheb_der          = _l_cheb.der
_cheb_der.restype  = None 
_cheb_der.argtypes = [
      np.ctypeslib.ndpointer(ctypes.c_double),
      np.ctypeslib.ndpointer(ctypes.c_double)
      ]

def der(v:np.array, dv:np.array) -> None:
   """
   Derivative over interval [lower, upper].
   """
   assert(_cheb_initialized)
   assert(v.size==n() and dv.size==n())
   _cheb_der(v,dv)
#=============================================================================
_cheb_filter          = _l_cheb.filter
_cheb_filter.restype  = None
_cheb_filter.argtypes = [
      np.ctypeslib.ndpointer(ctypes.c_double)
      ]

def filter(v:np.array) -> None:
   """
   Low pass filter in Chebyshev space.
   """
   assert(_cheb_initialized)
   assert(v.size==n())
   _cheb_filter(v)
