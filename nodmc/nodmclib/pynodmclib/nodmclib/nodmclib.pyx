###############################################################################
# "nodmclib.pyx"
# This is the NO_DMC_LIB wrapper for python.
#
# Author: C. Madsen
# Date: 2015-01-26
###############################################################################

"""This is the Python wrapper of the NO_DMC_LIB."""

import numpy as np
cimport numpy as np
cimport nodmclib.cno_dmc_data_access

# This is the version of this wrapper NOT the version of the NO_DMC_LIB!
__version__ = "1.0"

#-----------------------------------------------------------------------------

# Quite clause to the former piolib.InfoObject()
def getObjectInfo(object_name):
  """Allow to retrieve object information.
  This function is quite close to the former piolib.InfoObject().

  Args:
    object_name: the dmc object name (full path)

  Returns:
    the list of information associated to this object.
    In order (zero indexed): 
    0] TOItype
    1] Datatype
    2] BeginIndex
    3] EndIndex
    4] Author
    5] Date
    (optional) 6] Keyword list

  Exemple:
    >>> getObjectInfo("/m3gpfs3/datadmc/dmc/MISS03/DATA/calTOIs/02_143_1a_LFER6_JC_v64")
    ('TOI',
     'PIOFLOAT',
      1371347725,
      15147406766,
      240,
      27005,
      'opsman',
      '2013-12-07 10:57:55.123896')
  """
  cdef nodmclib.cno_dmc_data_access.metadata met

  # Call the corresponding C function
  res = nodmclib.cno_dmc_data_access.getMetadataFor(object_name, &met)

  # Test result
  if res != 0:
    raise Exception("Error: C lib return an error (see previous message)")

  # Construct the tuple to be returned
  return (met.PIOtype, met.Datatype, met.BeginIndex, met.EndIndex, met.BeginRing, met.EndRing, met.Author, met.Date)

#-----------------------------------------------------------------------------

def read_PIOBYTE(object_name, offset, nbsample):
  """Allow to read data of PIOBYTE type.

  Args:
    object_name (str): the dmc object name (full path)
    offset (int): the offset from which we start to read data
    nbsample (int): number of sample to be read

  Raises:
    Exception: in case of error.

  Returns:
    the read data if successful, otherwise raise an exception.
  """
  # Allocate required memory for data
  cdef np.ndarray data = np.zeros(nbsample, dtype=np.int8)

  # Call the corresponding C function
  res = nodmclib.cno_dmc_data_access.noDMC_readObject_PIOBYTE(object_name, offset, nbsample, <char *>data.data)

  # Test result
  if res != 0:
    raise Exception("Error: C lib return an error (see previous message)")

  # return data
  return data

#-----------------------------------------------------------------------------

def read_PIOFLAG(object_name, offset, nbsample):
  """Allow to read data of PIOFLAG type.

  Args:
    object_name (str): the dmc object name (full path)
    offset (int): the offset from which we start to read data
    nbsample (int): number of sample to be read

  Raises:
    Exception: in case of error.

  Returns:
    the read data if successful, otherwise raise an exception.
  """
  # Allocate required memory for data
  cdef np.ndarray data = np.zeros(nbsample, dtype=np.bool_)

  # Call the corresponding C function
  res = nodmclib.cno_dmc_data_access.noDMC_readObject_PIOFLAG(object_name, offset, nbsample, <unsigned char *>data.data)

  # Test result
  if res != 0:
    raise Exception("Error: C lib return an error (see previous message)")

  # return data
  return data

#-----------------------------------------------------------------------------

def read_PIOSHORT(object_name, offset, nbsample):
  """Allow to read data of PIOSHORT type.

  Args:
    object_name (str): the dmc object name (full path)
    offset (int): the offset from which we start to read data
    nbsample (int): number of sample to be read

  Raises:
    Exception: in case of error.

  Returns:
    the read data if successful, otherwise raise an exception.
  """
  # Allocate required memory for data
  cdef np.ndarray data = np.zeros(nbsample, dtype=np.int16)

  # Call the corresponding C function
  res = nodmclib.cno_dmc_data_access.noDMC_readObject_PIOSHORT(object_name, offset, nbsample, <short *>data.data)

  # Test result
  if res != 0:
    raise Exception("Error: C lib return an error (see previous message)")

  # return data
  return data

#-----------------------------------------------------------------------------

def read_PIOINT(object_name, offset, nbsample):
  """Allow to read data of PIOINT type.

  Args:
    object_name (str): the dmc object name (full path)
    offset (int): the offset from which we start to read data
    nbsample (int): number of sample to be read

  Raises:
    Exception: in case of error.

  Returns:
    the read data if successful, otherwise raise an exception.
  """
  # Allocate required memory for data
  cdef np.ndarray data = np.zeros(nbsample, dtype=np.int32)

  # Call the corresponding C function
  res = nodmclib.cno_dmc_data_access.noDMC_readObject_PIOINT(object_name, offset, nbsample, <int *>data.data)

  # Test result
  if res != 0:
    raise Exception("Error: C lib return an error (see previous message)")

  # return data
  return data

#-----------------------------------------------------------------------------

def read_PIOLONG(object_name, offset, nbsample):
  """Allow to read data of PIOLONG type.

  Args:
    object_name (str): the dmc object name (full path)
    offset (int): the offset from which we start to read data
    nbsample (int): number of sample to be read

  Raises:
    Exception: in case of error.

  Returns:
    the read data if successful, otherwise raise an exception.
  """
  # Allocate required memory for data
  cdef np.ndarray data = np.zeros(nbsample, dtype=np.int64)

  # Call the corresponding C function
  res = nodmclib.cno_dmc_data_access.noDMC_readObject_PIOLONG(object_name, offset, nbsample, <long *>data.data)

  # Test result
  if res != 0:
    raise Exception("Error: C lib return an error (see previous message)")

  # return data
  return data

#-----------------------------------------------------------------------------

def read_PIOFLOAT(object_name, offset, nbsample):
  """Allow to read data of PIOFLOAT type.

  Args:
    object_name (str): the dmc object name (full path)
    offset (int): the offset from which we start to read data
    nbsample (int): number of sample to be read

  Raises:
    Exception: in case of error.

  Returns:
    the read data if successful, otherwise raise an exception.
  """
  # Allocate required memory for data
  cdef np.ndarray data = np.zeros(nbsample, dtype=np.float32)

  # Call the corresponding C function
  res = nodmclib.cno_dmc_data_access.noDMC_readObject_PIOFLOAT(object_name, offset, nbsample, <float *>data.data)

  # Test result
  if res != 0:
    raise Exception("Error: C lib return an error (see previous message)")

  # return data
  return data

#-----------------------------------------------------------------------------

def read_PIODOUBLE(object_name, offset, nbsample):
  """Allow to read data of PIODOUBLE type.

  Args:
    object_name (str): the dmc object name (full path)
    offset (int): the offset from which we start to read data
    nbsample (int): number of sample to be read

  Raises:
    Exception: in case of error.

  Returns:
    the read data if successful, otherwise raise an exception.
  """
  # Allocate required memory for data
  cdef np.ndarray data = np.zeros(nbsample, dtype=np.float64)

  # Call the corresponding C function
  res = nodmclib.cno_dmc_data_access.noDMC_readObject_PIODOUBLE(object_name, offset, nbsample, <double *>data.data)

  # Test result
  if res != 0:
    raise Exception("Error: C lib return an error (see previous message)")

  # return data
  return data

#-----------------------------------------------------------------------------

def read_PIOSTRING(object_name, offset, nbsample):
  """Allow to read data of PIOSTRING type.

  Args:
    object_name (str): the dmc object name (full path)
    offset (int): the offset from which we start to read data
    nbsample (int): number of sample to be read

  Raises:
    Exception: in case of error.

  Returns:
    the read data if successful, otherwise raise an exception.
  """
  # Allocate required memory for data
  cdef np.ndarray data = np.zeros(nbsample, dtype='S256')

  # Call the corresponding C function
  res = nodmclib.cno_dmc_data_access.noDMC_readObject_PIOSTRING(object_name, offset, nbsample, <char **>data.data)

  # Test result
  if res != 0:
    raise Exception("Error: C lib return an error (see previous message)")

  # return data
  return data

#-----------------------------------------------------------------------------

def readFlag_Written(object_name, offset, nbsample):
  """Allow to read "Written flags" associated to the specified object (any type).

  Args:
    object_name (str): the dmc object name (full path)
    offset (int): the offset from which we start to read flags
    nbsample (int): number of sample to be read

  Raises:
    Exception: in case of error.

  Returns:
    the read flags if successful, otherwise raise an exception.
  """
  # Allocate required memory for data
  cdef np.ndarray data = np.zeros(nbsample, dtype=np.bool_)

  # Call the corresponding C function
  res = nodmclib.cno_dmc_data_access.noDMC_readFlagWritten(object_name, offset, nbsample, <unsigned char *>data.data)

  # Test result
  if res != 0:
    raise Exception("Error: C lib return an error (see previous message)")

  # return data
  return data
