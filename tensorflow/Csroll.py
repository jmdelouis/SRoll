import ctypes
from numpy.ctypeslib import ndpointer


class Csroll:
    def __init__(self,path):
        so_file = "%s/libSRoll4.so"%(path)

        my_functions = ctypes.CDLL(so_file)

        self.finit=my_functions.init_buffer
        self.finit.restype = ctypes.c_long
        self.finit.argtypes = [ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,ctypes.c_int,
                               ctypes.c_int,ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),ctypes.c_int,ctypes.c_int]

        self.finitSpline1D=my_functions.InitSpline1D
        self.finitSpline1D.argtypes = [ctypes.c_long,ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double]

        self.fallocbolo = my_functions.allocbolo
        self.fallocbolo.argtypes = [ctypes.c_long]
        self.fcount_data_per_pixel = my_functions.count_data_per_pixel
        self.fcount_data_per_pixel.argtypes = [ctypes.c_long,
                                               ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                               ctypes.c_int]

        self.falloccorrection = my_functions.alloccorrection
        self.falloccorrection.argtypes = [ctypes.c_long,
                                          ctypes.c_int,
                                          ctypes.c_int]

        self.fsetcorrection = my_functions.setcorrection
        self.fsetcorrection.argtypes = [ctypes.c_long,
                                        ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                        ctypes.c_int,
                                        ctypes.c_int]

        self.fgetdataidx = my_functions.getdataidx
        self.fgetdataidx.argtypes = [ctypes.c_long,
                                     ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
                                     ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]
        
        self.fgetdata = my_functions.getdata
        self.fgetdata.argtypes = [ctypes.c_long,
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ctypes.c_int,
                                  ctypes.c_int,
                                  ctypes.c_int]
        
        self.fgetidxdata = my_functions.getidxdata
        self.fgetidxdata.argtypes = [ctypes.c_long,
                                     ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                     ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
                                     ctypes.c_int]

        self.getndata = my_functions.getndata
        self.getndata.argtypes = [ctypes.c_long,
                                  ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]

        self.fsetbolo = my_functions.setbolo
        self.fsetbolo.argtypes = [ctypes.c_long,
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                  ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                                  ctypes.c_int,
                                  ctypes.c_int]

        self.fsetbolo_gain = my_functions.setbolo_gain
        self.fsetbolo_gain.argtypes = [ctypes.c_long,
                                       ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                       ctypes.c_int,
                                       ctypes.c_int]

        self.frelativ = my_functions.setrelativ
        self.frelativ.argtypes = [ctypes.c_long,
                                  ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]
        self.finitgrad = my_functions.initgrad
        self.finitgrad.argtypes = [ctypes.c_long]

        self.projdata=my_functions.projdata
        self.projdata.argtypes = [ctypes.c_long,ctypes.c_double,
                                  ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

        self.calc_A0=my_functions.calc_A0
        self.calc_A0.argtypes = [ctypes.c_long,ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

        self.proj=my_functions.proj
        self.proj.argtypes = [ctypes.c_long,
                              ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                              ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

        self.domap=my_functions.domap
        self.domap.argtypes = [ctypes.c_long,
                               ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                               ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
                               ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
                               ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                               ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]

        self.ang2pix=my_functions.ang2pix_ring_table
        self.ang2pix.argtypes = [ctypes.c_long,
                                 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                 ndpointer(ctypes.c_long, flags="C_CONTIGUOUS"),
                                 ctypes.c_long]

        self.pix2ang=my_functions.pix2ang_ring_table
        self.pix2ang.argtypes = [ctypes.c_long,
                                 ndpointer(ctypes.c_long, flags="C_CONTIGUOUS"),
                                 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                 ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                                 ctypes.c_long]

        self.writemap = my_functions.write_healpix_map
        self.writemap.argtypes = [ndpointer(ctypes.c_float, flags="C_CONTIGUOUS"),
                                  ctypes.c_long,
                                  ctypes.c_char_p,
                                  ctypes.c_int8,
                                  ctypes.c_char_p]

        self.UNSEEN = -1.6375e+30

        self.nocalib = my_functions.nocalib 
        self.nocalib.argtypes = [ctypes.c_long,
                                 ndpointer(ctypes.c_int, flags="C_CONTIGUOUS")]

        self.docalib = my_functions.docalib 
        self.docalib.argtypes = [ctypes.c_long]
