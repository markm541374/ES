import scipy as sp
import ctypes as ct
libGP = ct.cdll.LoadLibrary("./libGPc.so")
libGP.k.restype = ct.c_double

class GP_LKonly:
    def __init__(self, X_s, Y_s, S_s, D_s, kf):
        [n ,D] = X_s.shape
        R = ct.c_double()
        Dc = [0 if sp.isnan(x[0]) else int(sum([8**i for i in x])) for x in D_s]
        libGP.newGP_LKonly(ct.c_int(D),ct.c_int(n),X_s.ctypes.data_as(ct.POINTER(ct.c_double)),Y_s.ctypes.data_as(ct.POINTER(ct.c_double)),S_s.ctypes.data_as(ct.POINTER(ct.c_double)),(ct.c_int*len(Dc))(*Dc),kf.hyp.ctypes.data_as(ct.POINTER(ct.c_double)),ct.byref(R))
        self.l = R.value
        return
    
    def llk(self):
        return self.l


class GPcore:
    def __init__(self, X_s, Y_s, S_s, D_s, kf):
        
        [self.n ,self.D] = X_s.shape
        self.s = libGP.newGP(ct.c_int(self.D),ct.c_int(self.n))
        self.Y_s=Y_s
        libGP.set_X(self.s,X_s.ctypes.data_as(ct.POINTER(ct.c_double)))
        libGP.set_Y(self.s,Y_s.ctypes.data_as(ct.POINTER(ct.c_double)))
        libGP.set_S(self.s,S_s.ctypes.data_as(ct.POINTER(ct.c_double)))
        D = [0 if sp.isnan(x[0]) else int(sum([8**i for i in x])) for x in D_s]
        
        libGP.set_D(self.s,(ct.c_int*len(D))(*D))
        libGP.set_hyp(self.s,kf.hyp.ctypes.data_as(ct.POINTER(ct.c_double)))
        libGP.build_K(self.s)
        libGP.fac(self.s)
        libGP.presolv(self.s)
        return
    
    def __del__(self):
        libGP.killGP(self.s)
        return
    
    def printc(self):
        libGP.ping(self.s)
        return
    
    def infer_m(self,X_i,D_i):
        ns=X_i.shape[0]
        D = [0 if sp.isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.matrix(sp.empty(ns)).T
        libGP.infer_m(self.s,ns,X_i.ctypes.data_as(ct.POINTER(ct.c_double)),(ct.c_int*len(D))(*D),R.ctypes.data_as(ct.POINTER(ct.c_double)))
        return R
    
    def infer_full(self,X_i,D_i):
        ns=X_i.shape[0]
        D = [0 if sp.isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.matrix(sp.empty([ns+1,ns]))
        libGP.infer_full(self.s,ns,X_i.ctypes.data_as(ct.POINTER(ct.c_double)),(ct.c_int*len(D))(*D),R.ctypes.data_as(ct.POINTER(ct.c_double)))
        m = R[0,:].T
        V = R[1:,:]
        return [m,V]
    
    def infer_diag(self,X_i,D_i):
        ns=X_i.shape[0]
        D = [0 if sp.isnan(x[0]) else int(sum([8**i for i in x])) for x in D_i]
        R=sp.matrix(sp.empty([2,ns]))
        libGP.infer_diag(self.s,ns,X_i.ctypes.data_as(ct.POINTER(ct.c_double)),(ct.c_int*len(D))(*D),R.ctypes.data_as(ct.POINTER(ct.c_double)))
        m = R[0,:].T
        V = R[1,:].T
        return [m,V]

    def llk(self):
        R = ct.c_double()
        libGP.llk(self.s,ct.byref(R))
        return R.value
#kf = gen_sqexp_k_d([1.,0.3])




class gen_sqexp_k_d():
    def __init__(self,theta):
        self.dim = len(theta)-1
        self.hyp = sp.array(theta)
        self.hypinv = sp.array([1./x**2 for x in theta])
        self.hypinv[0] = theta[0]**2
        return
    def __call__(self,x1, x2, d1=[sp.NaN], d2=[sp.NaN]):
        D1 = 0 if sp.isnan(d1[0]) else int(sum([8**x for x in d1]))
        D2 = 0 if sp.isnan(d2[0]) else int(sum([8**x for x in d2]))
        r=libGP.k(x1.ctypes.data_as(ct.POINTER(ct.c_double)),x2.ctypes.data_as(ct.POINTER(ct.c_double)), ct.c_int(D1),ct.c_int(D2),ct.c_int(self.dim),self.hypinv.ctypes.data_as(ct.POINTER(ct.c_double)))
        return r
    

