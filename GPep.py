
# Gaussian process with expectation propagation on inequality constraints

import scipy as sp
from functools import partial
from matplotlib import pyplot as plt
import scipy.stats as sps
import scipy.linalg as spl
import GPd


def catObs(O):
    return GPd.catObs(O)


def plot1(g, llimit, ulimit):
    return GPd.plot1(g, llimit, ulimit)


def gen_sqexp_k_d(theta):
    k = partial(GPd.sqexp_k_d, theta)
    return k


def PhiR(x):
    # return sp.exp(sps.norm.logpdf(x) - sps.norm.logcdf(x))
    return sps.norm.pdf(x)/sps.norm.cdf(x)


class GPcore:
    def __init__(self, X_c, Y_c, S_c, D_c, X_z, D_z, I_z, G_z, N_z, kf):

        # X,Y,S,D_c are the x,y,var and derivative that have been observed
        # X,I,G_z are the x,y sign and noise of inequality constraints
        self.X_c = X_c
        self.Y_c = Y_c
        self.S_c = S_c
        self.D_c = D_c
        self.X_z = X_z
        self.D_z = D_z
        self.I_z = I_z
        self.G_z = G_z
        self.N_z = N_z
        self.kf = kf
        self.n_z = len(G_z)
        self.n_c = len(D_z)

        self.X_s = sp.array(X_c).flatten()
        self.Y_s = sp.array(Y_c).flatten()
        self.D_s = D_c

        self.nloops = 10
        self.invalidflag = True
        return

    def runEP(self, plot='none'):
        g = GPd.GPcore(self.X_c, self.Y_c, self.S_c, self.D_c, self.kf)
        # start by making the full inference at the inequality locations
        [m0, V0] = g.infer_full(self.X_z, self.D_z)
        V0Im0 = spl.cho_solve(spl.cho_factor(V0), m0)
        V0I = V0.I  # ------------------------------explicit inverse, not good
        # create the ep observations
        yt = sp.matrix(sp.zeros([self.n_z, 1]))
        St = 10**10*sp.matrix(sp.ones([self.n_z, 1]))

        ytR = yt.copy()
        StR = St.copy()
        m_R = [[], []]
        v_R = [[], []]
        for it in xrange(self.nloops):
            for i in xrange(self.n_z):
                # update the inference at z with the ep observations
                StIyt = sp.matrix(sp.zeros([self.n_z, 1]))
                StI = sp.matrix(sp.zeros([self.n_z, self.n_z]))
                VtIV0I = V0I.copy()
                # print V0I
                for j in xrange(self.n_z):
                    StIyt[j, 0] = yt[j, 0]/float(St[j, 0])
                    VtIV0I[j, j] += 1./float(St[j, 0])

                Vp = VtIV0I.I.copy()
                mp = Vp*(StIyt+V0Im0)
                # get the inference at the specific inequality currently being updates
                v_ = 1/(1/(Vp[i, i]) - 1/(St[i, 0]))
                m_ = v_*(mp[i, 0]/Vp[i, i]-yt[i, 0]/St[i, 0])

                m_R[i].append(m_)
                v_R[i].append(v_)
                # find the new ep obs

                alpha = (m_+((-1)**self.G_z[i])*self.I_z[i, 0]) / (sp.sqrt(v_+self.N_z[i]))

                pr = PhiR(alpha)
                if sp.isnan(pr):
                    pr = -alpha
                beta = pr*(pr+alpha)/(v_+self.N_z[i]**2)
                kappa = ((-1)**self.G_z[i])*(pr+alpha) / (sp.sqrt(v_+self.N_z[i]**2))

                # replace with the new ep obs
                yt[i, 0] = m_-1./kappa
                St[i, 0] = 1./beta - v_
                if (1./beta) == sp.inf:
                    yt[i, 0] = m_
                    St[i, 0] = 10**100

#            if plot:
#                print 'iter'
#                print yt
#                print St
            ytR = sp.hstack([ytR, yt])
            StR = sp.hstack([StR, St])
#        if plot:
#            plt.figure()
#            sign = sp.zeros([self.n_z, self.nloops])
#            for i in xrange(self.n_z):
#                for j in xrange(self.nloops+1):
#                    if StR[i, j] <= 0:
#                        sign[i, j] = 1
#                    StR[i, j] = sp.log10(StR[i, j])
#
#            for j in xrange(self.n_z):
#                plt.plot(sp.array(StR[j, :]).flatten())

        Y_z = yt.copy()
        # Dz already defined
        S_z = St.copy()
#        if plot=='verbose':
#            print 'EP observations:'
#            print 'x '+str(self.X_z)
#            print 'y '+str(Y_z)
#            print 's '+str(S_z)
#            print 'de '+str(self.D_z)
        [Xep, Yep, Sep, Dep] = GPd.catObs([[self.X_c, self.Y_c, self.S_c, self.D_c], [self.X_z, Y_z, S_z, self.D_z]])
        self.gep = GPd.GPcore(Xep, Yep, Sep, Dep, self.kf)

        self.invalidflag = False
        return

    def infer_m(self, X_i, D_i):
        if self.invalidflag:
            self.runEP()
        return self.gep.infer_m(X_i, D_i)

    def infer_full(self, X_i, D_i):
        # print self.invalidflag
        if self.invalidflag:
            self.runEP()
        # print self.invalidflag
        return self.gep.infer_full(X_i, D_i)

    def infer_diag(self, X_i, D_i):
        if self.invalidflag:
            self.runEP()
        return self.gep.infer_diag(X_i, D_i)
