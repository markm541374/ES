import scipy as sp
from matplotlib import pyplot as plt
import scipy.stats as sps
import scipy.linalg as spl
from scipy.stats import norm as norm
import GPep
import GPd
# from scipy.optimize import minimize as mnz
import DIRECT
from functools import partial
# import time
from multiprocessing import pool
from tools import *
import sys

def searchEI(D, kf, lower, upper):
    Xo = D[0]
    Yo = D[1]
    So = D[2]
    Do = D[3]
    best = sp.amin(Yo)
    g1 = GPd.GPcore(Xo, Yo, So, Do, kf)
    ee = lambda post, best: -EI(best, post[0][0, 0], sp.sqrt(post[1][0, 0]))
    ef = lambda g, x, best: ee(g.infer_diag(sp.matrix(x), [[sp.NaN]]), best)
    eg = lambda x, y: (ef(g1, x, best), 0)
    [x, EImin, ierror] = DIRECT.solve(eg, lower, upper, user_data=[], algmethod=1, maxf=400)
    return x[0]


def searchEnt(D, kf, lower, upper):
    Xo = D[0]
    Yo = D[1]
    So = D[2]
    Do = D[3]
    Si = [0.0001]
    GG = makedraws([Xo, Yo, So, Do], kf, nd=100)
    ee = lambda x: -inferHmulti(GG, sp.matrix(sp.matrix(x)).T, Si)
    ef = lambda x, y: (ee(x), 0)
    [x, EImin, ierror] = DIRECT.solve(ef, lower, upper, user_data=[], algmethod=1, maxf=400)

    return x[0]


def EI_1dfixk(D, kf, xs, best):
    Xo = D[0]
    Yo = D[1]
    So = D[2]
    Do = D[3]
    g1 = GPd.GPcore(Xo, Yo, So, Do, kf)
    nh = len(xs)
    X_h = sp.matrix(sp.linspace(-1, 1, nh)).T
    D_h = [[sp.NaN]]*nh
    mh, vh = g1.infer_diag(X_h, D_h)
    EIi = sp.zeros(nh)
    for i in xrange(nh):
        EIi[i] = EI(best, mh[i, 0], sp.sqrt(vh[i, 0]))
    return EIi


def predict_1dfixk(D, kf, xs, ss, nd=400, nx_inner=101):
    Xo = D[0]
    Yo = D[1]
    So = D[2]
    Do = D[3]
    g1 = GPd.GPcore(Xo, Yo, So, Do, kf)

    nh = len(xs)
    ns = len(ss)

    X_h = sp.matrix(sp.linspace(-1, 1, nh)).T
    X_x = sp.matrix(sp.linspace(-1, 1, nx_inner)).T

    D_x = [[sp.NaN]]*nx_inner
    mh, Vh = g1.infer_full(X_x, D_x)
    Vh_cho = spl.cholesky(Vh, lower=True)
    H = sp.zeros([nh, ns])
    # ss=sp.logspace(4,-4,ns)

    for i in xrange(nd):
        if i % 25 == 0:
            print i
        dr = mh+Vh_cho*sp.matrix(sp.random.normal(size=nx_inner)).T
        xsi = dr.argmin()
        xs = X_x[xsi, :][0, 0]
    # eq constraint

        Xg = sp.matrix([[xs]])
        Yg = sp.matrix([[0.]])
        Sg = sp.matrix([[0.000000001]])
        Dg = [[0]]

        [Xc, Yc, Sc, Dc] = GPep.catObs([[Xo, Yo, So, Do], [Xg, Yg, Sg, Dg]])

    # ineq constraints

        Xz = sp.matrix([[xs], [xs]])

        Dz = [[sp.NaN], [0, 0]]
    # the inequality
        Iz = sp.matrix([[Yo[Yo.argmin(), :][0, 0]], [0.]])
    # sign of the inequality
        Gz = [0, 0]
        Nz = [0.1, 0.]

        g = GPep.GPcore(Xc, Yc, Sc, Dc, Xz, Dz, Iz, Gz, Nz, kf)

        for j in xrange(nh):
            Xt = sp.matrix([[X_h[j, 0]], [xs]])
            Dt = [[sp.NaN], [sp.NaN]]

            m, V = g.infer_full(Xt, Dt)
            s = V[0, 0] + V[1, 1] - V[0, 1] - V[1, 0]
            mu = m[1, 0] - m[0, 0]
            alpha = mu / sp.sqrt(s)
            beta = sps.norm.pdf(alpha) / sps.norm.cdf(alpha)

            vnxxs = V[0, 0] - beta * (beta + alpha) * (1. / s) * (V[0, 0] - V[0, 1])**2
            for k in xrange(ns):
                Hydxxs = 0.5 * sp.log(2 * sp.pi * sp.e * (vnxxs + ss[k]))
                H[j, k] += Hydxxs
    H = H / float(nd)

    for i in xrange(nh):
        m, v = g1.infer_diag(sp.matrix([[X_h[i, 0]]]), [[sp.NaN]])

        for k in xrange(ns):
            Hydx = 0.5 * sp.log(2 * sp.pi * sp.e * (v[0, 0]+ss[k]))
            H[i, k] = Hydx - H[i, k]
    return H


def EI(ER, mu, sigma):
    alpha = (-ER + mu) / sigma

    Z = norm.cdf(-alpha)

    if Z == 0.0:
        return sp.matrix(0.0)

    E = -mu + norm.pdf(alpha) * sigma / Z + ER
    EI = Z * E
    if sp.isfinite(EI):
        return sp.matrix(EI)[0.0]
    else:
        return sp.matrix(0.0)[0, 0]


def makedraws(D, kf, nd=400, nx_inner=101, mode='uniform', fig=False):
    Xo = D[0]
    Yo = D[1]
    So = D[2]
    Do = D[3]
    g1 = GPd.GPcore(Xo, Yo, So, Do, kf)
    g1.infer_full(sp.matrix([0]),[[sp.NaN]])
    # nh=len(xs)
    # ns=len(ss)
    # X_h = sp.matrix(sp.linspace(-1,1,nh)).T
    
    if mode == 'uniform':
        X_x = sp.matrix(sp.linspace(-1, 1, nx_inner)).T
    elif mode == 'slice':
        nacc = 0
        X_tmp = []
        mn = min(Yo)
        mx = max(Yo)
        while nacc<nx_inner:
            X_prop = sp.matrix(sp.random.uniform(-1, 1)).T
            Y_prop = g1.infer_m(X_prop, [[sp.NaN]])
            # [m_, v_] = g1.infer_diag(X_prop, [[sp.NaN]])
            # print m_
            # print v_
            # theta = -(m_+sp.sqrt(v_)-mx)/(mx-mn)
            theta = -(Y_prop[0,0]-mx)/(mx-mn)
            p = norm.cdf(2*theta-1.)
            
            if sp.random.uniform(0,1)<=p:
                nacc+=1
                X_tmp.append(X_prop)
                #print X_prop
                #print nacc
        X_x  = sp.vstack(X_tmp)
        # raise ValueError
    else:
        print 'invalid mode in makedraws'
        print mode
        raise ValueError

    D_x = [[sp.NaN]] * nx_inner
    mh, Vh = g1.infer_full(X_x, D_x)
    Vh_cho = spl.cholesky(Vh, lower=True)

    # ss=sp.logspace(4,-4,ns)
    allG = []
    for i in xrange(nd):
        # if i%25==0:
        #    print i
        dr = mh+Vh_cho*sp.matrix(sp.random.normal(size=nx_inner)).T
        if not fig == False:
            fig.plot(sp.array(X_x).flatten(), sp.array(dr).flatten(),'.')
        xsi = dr.argmin()
        xs = X_x[xsi, :][0, 0]
        # print xs
    # eq constraint

        Xg = sp.matrix([[xs]])
        Yg = sp.matrix([[0.]])
        Sg = sp.matrix([[0.000000001]])
        Dg = [[0]]

        [Xc, Yc, Sc, Dc] = GPep.catObs([[Xo, Yo, So, Do], [Xg, Yg, Sg, Dg]])
# ineq constraints

        Xz = sp.matrix([[xs], [xs]])
        Dz = [[sp.NaN], [0, 0]]
    # the inequality
        Iz = sp.matrix([[Yo[Yo.argmin(), :][0, 0]], [0.]])
    # sign of the inequality
        Gz = [0, 0]
        Nz = [0.1, 0.]

        g = GPep.GPcore(Xc, Yc, Sc, Dc, Xz, Dz, Iz, Gz, Nz, kf)
        try:
            g.runEP()
        except:
            
            
            import pickle
            object = [[Xo, Yo, So, Do], [Xg, Yg, Sg, Dg],[Xz, Dz, Iz, Gz, Nz],kf]
            file_pi = open('GPEPfail.obj', 'wb') 
            pickle.dump(object, file_pi)
            raise ValueError()
        allG.append([g, xs])
    return [[g1], allG]


def singleG(X_h, ss, G):
    g = G[0]
    xs = G[1]
    nh = len(sp.array(X_h).flatten())

    ns = len(ss)
    H = sp.zeros([nh, ns])

    # VV=[]
    # Vr=[]
    # C=[]
    # Al = []
    # Be = []
    # S =[]
    # Mu =[]
    # Sq=[]
    for j in xrange(nh):

        Xt = sp.matrix([[X_h[j, 0]], [xs]])
        Dt = [[sp.NaN], [sp.NaN]]

        m, V = g.infer_full(Xt, Dt)

        s = V[0, 0]+V[1, 1]-V[0, 1]-V[1, 0]
        
        
        if s<=10**-10:
            print 's<10**10'
            print s
            
        mu = m[1, 0]-m[0, 0] 
        alpha = - mu / sp.sqrt(s) # sign difference compared to the paper because I am minimizing
        
        beta = sp.exp(sps.norm.logpdf(alpha) - sps.norm.logcdf(alpha))
        coef = beta*(beta+alpha)*(1./s)*(V[0, 0]-V[0, 1])**2
        # C.append(coef)
        # S.append(s)
        # Al.append(alpha)
        # try:
        #    Be.append(sp.log10(beta))
        # except:
        #     print '!!!'
        #     print beta
        # Mu.append(mu)
        # Sq.append((V[0, 0]-V[0, 1])**2)
        
        vnxxs = V[0, 0]-coef
        # VV.append(vnxxs)
        # Vr.append(V[0,0])
        for k in xrange(ns):
            Hydxxs = 0.5*sp.log(2*sp.pi*sp.e*(vnxxs+ss[k]))
            if sp.iscomplex(Hydxxs):
                import pickle 
                object = [X_h, ss, G]
                file_pi = open('filename_pi.obj', 'w') 
                pickle.dump(object, file_pi)
                
            H[j, k] += Hydxxs
#    if plot:
#        print 'xxx'
#        plt.figure()
#        plt.plot(sp.linspace(-1,1,nh),VV)
#        plt.plot(xs,0,'ro')
#        plt.ylabel('vnxxs')
#        plt.savefig('foo0.png')
#        plt.figure()
#        plt.plot(sp.linspace(-1,1,nh),Vr)
#        plt.ylabel('V00')
#        plt.savefig('foo1.png')
#        plt.figure()
#        plt.plot(sp.linspace(-1,1,nh),C)
#        plt.ylabel('coef')
#        plt.savefig('foo3.png')
#        plt.figure()
#        plt.semilogy(sp.linspace(-1,1,nh),S)
#        plt.ylabel('s')
#        plt.savefig('foo4.png')
#        plt.figure()
#        plt.semilogy(sp.linspace(-1,1,nh),Al)
#        plt.ylabel('alpha')
#        plt.savefig('foo5.png')
#        try:
#            plt.figure()
#            plt.plot(sp.linspace(-1,1,nh),Be)
#            plt.ylabel('log10beta')
#            plt.savefig('foo6.png')
#        except:
#            print 'betaplot fail'
#            
#        plt.figure()
#        plt.plot(sp.linspace(-1,1,nh),Mu)
#        plt.ylabel('mu')
#        plt.savefig('foo7.png')
#        plt.figure()
#        plt.plot(sp.linspace(-1,1,nh),Sq)
#        plt.ylabel('sq')
#        plt.savefig('foo8.png')
    return H


def inferHmulti(G, X_h, ss, plot=False):
    nh = len(sp.array(X_h).flatten())
    ns = len(ss)
    Hp = partial(singleG, X_h, ss)
    g1 = G[0]
    p = pool.Pool(8)
    allh = p.map(Hp, G[1])
    ng = len(G[1])
    p.close()
    H = sp.zeros([nh, ns])
    for h in allh:
        for j in xrange(nh):
            for k in xrange(ns):
                H[j, k] += h[j, k]
    H = -H/float(ng)
    if plot:
        plt.figure()
        for k in xrange(ns):
            plt.plot(H[:, k])
        
    ing1 = 1./float(len(g1))
    # Vnx = [] #!!!!!
    # Htmp = sp.zeros([nh, ns])
    for gi in g1:

        for i in xrange(nh):
            m, v = gi.infer_diag(sp.matrix([[X_h[i, 0]]]), [[sp.NaN]])
            # Vnx.append(v[0,0]) # !!!!!!!!!!
            for k in xrange(ns):
                Hydx = 0.5*sp.log(2*sp.pi*sp.e*(v[0, 0]+ss[k]))
                H[i, k] += ing1*Hydx
                # Htmp[i, k] += ing1*Hydx
    # !!!!!!!!!!11
#    plt.figure()
#    plt.plot(Vnx)
#    plt.ylabel('vnx')
#    plt.savefig('foo2.png')
#    if plot:
#        plt.figure()
#        for k in xrange(ns):
#            plt.plot(Htmp[:, k])
#            
#    if plot:
#        plt.figure()
#        for k in xrange(ns):
#            plt.plot(H[:, k],'g')
    return H


def inferH(G, X_h, ss):

    nh = len(sp.array(X_h).flatten())
    ns = len(ss)
    H = sp.zeros([nh, ns])
    g1 = G[0]
    for Gi in G[1]:
        g = Gi[0]
        xs = Gi[1]
        for j in xrange(nh):

            Xt = sp.matrix([[X_h[j, 0]], [xs]])
            Dt = [[sp.NaN], [sp.NaN]]

            m, V = g.infer_full(Xt, Dt)
            # t1 = time.time()
            s = V[0, 0]+V[1, 1]-V[0, 1]-V[1, 0]
            mu = m[1, 0]-m[0, 0]
            alpha = mu/sp.sqrt(s)
            beta = sps.norm.pdf(alpha)/sps.norm.cdf(alpha)

            vnxxs = V[0, 0]-beta*(beta+alpha)*(1./s)*(V[0, 0]-V[0, 1])**2

            for k in xrange(ns):
                Hydxxs = 0.5*sp.log(2*sp.pi*sp.e*(vnxxs+ss[k]))
                H[j, k] += Hydxxs

    H = H/float(len(G[1]))

    for i in xrange(nh):
        m, v = g1.infer_diag(sp.matrix([[X_h[i, 0]]]), [[sp.NaN]])

        for k in xrange(ns):
            Hydx = 0.5*sp.log(2*sp.pi*sp.e*(v[0, 0]+ss[k]))
            H[i, k] = Hydx-H[i, k]
    return H


def slice_sample(dist, init, iters, sigma, step_out=True):
    """
    from http://isaacslavitt.com/2013/12/30/metropolis-hastings-and-
    slice-sampling/
    with some changes
    based on http://homepages.inf.ed.ac.uk/imurray2/teaching/09mlss/
    """

    # dist = joint_dist()

    # set up empty sample holder
    D = len(init)
    samples = sp.zeros((D, iters))

    # initialize
    xx = init.copy()

    for i in xrange(iters):
        sys.stdout.write('\r'+str(i))
        sys.stdout.flush()

        perm = range(D)
        sp.random.shuffle(perm)
        last_llh = dist.loglike(xx)

        for d in perm:
            llh0 = last_llh + sp.log(sp.random.rand())
            rr = sp.random.rand(1)
            x_l = xx.copy()
            x_l[d] = x_l[d] - rr * sigma[d]
            x_r = xx.copy()
            x_r[d] = x_r[d] + (1 - rr) * sigma[d]

            if step_out:
                llh_l = dist.loglike(x_l)
                while llh_l > llh0:
                    # print x_l
                    x_l[d] = x_l[d] - sigma[d]
                    llh_l = dist.loglike(x_l)
                llh_r = dist.loglike(x_r)
                while llh_r > llh0:
                    x_r[d] = x_r[d] + sigma[d]
                    llh_r = dist.loglike(x_r)

            x_cur = xx.copy()
            while True:
                xd = sp.random.rand() * (x_r[d] - x_l[d]) + x_l[d]
                x_cur[d] = xd.copy()
                last_llh = dist.loglike(x_cur)
                if last_llh > llh0:
                    xx[d] = xd.copy()
                    break
                elif xd > xx[d]:
                    x_r[d] = xd
                elif xd < xx[d]:
                    x_l[d] = xd
                else:
                    raise RuntimeError('Slice sampler shrank too far.')

        #if i % 1000 == 0:
        #    print 'iteration', i

        samples[:, i] = xx.copy().ravel()
    sys.stdout.write('\r                   \n')
    sys.stdout.flush()
    return samples



class EntPredictor():
    def __init__(self, D, lower, upper, kfgen, kfprior):
        self.D = D
        # function takes hyperparameters and returns a kernel fn
        self.kfgen = kfgen
        # function takes hyperparameters and returns ln(prior)
        self.kfprior = kfprior
        # g1 = GPd.GPcore(self.Xo, self.Yo, self.So, self.Do, self.kf)
        self.MLEflag = True
        self.HypSampleflag = True
        self.Predictorflag = True
        
        self.nx_inner = 101
        self.HYPsamplen = 100
        self.MLEsearchn = 800
        self.ENTsearchn = 500
        self.dmode = 'uniform'
        return

    def seekMLEkf(self):
        # searches for the MLE hyperparameters over +- 4 deacdes from unity
        print 'seeking MLE phyperparameters'

        self.dist = dist_obj(squaresllk, self.D, self.kfprior, self.kfgen)
        ee = lambda hyp, y: (-self.dist.loglike(hyp), 0)
        [xmintrue, miny, ierror] = DIRECT.solve(ee, [-4, -4], [4, 4], user_data=[], algmethod=1, maxf=self.MLEsearchn)
        print 'MLEhyperparameters: '+str([10**i for i in xmintrue])
        self.loghypMLE = xmintrue
        self.kfMLE = self.kfgen([10**i for i in xmintrue])
        self.gMLE = GPd.GPcore(self.D[0], self.D[1], self.D[2], self.D[3], self.kfMLE)
        self.MLEflag = False
        return

    def drawHypSamples(self, n, plot=False):

        if self.MLEflag:
            self.seekMLEkf()
        print 'Drawing '+str(n)+' hyperparameter samples'
        w_0 = self.loghypMLE
        sigma = 0.1*sp.ones(2)
        samples = slice_sample(self.dist, w_0, iters=n, sigma=sigma)
        self.kfSam = []
        for i in xrange(n):
            self.kfSam.append(self.kfgen([10**samples[0][i], 10**samples[1][i]]))

        if plot:
            plt.figure(figsize=(6,6))
            plt.xlabel('outputscale')
            plt.ylabel('lengthscale')
            for i in xrange(n):
                plt.loglog(10**samples[0][i], 10**samples[1][i], 'rx')
        self.HypSampleflag = False
        sys.stdout.write('\r              \n')
        sys.stdout.flush()
        return

    def initPredictor(self):
        if self.HypSampleflag:
            raise ValueError('no samples have been taken over the hyperparameters')
        print 'cov decomposition for inference over hyp samples'
        self.Pred = [[], []]
        # d = plt.figure(figsize=(8,8))
        # drawplt = d.add_subplot(111)
        
        tmp=[]
        for i, k in enumerate(self.kfSam):
            sys.stdout.write('\r'+str(i))
            sys.stdout.flush()
            
            g = makedraws(self.D, k, nd=1, nx_inner=self.nx_inner,mode=self.dmode, fig=False)
            self.Pred[0].append(g[0][0])
            self.Pred[1].append(g[1][0])
            # print g
            tmp.append(g[1][0][1])
            
            # except:
            #    print 'not using kf '+str(i)+' hyp: '+str(k)
        # print tmp
        plt.figure(figsize=(6, 6))
        # print tmp
        plt.hist(tmp, 20)
        # plt.axis([-1,1,0,len(tmp)])
        self.Predictorflag = False
        sys.stdout.write('\r           \n')
        sys.stdout.flush()
        return

    def predictGraph(self,X,S):
        if self.Predictorflag:
            raise ValueError('no predictor has been set up')
        n = len(X)
        print 'plotting entropy graph over ' + str(n) + ' points'

        H = [[]]*n
        for i in xrange(n):
            sys.stdout.write('\r'+str(i))
            sys.stdout.flush()
            H[i] = inferHmulti(self.Pred, sp.matrix(X[i]).T, S, plot=False)
        sys.stdout.write('\r             \n')
        sys.stdout.flush()
        
        return H

    def searchAtS(self, lower, upper, S):
        if self.HypSampleflag:
            self.drawHypSamples(self.HYPsamplen, plot=True)
        if self.Predictorflag:
            self.initPredictor()
        searchmax = self.ENTsearchn
        print 'searching over '+str(searchmax)+' iterations for maxEnt'
        global sc
        global me
        sc = 0
        me = 0

        def ee(x, y):
            global sc
            global me
            sys.stdout.write('\r'+str(sc)+' max found '+str(me))
            sys.stdout.flush()
            sc += 1
            Ent = inferHmulti(self.Pred, sp.matrix(x).T, S)
            if Ent > me:
                me = Ent
            return (-Ent, 0)

        [xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=[], algmethod=1, maxf=searchmax)
        del sc
        del me
        sys.stdout.write('\rMaxEnt at '+str(xmintrue)+'             ')
        return xmintrue

    def showEntGraph(self, Xi, S):
        nx = len(Xi)
        ns = len(S)
        if self.HypSampleflag:
            self.drawHypSamples(self.HYPsamplen, plot=True)
        if self.Predictorflag:
            self.initPredictor()
        H = self.predictGraph(Xi, S)
        Hplot = sp.zeros([nx, ns])
        for i in xrange(nx):
            for j in xrange(ns):
                Hplot[i, j] = H[i][0, j]

        plt.figure(figsize=(16, 16))
        a = GPd.plot1(self.gMLE, [-1], [1])
        for i in xrange(ns):
            a.plot(Xi, Hplot[:, i].flatten())
        return a
