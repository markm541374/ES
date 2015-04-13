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



def makedraws(D, kf, nd=400, nx_inner=101, mode='uniform', fig=False, plot='none'):
    Xo = D[0]
    Yo = D[1]
    So = D[2]
    Do = D[3]
    g1 = GPd.GPcore(Xo, Yo, So, Do, kf)
    g1.infer_full(sp.matrix([0]),[[sp.NaN]])
    # nh=len(xs)
    # ns=len(ss)
    # X_h = sp.matrix(sp.linspace(-1,1,nh)).T
    # print 'x0'
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
        #print 'x1'
    else:
        print 'invalid mode in makedraws'
        print mode
        raise ValueError
    # print 'x1'
    D_x = [[sp.NaN]] * nx_inner
    mh, Vh = g1.infer_full(X_x, D_x)
    Vh_cho = spl.cholesky(Vh, lower=True)
    # print 'x2'
    # ss=sp.logspace(4,-4,ns)
    allG = []
    for i in xrange(nd):
        # if i%25==0:
        #    print i
        dr = mh+Vh_cho*sp.matrix(sp.random.normal(size=nx_inner)).T
#        if not fig == False:
#            fig.plot(sp.array(X_x).flatten(), sp.array(dr).flatten(),'.')
        xsi = dr.argmin()
        xs = X_x[xsi, :][0, 0]
        # print xs
    # eq constraint

        Xg = sp.matrix([[xs]])
        Yg = sp.matrix([[0.]])
        Sg = sp.matrix([[0.001]])
        Dg = [[0]]

        [Xc, Yc, Sc, Dc] = GPep.catObs([[Xo, Yo, So, Do], [Xg, Yg, Sg, Dg]])
# ineq constraints

        Xz = sp.matrix([[xs], [xs]])
        Dz = [[sp.NaN], [0, 0]]
    # the inequality
        
        
        Iz = sp.matrix([[Yo[Yo.argmin(), :][0, 0]], [0.]])
    # sign of the inequality
        Gz = [0, 0]
        Nz = [So[Yo.argmin(), :], 0.] #!!!!!this value is important, should it e the sigma for hte min obs or the posterior at that piont??
        
        g = GPep.GPcore(Xc, Yc, Sc, Dc, Xz, Dz, Iz, Gz, Nz, kf)
        
        if plot=='verbose':
            gp0=GPd.GPcore(Xo, Yo, So, Do, kf)
            gp1=GPd.GPcore(Xc, Yc, Sc, Dc, kf)
            xt=sp.linspace(-1,1,500)
            [m0,v0] = gp0.infer_diag(sp.matrix(xt).T, [[sp.NaN]]*500)
            [m1,v1] = gp1.infer_diag(sp.matrix(xt).T, [[sp.NaN]]*500)
            plt.figure()
            plt.semilogy(xt, sp.array(v0).flatten())
            plt.semilogy(xt, sp.array(v1).flatten())
            plt.xlabel('with/without grad obs')
        
        
        try:
            #print 'x2'
            g.runEP(plot=plot)
            
        except:
            
            
            import pickle
            object = [[Xo, Yo, So, Do], [Xg, Yg, Sg, Dg],[Xz, Dz, Iz, Gz, Nz],kf]
            file_pi = open('GPEPfail.obj', 'wb') 
            pickle.dump(object, file_pi)
            raise ValueError()
        #if plot=='verbose':
        #    a = GPd.plot1(g.gep,[-1],[1])
            
        allG.append([g, xs])
    return [[g1], allG]


def singleG(X_h, ss, G, plot='none',obstype=[sp.NaN]):
    g = G[0]
    xs = G[1]
    nh = len(sp.array(X_h).flatten())

    ns = len(ss)
    H = sp.zeros([nh, ns])

    VV=[]
    Vr=[]
    C=[]
    Al = []
    Be = []
    S =[]
    Mu =[]
    Sq=[]
    for j in xrange(nh):

        Xt = sp.matrix([[X_h[j, 0]], [xs]])
        Dt = [obstype, [sp.NaN]]

        m, V = g.infer_full(Xt, Dt)

        s = V[0, 0]+V[1, 1]-V[0, 1]-V[1, 0]
        
        
        if s<=10**-10:
            print 's<10**10'
            print s
            
        mu = m[1, 0]-m[0, 0] 
        alpha = - mu / sp.sqrt(s) # sign difference compared to the paper because I am minimizing
        
        beta = sp.exp(sps.norm.logpdf(alpha) - sps.norm.logcdf(alpha))
        coef = beta*(beta+alpha)*(1./s)*(V[0, 0]-V[0, 1])**2
        C.append(coef)
        S.append(s)
        Al.append(alpha)
        Be.append(beta)
        Mu.append(mu)
        Sq.append((V[0, 0]-V[0, 1])**2)
        
        vnxxs = V[0, 0]-coef
        VV.append(vnxxs)
        Vr.append(V[0,0])
        for k in xrange(ns):
            Hydxxs = 0.5*sp.log(2*sp.pi*sp.e*(vnxxs+ss[k]))
            if sp.iscomplex(Hydxxs):
                import pickle 
                object = [X_h, ss, G]
                file_pi = open('filename_pi.obj', 'w') 
                pickle.dump(object, file_pi)
                
            H[j, k] += Hydxxs
#    if plot == 'verbose':
#        
#        plt.figure()
#        plt.semilogy(sp.linspace(-1,1,nh),VV)
#        plt.plot(xs,0,'ro')
#        plt.ylabel('vnxxs')
#        plt.savefig('foo0.png')
#        plt.figure()
#        plt.semilogy(sp.linspace(-1,1,nh),Vr)
#        plt.ylabel('V00')
#        plt.savefig('foo1.png')
#        
#        plt.figure()
#        try:            
#            plt.semilogy(sp.linspace(-1,1,nh),C)
#        except:
#            plt.plot(sp.linspace(-1,1,nh),C)
#        plt.ylabel('coef')
#        plt.savefig('foo3.png')
#        
#        plt.figure()
#        plt.semilogy(sp.linspace(-1,1,nh),S)
#        plt.ylabel('s')
#        plt.savefig('foo4.png')
#        plt.figure()
#        plt.semilogy(sp.linspace(-1,1,nh),Al)
#        plt.ylabel('alpha')
#        plt.savefig('foo5.png')
#        
#        plt.figure()
#        try:
#            plt.semilogy(sp.linspace(-1,1,nh),Be)
#        except:
#            plt.plot(sp.linspace(-1,1,nh),Be)
#        plt.ylabel('beta')
#        plt.savefig('foo6.png')
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


def inferHmulti(G, X_h, ss, plot='none', obstype=[sp.NaN]):
    
    nh = len(sp.array(X_h).flatten())
    ns = len(ss)
    if nh<2:
        plot='none'
    Hp = partial(singleG, X_h, ss, plot=plot)
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
    if plot == 'verbose':
        f=plt.figure(figsize=(8, 8))
        ax0=f.add_subplot(111)
        for k in xrange(ns):
            ax0.plot(sp.array(X_h).flatten(),-H[:, k],'b')
    
    ing1 = 1./float(len(g1))
    Vnx = [] #!!!!!
    Htmp = sp.zeros([nh, ns])
    for gi in g1:

        for i in xrange(nh):
            m, v = gi.infer_diag(sp.matrix([[X_h[i, 0]]]), [obstype])
            Vnx.append(v[0,0]) # !!!!!!!!!!
            for k in xrange(ns):
                Hydx = 0.5*sp.log(2*sp.pi*sp.e*(v[0, 0]+ss[k]))
                H[i, k] += ing1*Hydx
                Htmp[i, k] += ing1*Hydx
    
    
    if plot == 'verbose':
        for k in xrange(ns):
            ax0.plot(sp.array(X_h).flatten(), Htmp[:, k],'r')
            plt.ylabel('b-Hmin r-Hexist')
            plt.savefig('foo2.png')
#    if plot:
#        plt.figure()
#        for k in xrange(ns):
#            plt.plot(H[:, k],'g')
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
        
        self.lb=lower
        self.ub=upper
        # g1 = GPd.GPcore(self.Xo, self.Yo, self.So, self.Do, self.kf)
        self.MLEflag = True
        self.HypSampleflag = True
        self.Predictorflag = True
        
        self.nx_inner = 101
        self.HYPsamplen = 100
        self.MLEsearchn = 800
        self.ENTsearchn = 500
        self.dmode = 'slice'
        return

    def seekMLEkf(self):
        # searches for the MLE hyperparameters over +- 4 deacdes from unity
        print 'seeking MLE phyperparameters'

        self.dist = dist_obj(squaresllk, self.D, self.kfprior, self.kfgen)
        ee = lambda hyp, y: (-self.dist.loglike(hyp), 0)
        [xmintrue, miny, ierror] = DIRECT.solve(ee, [-3, -3], [5, 2], user_data=[], algmethod=1, maxf=self.MLEsearchn)
        print 'MLEhyperparameters: '+str([10**i for i in xmintrue])
        self.loghypMLE = xmintrue
        self.kfMLE = self.kfgen([10**i for i in xmintrue])
        self.gMLE = GPd.GPcore(self.D[0], self.D[1], self.D[2], self.D[3], self.kfMLE)
        self.MLEflag = False
        return

    def drawHypSamples(self, n, plot='none'):

        if self.MLEflag:
            self.seekMLEkf()
        print 'Drawing '+str(n)+' hyperparameter samples'
        w_0 = self.loghypMLE
        sigma = 0.05*sp.ones(2)
        samples = slice_sample(self.dist, w_0, iters=n, sigma=sigma)
        # print samples
        self.kfSam = []
        self.hySam = []
        for i in xrange(n):
            self.kfSam.append(self.kfgen([10**samples[0][i], 10**samples[1][i]]))
            self.hySam.append([10**samples[0][i], 10**samples[1][i]])
        if plot == 'verbose' or plot == 'basic':
            plt.figure(figsize=(6,6))
            plt.xlabel('outputscale')
            plt.ylabel('lengthscale')
            for i in xrange(n):
                plt.loglog(10**samples[0][i], 10**samples[1][i], 'rx')
        self.HypSampleflag = False
        sys.stdout.write('\r              \n')
        sys.stdout.flush()
        return

    def initPredictor(self, plot='none'):
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
            try:
                g = makedraws(self.D, k, nd=1, nx_inner=self.nx_inner,mode=self.dmode, fig=False, plot=plot)
                self.Pred[0].append(g[0][0])
                self.Pred[1].append(g[1][0])
            # print g
                tmp.append(g[1][0][1])
            
            except:
                print 'not using kf '+str(i)+' hyp: '+str(self.hySam[i])
                
        
        if plot == 'verbose' or plot == 'basic':
            plt.figure(figsize=(6, 6))
            try:
                plt.hist(tmp, 20)
            except:
                print 'histogram error'
                print tmp
        # plt.axis([-1,1,0,len(tmp)])
        self.Predictorflag = False
        sys.stdout.write('\r           \n')
        sys.stdout.flush()
        return

    def predictGraph(self,X,S, plot='none',obstype=[sp.NaN]):
        if self.Predictorflag:
            raise ValueError('no predictor has been set up')
        n = len(X)
        print 'plotting entropy graph over ' + str(n) + ' points'

        H = [[]]*n
        for i in xrange(n):
            sys.stdout.write('\r'+str(i))
            sys.stdout.flush()
            H[i] = inferHmulti(self.Pred, sp.matrix(X[i]).T, S, plot=plot, obstype=obstype)
        sys.stdout.write('\r             \n')
        sys.stdout.flush()
        
        return H

    def searchAtS(self, lower, upper, S, plot='none', obstype=[sp.NaN]):
        if self.HypSampleflag:
            self.drawHypSamples(self.HYPsamplen, plot=plot)
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
            Ent = inferHmulti(self.Pred, sp.matrix(x).T, S, obstype=obstype)
            if Ent > me:
                me = Ent
            return (-Ent, 0)

        [xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=[], algmethod=1, maxf=searchmax)
        del sc
        del me
        sys.stdout.write('\rMaxEnt at '+str(xmintrue)+'             ')
        return xmintrue

    def searchOverS(self, lower, upper, S, U, plot='none', obstype=[sp.NaN]):
        if self.HypSampleflag:
            self.drawHypSamples(self.HYPsamplen, plot=plot)
        if self.Predictorflag:
            self.initPredictor()
        searchmax = self.ENTsearchn
        print 'searching over '+str(searchmax)+' iterations for maxEnt'
        global sc
        global me
        sc = 0
        me = 0
        ns=len(S)
        def ee(x, y):
            global sc
            global me
            sys.stdout.write('\r'+str(sc)+' max found '+str(me))
            sys.stdout.flush()
            sc += 1
            Ent = inferHmulti(self.Pred, sp.matrix(x).T, S, obstype=obstype)
            EU=sp.zeros(ns)
            for i in xrange(ns):
                #print Ent[0]
                EU[i]=(Ent[0][i])/float(U[i])
            EUmax=max(EU)
            if EUmax > me:
                me = EUmax
            return (-EUmax, 0)

        [xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=U, algmethod=1, maxf=searchmax)
        del sc
        del me
        
        Ent = inferHmulti(self.Pred, sp.matrix(xmintrue).T, S, obstype=obstype)
        EU=sp.zeros(ns)
        for i in xrange(ns):
            EU[i]=Ent[0][i]/float(U[i])
        j=sp.argmax(EU)
        
        sys.stdout.write('\rMaxEU at '+str(xmintrue)+'with s= '+str(S[j]))
        return [xmintrue,j]
        
    def showEntGraph(self, Xi, S,U='none', plot='basic', obstype=[sp.NaN]):
        if U=='none':
            uflag=False
            U=sp.ones(len(S))
        else:
            uflag=True
        nx = len(Xi)
        ns = len(S)
        if self.HypSampleflag:
            self.drawHypSamples(self.HYPsamplen, plot=plot)
        if self.Predictorflag:
            self.initPredictor(plot=plot)
        H = self.predictGraph(Xi, S, obstype=obstype)
        Hplot = sp.zeros([nx, ns])
        for i in xrange(nx):
            for j in xrange(ns):
                Hplot[i, j] = (H[i][0, j])/float(U[j])
        
        plt.figure(figsize=(16, 16))
        a = GPd.plot1(self.gMLE, [-1], [1])
        a2=a.twinx()
        for i in xrange(ns):
            try:
                a2.semilogy(Xi, Hplot[:, i].flatten())
            except:
                print 'failed to plot entropy for s= '+str(S[i])
                print Hplot[:, i].flatten()
                
            
        if uflag:
            uopt=sp.zeros(nx)
            for i in xrange(nx):
                j=sp.argmax(Hplot[i,:].flatten())
                uopt[i]=U[j]
            plt.figure()
            #print uopt
            plt.semilogy(Xi,uopt)
        return a

    def inferPostGP(self, x, d, plot='none'):
        if self.HypSampleflag:
            self.drawHypSamples(self.HYPsamplen, plot=plot)
        if self.Predictorflag:
            self.initPredictor(plot=plot)
        np=len(d)
        acc0=sp.zeros(np)
        acc1=sp.zeros(np)
        vp=sp.zeros(np)
        ng=len(self.Pred[0])
        for g in self.Pred[0]:
            [m,v]=g.infer_diag(sp.matrix(x).T,d)
            for i in xrange(np):
                acc0[i]+=1./v[i]
                acc1[i]+=m[i]/v[i]
                vp[i] += v[i]
        mp=sp.zeros(np)
        
        for i in xrange(np):
            mp[i]=acc1[i]/(acc0[i])
            vp[i]=vp[i]/float(ng)
        
        return [mp,vp]
    
    def searchYminEst(self):
        def ee(x,y):
            [m,v] = self.inferPostGP(x,[[sp.NaN]])
            return (m,0)
            
        [xmintrue, miny, ierror] = DIRECT.solve(ee, self.lb, self.ub, user_data=[], algmethod=1, maxf=800)
        return [xmintrue,miny]
        
    def plotPostGP(self, x, d, plotEI=False):
        [m, v] = self.inferPostGP(x, d)
        np = len(d)
        ub = sp.zeros(np)
        lb = sp.zeros(np)
        for i in xrange(np):
            ub[i] = m[i]+2*sp.sqrt(v[i])
            lb[i] = m[i]-2*sp.sqrt(v[i])

        f0 = plt.figure(figsize=(8, 8))
        a0 = f0.add_subplot(111)
        a0.plot(x, m)
        a0.fill_between(x, lb, ub, facecolor='lightskyblue', alpha=0.5)

        xs = sp.array(self.D[0]).flatten()
        ys = sp.array(self.D[1]).flatten()
        for i in xrange(len(self.D[3])):
            if (sp.isnan(self.D[3][i])).any():
                a0.plot(xs, ys, 'rx')
                
        if plotEI:
            best = min(self.D[1])[0,0]
            ei = sp.zeros(np)
            for i in xrange(np):
                ei[i]=EI(best,m[i],sp.sqrt(v[i]))
            a1=a0.twinx()
            a1.plot(x,ei,'r')
        return a0

    def plotMLEGP(self):
        if self.MLEflag:
            self.seekMLEkf()
        a = GPd.plot1(self.gMLE, [-1], [1])
        return
        
    def searchEI(self, lower, upper, plot='none'):
        searchmax=400
        best = min(self.D[1])[0,0]
        
        global sc
        global me
        sc = 0
        me = 0
        print 'searching over '+str(searchmax)+' iterations for maxEI'
        def ee(x, y):
            [m,v] = self.inferPostGP(x,[[sp.NaN]])
            ei = EI(best,m[0],sp.sqrt(v[0]))
            
            global sc
            global me
            sc+=1
            if ei[0,0]>me:
                me=ei[0,0]
            sys.stdout.write('\r'+str(sc)+' max found '+str(me)+'      ')
            sys.stdout.flush()
            
            return (-ei, 0)

        [xmintrue, miny, ierror] = DIRECT.solve(ee, lower, upper, user_data=[], algmethod=1, maxf=searchmax)
        del sc
        del me
        sys.stdout.write('\rMaxEI at '+str(xmintrue)+'             ')
        
        return xmintrue
        
class Optimizer():
    def __init__(self, f, kfGen, kfPrior, lb, ub):
        self.f = f
        self.kfGen = kfGen
        self.kfPrior = kfPrior
        self.lb = lb
        self.ub = ub
        
        self.MLEsearchn = 200
        self.HYPsamplen = 30
        self.ENTsearchn = 160
        self.result=[]
        self.Uo=[]
        return
        
    def initrandobs(self, n_init, s_init):
        x = sp.random.uniform(self.lb, self.ub, n_init)
        y = map(self.f, x)+sp.random.normal(scale=sp.sqrt(s_init), size=n_init)
        self.Xo = sp.matrix(x).T
        self.Yo = sp.matrix(y).T
        self.So = sp.matrix([[s_init]]*n_init)
        self.Do = [[sp.NaN]]*n_init
        
        return
        
    def initspecobs(self, Xo, Yo, So, Do):
        self.Xo = Xo
        self.Yo = Yo
        self.So = So
        self.Do = Do
        
        return
        
    def searchats(self, s, obstype=[sp.NaN], method='Ent'):
        e = EntPredictor([self.Xo, self.Yo, self.So, self.Do], self.lb, self.ub, self.kfGen, self.kfPrior)
        e.MLEsearchn = self.MLEsearchn
        e.HYPsamplen = self.HYPsamplen
        e.ENTsearchn = self.ENTsearchn
        
        n = 150
        Xi = sp.linspace(-1, 1, n)
        if method=='Ent':
            a = e.showEntGraph(Xi, [s], plot='basic', obstype=obstype)
            xmin = e.searchAtS(self.lb, self.ub, [s], obstype=obstype)
        if method=='EI':
            a = e.plotPostGP(Xi,[[sp.NaN]]*n, plotEI=True)
            xmin = e.searchEI(self.lb, self.ub)
        a.plot(xmin, [0], 'go')
        a.plot(xmin, self.f(xmin[0]), 'ro')
        xIn = xmin[0]
        yIn = self.f(xIn)+sp.random.normal(scale=sp.sqrt(s))

        self.Xo = sp.vstack([self.Xo, xIn])
        self.Yo = sp.vstack([self.Yo, yIn])
        self.So = sp.vstack([self.So, s])
        self.Do.append([sp.NaN])
        plt.show()
        print 'getting current min est'
        res = e.searchYminEst()
        print res
        self.result.append(res)
        return
        
    def searchOvers(self, s, u, obstype=[sp.NaN], method='Ent'):
        e = EntPredictor([self.Xo, self.Yo, self.So, self.Do], self.lb, self.ub, self.kfGen, self.kfPrior)
        e.MLEsearchn = self.MLEsearchn
        e.HYPsamplen = self.HYPsamplen
        e.ENTsearchn = self.ENTsearchn
        
        n = 150
        Xi = sp.linspace(-1, 1, n)
        a = e.showEntGraph(Xi, s,U=u, plot='basic', obstype=obstype)
        [xmin,j] = e.searchOverS(self.lb, self.ub, s,u, obstype=obstype)
        a.plot(xmin, [0], 'go')
        a.plot(xmin, self.f(xmin[0]), 'ro')
        xIn = xmin[0]
        yIn = self.f(xIn)+sp.random.normal(scale=sp.sqrt(s[j]))
        self.Uo.append(u[j])
        self.Xo = sp.vstack([self.Xo, xIn])
        self.Yo = sp.vstack([self.Yo, yIn])
        self.So = sp.vstack([self.So, s[j]])
        self.Do.append([sp.NaN])
        plt.show()
        print 'getting current min est'
        res = e.searchYminEst()
        print res
        self.result.append(res)
        return
        
    def showinov(self, s, obstype=[sp.NaN]):
        e = EntPredictor([self.Xo, self.Yo, self.So, self.Do], self.lb, self.ub, self.kfGen, self.kfPrior)
        e.MLEsearchn = self.MLEsearchn
        e.HYPsamplen = self.HYPsamplen
        e.ENTsearchn = self.ENTsearchn
        
        n = 300
        Xi = sp.linspace(-1, 1, n)
        a = e.showEntGraph(Xi, [s], plot='basic', obstype=obstype)
        return
        
