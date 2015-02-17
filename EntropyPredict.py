import scipy as sp
#from matplotlib import pyplot as plt
import scipy.stats as sps
import scipy.linalg as spl
from scipy.stats import norm as norm
import GPep
import GPd
from scipy.optimize import minimize as mnz
import DIRECT
from functools import partial
import time
from multiprocessing import pool

def searchEI(D,kf,lower,upper):
    Xo=D[0]
    Yo=D[1]
    So=D[2]
    Do=D[3]
    best=sp.amin(Yo)
    g1=GPd.GPcore(Xo,Yo,So,Do,kf)
    ee = lambda post,best : -EI(best,post[0][0,0],sp.sqrt(post[1][0,0]))
    ef = lambda g,x,best: ee(g.infer_diag(sp.matrix(x),[[sp.NaN]]),best)
    eg = lambda x,y:(ef(g1,x,best),0)
    [x,EImin,ierror]=DIRECT.solve(eg,lower,upper,user_data=[],algmethod=1,maxf=400)
    
    return x[0]

def searchEnt(D,kf,lower,upper):
    Xo=D[0]
    Yo=D[1]
    So=D[2]
    Do=D[3]
    Si=[0.01]
    GG=makedraws([Xo,Yo,So,Do],kf,nd=100)

    ee = lambda x:-inferHmulti(GG,sp.matrix(sp.matrix(x)).T,Si)
    ef = lambda x,y:(ee(x),0)
    [x,EImin,ierror]=DIRECT.solve(ef,lower,upper,user_data=[],algmethod=1,maxf=400)

    return x[0]

def EI_1dfixk(D,kf,xs,best):
    Xo=D[0]
    Yo=D[1]
    So=D[2]
    Do=D[3]
    g1=GPd.GPcore(Xo,Yo,So,Do,kf)
    nh=len(xs)
    X_h = sp.matrix(sp.linspace(-1,1,nh)).T
    D_h = [[sp.NaN]]*nh
    mh,vh=g1.infer_diag(X_h,D_h)
    EIi=sp.zeros(nh)
    for i in xrange(nh):
        EIi[i]=EI(best,mh[i,0],sp.sqrt(vh[i,0]))
    return EIi

def predict_1dfixk(D,kf,xs,ss,nd=400,nx_inner=101):
    Xo=D[0]
    Yo=D[1]
    So=D[2]
    Do=D[3]
    g1=GPd.GPcore(Xo,Yo,So,Do,kf)

    nh=len(xs)
    ns=len(ss)

    X_h = sp.matrix(sp.linspace(-1,1,nh)).T
    X_x = sp.matrix(sp.linspace(-1,1,nx_inner)).T

    D_x = [[sp.NaN]]*nx_inner
    mh,Vh=g1.infer_full(X_x,D_x)
    Vh_cho = spl.cholesky(Vh,lower=True)
    H=sp.zeros([nh,ns])
    #ss=sp.logspace(4,-4,ns)

    for i in xrange(nd):
        if i%25==0:
            print i
        dr = mh+Vh_cho*sp.matrix(sp.random.normal(size=nx_inner)).T
        xsi=dr.argmin()
        xs=X_x[xsi,:][0,0]
    #eq constraint
    
        Xg=sp.matrix([[xs]])
        Yg=sp.matrix([[0.]])
        Sg=sp.matrix([[0.000000001]])
        Dg=[[0]]
    
        [Xc,Yc,Sc,Dc]=GPep.catObs([[Xo,Yo,So,Do],[Xg,Yg,Sg,Dg]])
    
    #ineq constraints
    
        Xz=sp.matrix([[xs],[xs]])
    
        Dz=[[sp.NaN],[0,0]]
    #the inequality
        Iz=sp.matrix([[Yo[Yo.argmin(),:][0,0]],[0.]])
    #sign of the inequality
        Gz=[0,0]
        Nz=[0.1,0.]
    
        g=GPep.GPcore(Xc,Yc,Sc,Dc,Xz,Dz,Iz,Gz,Nz,kf)
    
        for j in xrange(nh):
            Xt=sp.matrix([[X_h[j,0]],[xs]])
            Dt=[[sp.NaN],[sp.NaN]]
        
            m,V=g.infer_full(Xt,Dt)
            s=V[0,0]+V[1,1]-V[0,1]-V[1,0]
            mu=m[1,0]-m[0,0]
            alpha=mu/sp.sqrt(s)
            beta=sps.norm.pdf(alpha)/sps.norm.cdf(alpha)
        
            vnxxs=V[0,0]-beta*(beta+alpha)*(1./s)*(V[0,0]-V[0,1])**2
            for k in xrange(ns):
                Hydxxs=0.5*sp.log(2*sp.pi*sp.e*(vnxxs+ss[k]))
                H[j,k]+=Hydxxs
    H=H/float(nd)
    
    for i in xrange(nh):
        m,v=g1.infer_diag(sp.matrix([[X_h[i,0]]]),[[sp.NaN]])
    
        for k in xrange(ns):
            Hydx=0.5*sp.log(2*sp.pi*sp.e*(v[0,0]+ss[k]))
            H[i,k]=Hydx-H[i,k]
    return H

def EI(ER,mu,sigma):
    alpha=(-ER+mu)/sigma

    Z=norm.cdf(-alpha)

    if Z==0.0:
        return sp.matrix(0.0)
#print "alpha: "+str(alpha)
#print "Z: "+str(Z)
    E=-mu+norm.pdf(alpha)*sigma/Z+ER
    EI=Z*E
    if sp.isfinite(EI):
        return sp.matrix(EI)[0.0]
    else:
        return sp.matrix(0.0)[0,0]

def makedraws(D,kf,nd=400,nx_inner=101):
    Xo=D[0]
    Yo=D[1]
    So=D[2]
    Do=D[3]
    g1=GPd.GPcore(Xo,Yo,So,Do,kf)

    #nh=len(xs)
    #ns=len(ss)

    #X_h = sp.matrix(sp.linspace(-1,1,nh)).T
    X_x = sp.matrix(sp.linspace(-1,1,nx_inner)).T

    D_x = [[sp.NaN]]*nx_inner
    mh,Vh=g1.infer_full(X_x,D_x)
    Vh_cho = spl.cholesky(Vh,lower=True)
    
    #ss=sp.logspace(4,-4,ns)
    allG=[]
    for i in xrange(nd):
        if i%25==0:
            print i
        dr = mh+Vh_cho*sp.matrix(sp.random.normal(size=nx_inner)).T
        xsi=dr.argmin()
        xs=X_x[xsi,:][0,0]
    #eq constraint
    
        Xg=sp.matrix([[xs]])
        Yg=sp.matrix([[0.]])
        Sg=sp.matrix([[0.000000001]])
        Dg=[[0]]
    
        [Xc,Yc,Sc,Dc]=GPep.catObs([[Xo,Yo,So,Do],[Xg,Yg,Sg,Dg]])
    
    #ineq constraints
    
        Xz=sp.matrix([[xs],[xs]])
    
        Dz=[[sp.NaN],[0,0]]
    #the inequality
        Iz=sp.matrix([[Yo[Yo.argmin(),:][0,0]],[0.]])
    #sign of the inequality
        Gz=[0,0]
        Nz=[0.1,0.]
    
        g=GPep.GPcore(Xc,Yc,Sc,Dc,Xz,Dz,Iz,Gz,Nz,kf)
        g.runEP()
        allG.append([g,xs])
    return [g1,allG]

def singleG(X_h,ss,G):
    g=G[0]
    xs=G[1]
    nh=len(sp.array(X_h).flatten())
    
    ns=len(ss)
    H=sp.zeros([nh,ns])
    
    for j in xrange(nh):
        
        Xt=sp.matrix([[X_h[j,0]],[xs]])
        Dt=[[sp.NaN],[sp.NaN]]
            
        m,V=g.infer_full(Xt,Dt)
        
        s=V[0,0]+V[1,1]-V[0,1]-V[1,0]
        mu=m[1,0]-m[0,0]
        alpha=mu/sp.sqrt(s)
        beta=sps.norm.pdf(alpha)/sps.norm.cdf(alpha)
        
        vnxxs=V[0,0]-beta*(beta+alpha)*(1./s)*(V[0,0]-V[0,1])**2
        
        for k in xrange(ns):
            Hydxxs=0.5*sp.log(2*sp.pi*sp.e*(vnxxs+ss[k]))
            H[j,k]+=Hydxxs
        
    return H
def inferHmulti(G,X_h,ss):
    nh=len(sp.array(X_h).flatten())
    ns=len(ss)
    Hp=partial(singleG,X_h,ss)
    g1=G[0]
    p=pool.Pool(8)
    allh=p.map(Hp,G[1])
    ng=len(G[1])
    p.close()
    H=sp.zeros([nh,ns])
    for h in allh:
        for j in xrange(nh):
            for k in xrange(ns):
                H[j,k]+=h[j,k]
        
    H=H/float(ng)
    
    for i in xrange(nh):
        m,v=g1.infer_diag(sp.matrix([[X_h[i,0]]]),[[sp.NaN]])
    
        for k in xrange(ns):
            Hydx=0.5*sp.log(2*sp.pi*sp.e*(v[0,0]+ss[k]))
            H[i,k]=Hydx-H[i,k]
    return H

def inferH(G,X_h,ss):
    
    nh=len(sp.array(X_h).flatten())
    
    ns=len(ss)
    H=sp.zeros([nh,ns])
    g1=G[0]
    for Gi in G[1]:
        g=Gi[0]
        
        xs=Gi[1]
        for j in xrange(nh):
            
            Xt=sp.matrix([[X_h[j,0]],[xs]])
            Dt=[[sp.NaN],[sp.NaN]]
            
            m,V=g.infer_full(Xt,Dt)
            t1=time.time()
            s=V[0,0]+V[1,1]-V[0,1]-V[1,0]
            mu=m[1,0]-m[0,0]
            alpha=mu/sp.sqrt(s)
            beta=sps.norm.pdf(alpha)/sps.norm.cdf(alpha)
        
            vnxxs=V[0,0]-beta*(beta+alpha)*(1./s)*(V[0,0]-V[0,1])**2
            
            for k in xrange(ns):
                Hydxxs=0.5*sp.log(2*sp.pi*sp.e*(vnxxs+ss[k]))
                H[j,k]+=Hydxxs
            
            
    H=H/float(len(G[1]))
    
    for i in xrange(nh):
        m,v=g1.infer_diag(sp.matrix([[X_h[i,0]]]),[[sp.NaN]])
    
        for k in xrange(ns):
            Hydx=0.5*sp.log(2*sp.pi*sp.e*(v[0,0]+ss[k]))
            H[i,k]=Hydx-H[i,k]
    return H
