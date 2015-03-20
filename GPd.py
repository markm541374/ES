import scipy as sp
from scipy import linalg as spl
from matplotlib import pyplot as plt
from numpy.linalg import slogdet as slogdet
# conbine sets of observations


def catObs(O):
    i = True
    for Os in O:
        if i:
            i = False
            X = Os[0]
            Y = Os[1]
            S = Os[2]
            D = Os[3]
        else:
            X = sp.vstack([X, Os[0]])
            Y = sp.vstack([Y, Os[1]])
            S = sp.vstack([S, Os[2]])
            D = D+Os[3]
    return [X, Y, S, D]

# given a 1d input return a matrix with the input on the diagonal


def vec2trace(v):
    (l, _) = v.shape
    M = sp.eye(l)
    for i in range(l):
        M[i, i] = v[i]
    return M


def plot1(g, llimit, ulimit):
    n_p = 120
    X_p = sp.matrix(sp.linspace(llimit[0], ulimit[0], n_p)).T   # plot points
    D_p = [[sp.NaN]]*n_p

    [m_p, V_p] = g.infer_diag(X_p, D_p)

    ub = sp.zeros(n_p)
    lb = sp.zeros(n_p)
    for i in xrange(n_p):
        ub[i] = m_p[i, 0]+2*sp.sqrt(V_p[i, 0])
        lb[i] = m_p[i, 0]-2*sp.sqrt(V_p[i, 0])
    f0 = plt.figure(figsize=(8, 8))
    a0 = f0.add_subplot(111)
    a0.plot(sp.array(X_p).flatten(), sp.array(m_p).flatten())
    a0.fill_between(sp.array(X_p).flatten(), lb, ub, facecolor='lightskyblue', alpha=0.5)
    # print g
    # print dir(g)
    xs = sp.array(g.X_s).flatten()
    ys = sp.array(g.Y_s).flatten()
    for i in xrange(len(g.D_s)):
        if sp.isnan(g.D_s[i]):
            a0.plot(xs, ys, 'rx')
    return a0


# return xlogx but return zero for x=0 instead of raising a log0 error
def xlx0(x):
    if x == 0:
        return 0
    else:
        return x*sp.log(x)


# sqexp kernel generator
def gen_sqexp_k(theta):
    A = theta[0]
    d = sp.matrix(theta[1:]).T
    N = sp.matrix(d).shape[0]
    D = sp.eye(N)
    for i in range(N):
        D[i, i] = 1./(d[i, 0]**2)
    return lambda x, y: A*sp.exp(-0.5*(x-y)*D*(x-y).T)


# sqexp kernel generator
def sqexp_k_d(theta, x1, x2, d1=[sp.NaN], d2=[sp.NaN]):
    x1 = sp.matrix(x1)
    x2 = sp.matrix(x2)
    # print "xxxxx"
    A = theta[0]
    d = sp.matrix(theta[1:]).T

    N = sp.matrix(d).shape[0]
    D = sp.eye(N)
    for i in range(N):
        D[i, i] = 1./(d[i, 0]**2)
    X = x1-x2
    # print X
    # print D
    core = (A*sp.exp(-0.5*(X)*D*(X).T))[0, 0]
    # print core
    # d0 is the all the derivitive directions unsorted, then sorted
    d0 = []
    n1 = 0
    n2 = 0
    for i in d2:
        if not sp.isnan(i):
            # print 'x'+str(d2)
            n2 += 1
    # print 'y'+str(n2)
    sign = (-1)**(n2 % 2)
    # print sign
    for i in d1+d2:

        if not sp.isnan(i):
            d0.append(i)
  
    if len(d0) == 0:
        # print "basic kernel"
        # print core
        return core
    if len(d0)==1:
        #print "1st derivative"
        l=D[d0[0],d0[0]]
        
        #print dx
        x=X[0,d0[0]]
        coef = -l*x
        #print core
        return sign*coef*core
    #print d0
    d0.sort()
    #print d0
    im=d0[0]
    #d1 is a list of lists witht each sublist for  derivative direction
    d1=[[im]]
    j=0
    for i in d0[1:]:
        if i==im:
            d1[j].append(i)
            im=i
        else:
            d1.append([i])
            j+=1
            im=i
    #print d1
    #d2 is the number of repeated derivatives in each direction
    d2=[]
    for i in d1:
        d2.append(len(i))
    #print d2
    #d3 is the delta x for each direction
    d3=[]
    #d4 is the lengthscale in each direction
    d4=[]
    for i in d1:
        d3.append((x2[0,i[0]]-x1[0,i[0]]))
        d4.append(D[i[0],i[0]])
        
    

    
    
    d3=[x for (y,x) in sorted(zip(d2,d3))]
    d4=[x for (y,x) in sorted(zip(d2,d4))]
    d2.sort()

   
    #print d2
    #print d3
    #print d4
    
    
    if d2==[1,1]:
        #print "2nd derivative i!=j"
        ai=d1[0][0]
        aj=d1[1][0]
        xi=X[0,ai]
        xj=X[0,aj]
        li=D[ai,ai]
        lj=D[aj,aj]
        
        coef=xi*li*xj*lj
        return sign*coef*core
        
    elif d2==[2]:
        #print "2nd derivative i==j"
        ai=d1[0][0]
        xi=X[0,ai]
        li=D[ai,ai]
        
        coef = li*(li*xi**2-1)
        return sign*coef*core
        
    elif d2==[1,1,1]:
        #print "3nd derivative i!=j!=k"
        ai=d1[0][0]
        aj=d1[1][0]
        ak=d1[2][0]
        xi=X[0,ai]
        xj=X[0,aj]
        xk=X[0,ak]
        li=D[ai,ai]
        lj=D[aj,aj]
        lk=D[ak,ak]
        
        coef=-xi*li*xj*lj*xk*lk
        return sign*coef*core
        
    elif d2==[1,2]:
        #print "3nd derivative i!=j==k"
        ai=d1[0][0]
        aj=d1[1][0]
        
        xi=X[0,ai]
        xj=X[0,aj]
        
        li=D[ai,ai]
        lj=D[aj,aj]
        
        coef = -xi*li*lj*(lj*xj**2-1)
        print sign*coef*core
        
    elif d2==[3]:
        #print "3nd derivative i==j==k"
        ai=d1[0][0]
        xi=X[0,ai]
        li=D[ai,ai]
        
        coef=(li**2)*(3*xi-li*xi**3)
        return sign*coef*core
    
    elif d2==[1,1,1,1]:
        #print "4th derivative i!=j!=k!=l"
        ai=d1[0][0]
        aj=d1[1][0]
        ak=d1[2][0]
        al=d1[3][0]
        xi=X[0,ai]
        xj=X[0,aj]
        xk=X[0,ak]
        xl=X[0,al]
        li=D[ai,ai]
        lj=D[aj,aj]
        lk=D[ak,ak]
        ll=D[al,al]
        
        coef = xi*li*xj*lj*xk*lk*xl*ll
        return sign*coef*core
        
    elif d2==[1,1,2]:
        #print "4th derivative i!=j!=k==l"
        ai=d1[0][0]
        aj=d1[1][0]
        ak=d1[2][0]
        xi=X[0,ai]
        xj=X[0,aj]
        xk=X[0,ak]
        li=D[ai,ai]
        lj=D[aj,aj]
        lk=D[ak,ak]
        
        coef=lk*(lk*xk**2-1)*xi*li*xj*lj
        return sign*coef*core
    
    elif d2==[1,3]:
        #print "4th derivative i!=j==k==l"
        ai=d1[0][0]
        aj=d1[1][0]
        xi=X[0,ai]
        xj=X[0,aj]
        li=D[ai,ai]
        lj=D[aj,aj]
        
        coef=-li*xi*(lj**2)*(3*xj-lj*xj**3)
        return sign*coef*core
        
    elif d2==[2,2]:
        #print "4th derivative i==j!=k==l"
        ai=d1[0][0]
        aj=d1[1][0]
        xi=X[0,ai]
        xj=X[0,aj]
        li=D[ai,ai]
        lj=D[aj,aj]
        
        coef=li*(li*xi**2-1)*lj*(lj*xj**2-1)
        return sign*coef*core
    
    elif d2==[4]:
        #print "4th derivative i==j==k==l"
        ai=d1[0][0]
        xi=X[0,ai]
        li=D[ai,ai]
        
        coef = (li*xi)**4 - 6*(li**3)*(xi**2) + 3*li**2
        return sign*coef*core
    
    else:
        print d2
        raise ValueError("derivative combination not suported")
        
from functools import partial

def gen_sqexp_k_d(theta):
    k = partial(sqexp_k_d,theta)
    return k


#k builder for asymettric k
def buildKasym_d(kf,x1,x2,d1,d2):
        #x should be column vectors, returns a matrix with x1 rows, x2 columns
        (l1,_)=x1.shape
        (l2,_)=x2.shape
        K=sp.matrix(sp.empty([l1,l2]))
        for i in range(l1):
            for j in range(l2):
                K[i,j]=kf(x1[i,:],x2[j,:],d1=d1[i],d2=d2[j])
                
        return K
# builder for symmetric k. 10**-10 is added tot eh diagonal for numerical stabulity

# builder for symmetric k. 10**-10 is added tot eh diagonal for numerical stabulity
def buildKsym_d(kf,x,d):
        #x should be  column vector
        (l,_)=x.shape
        K=sp.matrix(sp.empty([l,l]))
        for i in range(l):
            K[i,i]=kf(x[i,:],x[i,:],d1=d[i],d2=d[i])+10**-10
            for j in range(i+1,l):
                K[i,j]=kf(x[i,:],x[j,:],d1=d[i],d2=d[j])
                
                K[j,i]=K[i,j]
                
        return K
    
class GPcore:
    def __init__(self,X_s,Y_s,S_s,D_s,kf):
        #print "GPcore init"
        #print X_s.T
        #print Y_s.T
        #print S_s.T
        self.D_s=D_s
        self.X_s=X_s
        self.Y_s=Y_s
        #self.S_s=S_s
        self.kf=kf
        
        K_ss = buildKsym_d(kf,X_s,D_s)
        (sign,self.logdet)=slogdet(K_ss)
        self.K_ss_cf = spl.cho_factor(K_ss+vec2trace(S_s))
        self.a = spl.cho_solve(self.K_ss_cf,Y_s)
        
        
        #print "----------"
        return
    
    def changeY(self,Y_s):
        #print "changeY"
        #print Y_s.T
        
        self.a = spl.cho_solve(self.K_ss_cf,Y_s)
        #print "a"
        #print self.a.T
        return
    
    def infer_m(self,X_i,D_i):
        K_si = buildKasym_d(self.kf,X_i,self.X_s,D_i,self.D_s)
        m = K_si*self.a
        return m
    
    def infer_full(self,X_i,D_i):
        K_ii = buildKsym_d(self.kf,X_i,D_i)
        K_si = buildKasym_d(self.kf,X_i,self.X_s,D_i,self.D_s)
        
        m = K_si*self.a
	V = K_ii - K_si*spl.cho_solve(self.K_ss_cf,K_si.T)
        
        return [m,V]
    
    def infer_diag(self,X_i,D_i):
        K_si = buildKasym_d(self.kf,X_i,self.X_s,D_i,self.D_s)
        
        m = K_si*self.a
        V=sp.matrix(sp.zeros(X_i.shape[0])).T
        for j,r in enumerate(K_si):
            
            V[j,0]=self.kf(X_i[j,:],X_i[j,:],d1=D_i[j],d2=D_i[j])-r*spl.cho_solve(self.K_ss_cf,r.T)
        return [m,V]

    def llk(self):
       return (-0.5*self.Y_s.T*self.a -0.5*self.logdet-0.5*len(self.D_s)*sp.log(2*sp.pi))[0,0]
#kf = gen_sqexp_k_d([1.,0.3])
