# -*- coding: utf-8 -*-
"""
Created on Wed May  6 14:51:54 2015

class to create a set of GPd andor GPep and make parallel inference over them
@author: mark
"""

from multiprocessing import Process, Pipe, Event
import GPd
import GPep
import time
import traceback
import scipy as sp
import scipy.linalg as spl
from scipy.stats import norm as norm

class multiGP:
    def __init__(self):
        self.processes = []
        self.conns = []
        self.exits = []
        return
    
    def addGPd(self, X_s, Y_s, S_s, D_s, kf):
        
        pconn,cconn=Pipe()
        ex = Event()
        p = mGPd(cconn,ex, X_s, Y_s, S_s, D_s, kf)
    
        self.processes.append(p)
        self.conns.append(pconn)
        self.exits.append(ex)
        
        p.start()
        return
        
    def addGPep(self, X_c, Y_c, S_c, D_c, X_z, D_z, I_z, G_z, N_z, kf):
        
        pconn,cconn=Pipe()
        ex = Event()
        p = mGPep(cconn,ex, X_c, Y_c, S_c, D_c, X_z, D_z, I_z, G_z, N_z, kf)
    
        self.processes.append(p)
        self.conns.append(pconn)
        self.exits.append(ex)
        
        p.start()
        return
        
    def infer_m(self, X_i, D_i):
        return self.infer_any(2, X_i, D_i)
        
    def infer_diag(self, X_i, D_i):
        return self.infer_any(1, X_i, D_i)
        
        
    def infer_full(self, X_i, D_i):
        return self.infer_any(0, X_i, D_i)
    
    def infer_any(self,code, X_i, D_i):
        for c in self.conns:
            c.send([code,X_i,D_i])

        result=[]
        for c in self.conns:
            while not c.poll():
                time.sleep(0.1)
            result.append(c.recv())
        return result
      
    def llk(self):
        return self.infer_any(3, [], [])
    
    def status(self):
        return self.infer_any(-1,[],[])
    
    def drawmins(self,n_points,bound):
        return self.infer_any(4,n_points,bound)
        
    def close(self):
        for ex in self.exits:
            ex.set()
        for c in self.conns:
            c.close()
        for proc in self.processes:
            proc.join()

        return


class mGPd(Process):
    def __init__(self,conn,exit_event, X_s, Y_s, S_s, D_s, kf):
        super(mGPd, self).__init__()
        self.conn=conn
        self.exit_event=exit_event
        self.GP = GPd.GPcore(X_s, Y_s, S_s, D_s, kf,precom=False)
        self.status=[0,0,0]
        return
        
    def run(self):
        try:
            self.GP.precompute()
        except:
            err = traceback.format_exc()
            self.status=[-2,err]
            while not self.exit_event.is_set():
                if not self.conn.poll():
                    time.sleep(0.1)
                    continue
                msg = self.conn.poll
                self.conn.send(self.status)
                
        while not self.exit_event.is_set():
            if not self.conn.poll():
                time.sleep(0.1)
                continue
            [code, X_i, D_i] = self.conn.recv()
            try:
                if code==-1:
                    res=self.status
                elif code==0:
                    res = self.GP.infer_full(X_i,D_i)
                elif code==1:
                    res = self.GP.infer_diag(X_i,D_i)
                elif code==2:
                    res = self.GP.infer_m(X_i,D_i)
                elif code==3:
                    res = self.GP.llk()
                elif code==4:
                    n_points=X_i
                    bound=D_i
                    res=self.drawmin(n_points,bound)
                else:
                    raise ValueError('code not supported')
                self.conn.send([0, res])
                self.status[1]+=1
            except:
                self.conn.send([-1,traceback.format_exc()])
                self.status[0]=-1
                self.status[1]+=1
                self.status[2]+=1
        return
        
    def drawmin(self,n_points,bound):
        ub=bound[1]
        lb=bound[0]
        mx = self.GP.Y_s.max()
        mn = self.GP.Y_s.min()
        nacc=0
        X_tmp=[]
        sp.random.seed()
        while nacc<n_points:
            X_prop = sp.matrix(sp.random.uniform(lb, ub)).T
            Y_prop = self.GP.infer_m(X_prop, [[sp.NaN]])
            theta = -(Y_prop[0,0]-mx)/(mx-mn)
            p = norm.cdf(2*theta-1.)
            
            if sp.random.uniform(0,1)<=p:
                nacc+=1
                X_tmp.append(X_prop)
        X_x  = sp.vstack(X_tmp)
        D_x = [[sp.NaN]] * n_points
        
        mh, Vh = self.GP.infer_full(X_x, D_x)
        Vh_cho = spl.cholesky(Vh, lower=True)
        dr = mh+Vh_cho*sp.matrix(sp.random.normal(size=n_points)).T
        xs = dr.argmin()
        return X_x[xs,:]
        
class mGPep(Process):
    def __init__(self,conn,exit_event,X_c, Y_c, S_c, D_c, X_z, D_z, I_z, G_z, N_z, kf):
        super(mGPep, self).__init__()
        self.conn=conn
        self.exit_event=exit_event
        self.GP = GPep.GPcore(X_c, Y_c, S_c, D_c, X_z, D_z, I_z, G_z, N_z, kf)
        self.status=[0,0,0]
        return
        
    def run(self):
        try:
            self.GP.runEP()
        except:
            err = traceback.format_exc()
            self.status=[-2,err]
            while not self.exit_event.is_set():
                if not self.conn.poll():
                    time.sleep(0.1)
                    continue
                msg = self.conn.poll
                
                self.conn.send(self.status)
                
        while not self.exit_event.is_set():
            if not self.conn.poll():
                time.sleep(0.1)
                continue
            [code, X_i, D_i] = self.conn.recv()
            try:
                if code==-1:
                    res=self.status
                elif code==0:
                    res = self.GP.infer_full(X_i,D_i)
                elif code==1:
                    res = self.GP.infer_diag(X_i,D_i)
                elif code==2:
                    res = self.GP.infer_m(X_i,D_i)
                elif code==3:
                    res = self.GP.llk()
                else:
                    raise ValueError('code not supported')
                self.conn.send([0, res])
                self.status[1]+=1
            except:
                self.conn.send([-1,traceback.format_exc()])
                self.status[0]=-1
                self.status[1]+=1
                self.status[2]+=1
        return