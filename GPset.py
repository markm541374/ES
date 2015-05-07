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
        p = mGPd(cconn,ex, X_c, Y_c, S_c, D_c, X_z, D_z, I_z, G_z, N_z, kf)
    
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
        return
        
    def run(self):
        try:
            self.GP.precompute()
        except:
            err = traceback.format_exc()
            while not self.exit_event.is_set():
                if not self.conn.poll():
                    time.sleep(0.1)
                    continue
                msg = self.conn.poll
                self.conn.send([-2,err])
                
        while not self.exit_event.is_set():
            if not self.conn.poll():
                time.sleep(0.1)
                continue
            [code, X_i, D_i] = self.conn.recv()
            try:
                
                if code==0:
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
            except:
                self.conn.send([-1,traceback.format_exc()])


        return
        
class mGPep(Process):
    def __init__(self,conn,exit_event,X_c, Y_c, S_c, D_c, X_z, D_z, I_z, G_z, N_z, kf):
        super(mGPd, self).__init__()
        self.conn=conn
        self.exit_event=exit_event
        self.GP = GPep.GPcore(X_c, Y_c, S_c, D_c, X_z, D_z, I_z, G_z, N_z, kf)
        return
        
    def run(self):
        try:
            self.GP.runEP()
        except:
            err = traceback.format_exc()
            while not self.exit_event.is_set():
                if not self.conn.poll():
                    time.sleep(0.1)
                    continue
                msg = self.conn.poll
                self.conn.send([-2,err])
                
        while not self.exit_event.is_set():
            if not self.conn.poll():
                time.sleep(0.1)
                continue
            [code, X_i, D_i] = self.conn.recv()
            try:
                
                if code==0:
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
            except:
                self.conn.send([-1,traceback.format_exc()])


        return