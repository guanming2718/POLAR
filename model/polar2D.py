#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 13 00:01:48 2020

@author: richard
"""
import configparser
import numpy as np

config = configparser.ConfigParser()
config.read('input.ini')

class Model:
    def __init__(self,_dict):
        config_params = ['nx','ny','Lx','Ly','init-config','BC_rho','BC_polar','BC_velocity']
        model_params = ['Ac','Bc','Kc','Ap','Bp','Kp','w','gamma','nu','zeta',
                       'alpha','Ks','Kv','xi']
        for p in config_params:
            if p in _dict:
                self.__dict__[p]=_dict[p]
            else:
                raise Exception(p +' is not specified' )
        if (_dict['init-config'] =='cluster') or (_dict['init-config'] == 'circular-wound'):
            self.__dict__['R'] = _dict['R']
        elif (_dict['init-config'] == 'two-side-wound') or (_dict['init-config'] == 'one-side-wound'):
            self.__dict__['wound-ratio'] = _dict['wound-ratio']
        for p in model_params:
            if p in _dict:
                self.__dict__[p]=_dict[p]
            else:
                raise Exception(p +' is not specified' )

        #initial condition 2D,only the interior points are calculated
        self.rho = np.ones(self.nx+2,self.ny+2)
        self.px = np.zeros(self.nx+2,self.ny+2)
        self.py = np.zeros(self.nx+2,self.ny+2)
        self.vx = np.zeros(self.nx+2,self.ny+2)
        self.vy = np.zeros(self.nx+2,self.ny+2)
        self.dx,self.dy= self.Lx/(self.nx+1.0),self.Ly/(self.ny+1.0)

    def configure(self):
        if self.config=='confluent':
            a= 1
        elif self.config=='cluster':
            self.rho = np.zeros(self.nx+2,self.ny+2)
            cx = int(self.nx/2)
            cy = int(self.ny/2)
            for i in range(0,int(self.R/self.dx)+1):
                for j in range(0,int(self.R/self.dx)+1):
                    if   i*i*self.dx*self.dx+j*j*self.dy*self.dy <self.R*self.R:
                        self.rho[cx+i,cy+j] = 1.0
                        self.rho[cx+i,cy-j] = 1.0
                        self.rho[cx-i,cy+j] = 1.0
                        self.rho[cx-i,cy-j] = 1.0
        elif self.config=='circular-wound':
            self.rho = np.ones(self.nx+2,self.ny+2)
            cx = int(self.nx/2)
            cy = int(self.ny/2)
            for i in range(0,int(self.R/self.dx)+1):
                for j in range(0,int(self.R/self.dx)+1):
                    if   i*i*self.dx*self.dx+j*j*self.dy*self.dy <self.R*self.R:
                        self.rho[cx+i,cy+j] = 0.0
                        self.rho[cx+i,cy-j] = 0.0
                        self.rho[cx-i,cy+j] = 0.0
                        self.rho[cx-i,cy-j] = 0.0
        elif self.config=='one-side-wound':
            self.rho = np.ones(self.nx+2,self.ny+2)
            r = self.Lx*self.wound_ratio
            for i in range(int(r/self.dx)+1):
                for j in range(self.ny+2):
                        self.rho[i,j] = 1.0
        elif self.config=='two-side-wound':
            self.rho = np.ones(self.nx+2,self.ny+2)
            l = 0.5*self.Lx*(1 - self.wound_ratio)
            r = self.Lx - 0.5*self.Lx*(1 - self.wound_ratio)
            for i in range(int(l/self.dx)-1,int(r/self.dx)+1):
                for j in range(self.ny+2):
                        self.rho[i,j] = 1.0
    def calculate(self):
        '''
        calculate the 1 free energy 2molecular field 3 stress 4 velocity
        '''
        drho_x,drho_y,d2rho_x,d2rho_y = self.grad(self.rho,component = 'scalar',order = 2,BC=self.BC_rho)
        dpx_x,dpx_y,d2px_x,d2px_y = self.grad(self.px,component = 'scalar',order = 2,BC=self.BC_polar)
        dpy_x,dpy_y,d2py_x,d2py_y = self.grad(self.py,component = 'scalar',order = 2,BC=self.BC_polar)
        p2 = self.px*self.px + self.py*self.py
        #free energy of density
        fd = 0.5*self.Ac*np.power(self.rho*(self.rho-1),2.0) + 0.5*self.Kc*(drho_x*drho_x + drho_y*drho_y)
        #free energy of polar field
        fp = 0.5*self.Ap*p2 + 0.25*self.Bp*p2*p2 + 0.5*self.Kp*(dpx_x*dpx_x+dpx_y*dpx_y + dpy_y*dpy_y + dpy_x*dpy_x)
        #coupling
        fdp = -self.w*(self.rho-1)*(dpx_x + dpy_y)
        #total free energy
        f = fd + fp + fdp
        # calculate the molecular field dleta f/ dleta p
        hx = self.Ap*self.px + self.Bp*p2*self.px - self.Kp*(d2px_x + d2px_y)- w*drho_x
        hy = self.Ap*self.py + self.Bp*p2*self.py - self.Kp*(d2py_x + d2py_y)- w*drho_y
        # calulate the chemical potential
        mu = self.Ac*self.rho*(self.rho-1.0)*(2.0*self.rho-1.0) - self.Kc*(d2rho_x+d2rho_y) - w*(dpx_x+dpy_y)
        
        
            
            
    def grad(self,f,order,component,BC):
        # Take derivatives with certain conditions
        # f is the field
        # order = 1 or 2 means first derivative or second derivatives
        # component = 'x' means f is the x component of a vector field
        #             'y' means f is the y component of a vector field
        #             'scalar' means f is a scalar field  
        ff= np.zeros(self.nx+2,self.ny+2)
        # peridoic beoundary condition
        # only compute the points in 0:nx,0:ny
        # f(0,y) = f(nx,y) & f(y=0) = f(y=Ly)
        ff= np.zeros(self.nx+2,self.ny+2)
        if BC=='periodic':
            ff= np.zeros(self.nx+3,self.ny+3)
            ff[1:-1,1:-1] = f[1:-1,1:-1]
            ff[0,1:-1] = f[-1,:]
            ff[-1,1:-1] = f[0,:]
            ff[1:-1,0] = f[:,-1]
            ff[1:-1,-1] = f[:,0]
        # for drho/dn = 0
        elif (BC=='no-flux') and (component=='scalar'):
            ff= np.zeros(self.nx+4,self.ny+4)
            ff[1:-1,1:-1] = f[1:-1,1:-1]
            ff[0,1:-1] = f[1,:]
            ff[-1,1:-1] = f[-2,:]
            ff[1:-1,0] = f[:,1]
            ff[1:-1,-1] = f[:,-2]
        # for rho(x=0 or Lx)= rho(y=0 or Ly)=0 
        elif (BC=='zero'):
            ff= np.zeros(self.nx+4,self.ny+4)
            ff[1:-1,1:-1] = f[1:-1,1:-1]
            ff[0,:] = np.zeros(self.ny+4)
            ff[-1,:] = np.zeros(self.ny+4)
            ff[:,0] = np.zeros(self.nx+4)
            ff[:,-1] = np.zeros(self.nx+4)
        # dv/dn = 0 and v dot n =0 
        elif (BC=='free-slip') and (component=='x'):
            #vx = 0 at x=0 and x=Lx
            ff= np.zeros(self.nx+4,self.ny+4)
            ff[1:-1,1:-1] = f[1:-1,1:-1]
            ff[0,:] = np.zeros(self.ny+4)
            ff[-1,:] = np.zeros(self.ny+4)
            ff[1:-1,0] = f[:,1]
            ff[1:-1,-1] = f[:,-2]
        elif (BC=='free-slip') and (component=='y'):
            #vx = 0 at x=0 and x=Lx
            ff= np.zeros(self.nx+4,self.ny+4)
            ff[1:-1,1:-1] = f[1:-1,1:-1]
            ff[:,0] = np.zeros(self.ny+4)
            ff[:,-1] = np.zeros(self.ny+4)
            ff[0,1:-1] = f[1,:]
            ff[-1,1:-1] = f[-2,:]
        # 1st order centre difference d/dx = 1/dx*(f(i+1)-f(i-1))
        df_dx = np.gradient(ff,axis=0)
        df_dy = np.gradient(ff,axis=1)
        if order ==1: 
            return [ df_dx[1:-1,1:-1], df_dy[1:-1,1:-1] ]
        if order ==2:# d2/dx2 = 1/dx^2 *(f(i+1)-2f(i)+f(i-1))
            df2_dx2 = np.diff(ff,axis=0) - np.diff(np.roll(ff,1,axis=0),axis=0)
            df2_dy2 = np.diff(ff,axis=1) - np.diff(np.roll(ff,1,axis=1),axis=1)
            return [ df_dx[1:-1,1:-1], df_dy[1:-1,1:-1], df2_dx2[1:,:],df2_dy2[:,1:] ]
        
        
        
        
