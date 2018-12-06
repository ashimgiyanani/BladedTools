# -*- coding: utf-8 -*-
"""
def  readBladed(filename)
%readBladed - reads Bladed format output *.wnd
%
% Syntax:
import readBladed as rb
filename = "d:\Windwise\Bladed\KaimalWind.wnd"
(u, v, w, x, y, z, t, U) = rb.readBladed(filename)
%
% Inputs:
filename - name of the Bladed file with the path
%
% Outputs:
    u,v,w - lonitudinal, lateral and vertical components of wind
    x,y,z - dimensions of the domain
        t - timestamp
        U - mean wind speed at hub height
% Example:

import sys
new_path = 'd:\Windwise\PostProcessing'
if new_path not in sys.path:
    sys.path.append(new_path)
    
import readBladed as rb
   
filename = "d:\Windwise\Bladed\KaimalWind.wnd"
(u, v, w, x, y, z, t, U) = rb.readBladed(filename)


See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: Ashim Giyanani, PhD candidate, TU Delft
% Wind energy department, TU Delft
% email address: ashimgiyanani@yahoo.com
% March 2018; Last revision: 06-04-2018
"""

# Script to read bladed format in python
def readBladed(filename):

#------------- BEGIN CODE --------------

    import numpy as np
    import sys
    
    with open(filename, "r") as f:
        New = np.fromfile(f, dtype=np.int16, count=1)
        Spec = np.fromfile(f,dtype=np.uint16, count=1)
    
        if (Spec==3) | (Spec==5) | (Spec==7):
            Nturb = 3
        else:
            Nturb = 1
    
        if (Spec==4):
          Nturb = np.fromfile(f,dtype=np.uint32, count=1)
          Lattitude = np.fromfile(f,dtype=np.float32, count=1)
          z0 = np.fromfile(f,dtype=np.float32, count=1)
          H = np.fromfile(f,dtype=np.float32, count=1)
          TI_u = np.fromfile(f,dtype=np.float32, count=1)
          TI_v = np.fromfile(f,dtype=np.float32, count=1)
          TI_w = np.fromfile(f,dtype=np.float32, count=1)
    
        if (Spec==7):  
            unknown_1 =np.fromfile(f,dtype=np.int32, count=1)
            unknown_2 =np.fromfile(f,dtype=np.int32, count=1)
    
        delta_z = np.fromfile(f,dtype=np.float32, count=1)
        delta_y = np.fromfile(f,dtype=np.float32, count=1)
        delta_x = np.fromfile(f,dtype=np.float32, count=1)
    
        Nx2 = np.fromfile(f, dtype=np.uint32, count=1)
        Nx = np.int(2*Nx2)
        
        U = np.fromfile(f, dtype=np.float32, count=1)
        Luz = np.fromfile(f, dtype=np.float32, count=1)
        Luy = np.fromfile(f, dtype=np.float32, count=1)
        Lux = np.fromfile(f, dtype=np.float32, count=1)
        
        Dummy = np.fromfile(f, dtype=np.float32, count=1)
        Seed = np.fromfile(f, dtype=np.float32, count=1)
    
        Nz = np.int(np.fromfile(f, dtype=np.uint32, count=1))
        Ny = np.int(np.fromfile(f, dtype=np.uint32, count=1))
    
        # dimensions
        x = (np.arange(0,(Nx),1))*delta_x
        # BLADED: first element at right bottom corner, i.e. maximum y
        y=-(np.arange(-(Ny-1)/2,(Ny-1)/2+1))*np.int(delta_y)
        z= (np.arange(-(Nz-1)/2,(Nz-1)/2+1))*np.int(delta_z)
        t=x/U
    
        if (Nturb==3) | (Spec==7):
          Lvz=np.fromfile(f,dtype=np.float32, count=1)
          Lvy=np.fromfile(f,dtype=np.float32, count=1)
          Lvx=np.fromfile(f,dtype=np.float32, count=1)
          Lwz=np.fromfile(f,dtype=np.float32, count=1)
          Lwy=np.fromfile(f,dtype=np.float32, count=1)
          Lwx=np.fromfile(f,dtype=np.float32, count=1)
        
        
        if (Spec==7):
            Alpha = np.fromfile(f,dtype=np.float32, count=1)
            Gamma = np.fromfile(f,dtype=np.float32, count=1)
            
        u=np.zeros(shape=(np.int(Nz),np.int(Ny),np.int(Nx)))
        v=np.zeros(shape=(np.int(Nz),np.int(Ny),np.int(Nx)))
        w=np.zeros(shape=(np.int(Nz),np.int(Ny),np.int(Nx)))
        Nplane=Ny*Nz*Nturb
        Ind1 = (np.arange(0,Nplane,3))
        Ind2 = (np.arange(1,Nplane,3))
        Ind3 = (np.arange(2,Nplane,3))
    
        for i in np.arange(0,Nx,1):
           uvw=np.fromfile(f,dtype=np.int16, count=np.int(Nplane))
           if (Nturb==3):
             uu = uvw[Ind1]
             vv = uvw[Ind2]
             ww = uvw[Ind3]
             v[:,:,i] = vv.reshape(Ny, Nz, order='C')
             w[:,:,i] = ww.reshape(Ny, Nz, order='C')
           elif (Nturb==1):
             uu=uvw
           else:
               sys.exit("wrong number of Nturb")
           u[:,:,i] = uu.reshape(Ny, Nz, order='C')  
    
        u = u/1000;
        v = v/1000;
        w = w/1000;
    return [u, v, w, x, y, z, t, U, delta_x, delta_y, delta_z, Nx, Ny, Nz, Nturb, uvw] ;

#------------- END OF CODE --------------
