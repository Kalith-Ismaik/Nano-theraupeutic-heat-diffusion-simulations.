# python pakage to process and save molecular dynamics simulation files.
import numpy as np
import pandas as pd
from collections import Counter

def read_pdb(filename):
    inputfile = open(filename, "r", encoding='utf-8')
    box    = []
    x_coords = []
    y_coords = []
    z_coords = []
    atomnm = []
    resnm  = []
    resnr  = []
    ele    = []
    ele_num = 0
    
    try:
        for line in inputfile:
            ele_num += 1
            if (line.find("ATOM") == 0 or
                line.find("HETATM") == 0):
                x = float(line[31:37])
                y = float(line[38:45])
                z = float(line[46:53])
                x_coords.append(x)
                y_coords.append(y)
                z_coords.append(z)
                atomnm.append(line[12:16])
                resnm.append(line[17:21])
                resnr.append(ele_num)
                elem = atomnm
                   
            elif (line.find("CRYST1") == 0):
                box.append(float(line[7:15]))
                box.append(float(line[16:24]))
                box.append(float(line[25:33]))
                
        items = list(Counter(resnm).keys())
        l = len(items)+1
        
        for x in items: 
            for y in resnm:
                if x == y:
                    ele.append(items.index(x))
        
    finally:
        inputfile.close()
        
    df = pd.DataFrame(list(zip(x_coords,y_coords,z_coords,ele)),
               columns =['x_coords','y_coords','z_coords','ele'])
    x = df.to_numpy()
    np.savetxt('Initial_conf.out', x)
    return df, box, atomnm, resnm, resnr, elem

def write_out(dataframe_coords,name):
    xx = df.to_numpy()
    np.savetxt(name, x)    


        
def boundary(box, resnr, coords):
    N = len(coords/4)
    for x in range(N):
        for y in range(3):
            coords.iloc[x,y] = (box[y]/2) + coords.iloc[x,y]
    x = coords.to_numpy()
    np.savetxt('positive_conf.out',x)
    return coords
        
#Author      : Kalith M Ismail.
#Objective   : Molecular Dynamic simulation of heat distribution from nanomaterial to surrounding environment by photothermal excitation. 
#Organization: NRNU MEPhI___PhysBIo___Moscow__Russian Federation.
#Date        : 27/11/2021.
#Mentor      : Prof.Dr.Ivanov Dmitry. [University Of Kassel]
