import numpy as np
import matplotlib.pyplot as plt


title = "90 Periodic"
file = "./results/periodic/90_PER.dat"
nnx = 11
nny = 11

def plot():
        extendedC = readFromFile()
        c = np.zeros(( nnx , nny ))
        c = np.array(extendedC[extendedC!=0].tolist()).reshape(-1,nnx)
        
        fig = plt.figure(figsize=(100,100))
        plt.imshow(c, cmap='bwr')
        plt.title(title)
        plt.colorbar()
        plt.show() 

def readFromFile():
    raw = np.genfromtxt(file)
    xrange = ( raw[:,0].min() , raw[:,0].max() )
    yrange = ( raw[:,1].min() , raw[:,1].max() )
    xidx = (raw[:,0]-xrange[0]).astype(int)
    yidx = (raw[:,1]-yrange[0]).astype(int)
    extendedC = np.zeros(( len(xidx) , len(yidx) ))
    extendedC[xidx, yidx] = raw[:,2]
    return extendedC
        
plot()
        