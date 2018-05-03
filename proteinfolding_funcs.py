import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def plot_aa(aa):
    fig, ax = plt.subplots(1,1)
    ax.grid(which='both')
    ax.set_aspect('equal')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.axis([ min(aa.x)-1, max(aa.x)+1, min(aa.y)-1, max(aa.y)+1 ])
    #ax.axis([ -aa.length-1, aa.length+1, -aa.length-1, aa.length+1 ])
    xticks = np.arange(min(aa.x)-1, max(aa.x)+1, 1)
    yticks = np.arange(min(aa.y)-1, max(aa.y)+1, 1)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_yticklabels([])
    ax.set_xticklabels([])  
    
    ax.scatter(aa.x, aa.y, label='Polar/Hydrophillic', color='darkorange')
    ax.scatter(aa.x[aa.H], aa.y[aa.H], label='Hydrophobic', color='darkcyan')
    for i in range(0, aa.length-1):
        ax.plot([ aa.x[i], aa.x[i+1] ], [aa.y[i], aa.y[i+1] ], '-k')
    plt.title(("Energy: %i"%aa.measure_energy(aa.x, aa.y)))
    plt.legend()
    plt.draw()

def save_im(aa, num, save_path):
    fig, ax = plt.subplots(1,1)
    ax.grid(which='both')
    ax.set_aspect('equal')
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.axis([ -1, aa.length+1, -aa.length//2-1, aa.length//2+1 ])
    xticks = np.arange(-1, aa.length+1, 1)
    yticks = np.arange(-aa.length//2-1, aa.length//2+1, 1)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_yticklabels([])
    ax.set_xticklabels([])  
    
    ax.scatter(aa.temp_x, aa.temp_y, label='Polar/Hydrophillic')
    ax.scatter(aa.temp_x[aa.H], aa.temp_y[aa.H], label='Hydrophobic')
    for i in range(0, aa.length-1):
        ax.plot([ aa.temp_x[i], aa.temp_x[i+1] ], [aa.temp_y[i], aa.temp_y[i+1] ], '-k')
    plt.title(("Energy: %i"%aa.measure_energy(aa.temp_x, aa.temp_y)))
    
    plt.savefig('%s/image_%i.png'%(save_path,num), bbox_inches='tight')
    plt.close(fig)