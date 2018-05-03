
from matplotlib.widgets import Cursor
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from sys import exit

plt.close('all')

grid_size = 10

fig, ax = plt.subplots(1,1)
ax.grid(which='both')
ax.set_aspect('equal')
ax.xaxis.set_major_locator(MaxNLocator(integer=True))
ax.yaxis.set_major_locator(MaxNLocator(integer=True))
ax.axis([ -1, grid_size+1, -grid_size//2-1, grid_size//2+1 ])
xticks = np.arange(-1, grid_size+1, 1)
yticks = np.arange(-grid_size//2-1, grid_size//2+1, 1)
ax.set_xticks(xticks)
ax.set_yticks(yticks)
ax.set_yticklabels([])
ax.set_xticklabels([])  

class AminoAcidChain:
    def __init__(self, fig, ax):
        self.x = []
        self.y = []

        self.hx =[]
        self.hy = []

        self.fig = fig
        self.ax = ax
        
        self.cid = fig.canvas.mpl_connect('button_press_event', self)
        self.bid = fig.canvas.mpl_connect('key_press_event', self.__press__)
        
        self.key_count = 0

        
    def __press__(self, event):
        self.key_count += 1
        if self.key_count == 1:
            print('Now, please define which amino acids are hydrophobic:')
        elif self.key_count >= 2:
            print('The amino acid chain has been recorded:')


    def __call__(self, event):

        def add_point(self, xclick, yclick):
            self.x.append(xclick)
            self.y.append(yclick)
            self.ax.scatter(xclick, yclick, color='darkorange')            
            print("New amino acid added at location:",[int(xclick),int(yclick)])
        
        def add_h(self, xclick, yclick):
            self.hx.append(xclick)
            self.hy.append(yclick)
            self.ax.scatter(xclick, yclick, color='darkcyan')            
            print("Amino acid added at location:",[int(xclick),int(yclick)],"noted as hydrophobic")
        
        def in_list(self, x_list, y_list, xclick, yclick):
            check_x = np.where(x_list == xclick)
            check_y = np.where(y_list == yclick)
            check = len(np.intersect1d(check_x, check_y))
            return not(check == 0)

        xclick = np.round(event.xdata)
        yclick = np.round(event.ydata)

        ## This is the time when the chain is made, but H/P is not declared
        if self.key_count == 0:
            if not self.x:

                add_point(self, xclick, yclick)
                self.fig.canvas.draw()

            else:
                
                r = (xclick - self.x[-1])**2 + (yclick - self.y[-1])**2
                if r > 1.0:
                    return
                else:
                    inchain = in_list(self, self.x, self.y, xclick, yclick)
                    
                    if not(inchain):
                        add_point(self, xclick, yclick)
                        self.ax.plot(self.x[-2:], self.y[-2:],color='black')
                        self.fig.canvas.draw()

                    else:
                        return

        elif self.key_count == 1:
            inchain = in_list(self, self.x, self.y, xclick, yclick)
            if inchain:
                inhydro = in_list(self, self.hx, self.hy, xclick, yclick)
                if not(inhydro):
                    add_h(self, xclick, yclick)
                    self.fig.canvas.draw()

        else:
            return

chain = AminoAcidChain(fig, ax)
plt.show()
