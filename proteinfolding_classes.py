import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from random import random, randint, sample

class ChainBuilder:
    def __init__(self, grid_size):

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


        self.x = []
        self.y = []

        self.hx =[]
        self.hy = []
        self.h = []

        self.fig = fig
        self.ax = ax
        
        self.cid = fig.canvas.mpl_connect('button_press_event', self)
        self.bid = fig.canvas.mpl_connect('key_press_event', self.__press__)
        
        self.key_count = 0

        
    def __press__(self, event):
        self.key_count += 1
        if self.key_count == 1:
            print('Now, please define which amino acids are hydrophobic.')
        elif self.key_count == 2:
            print('The amino acid chain has been recorded.')
        else:
            return

    def __call__(self, event):

        def add_point(self, xclick, yclick):
            self.x.append(xclick)
            self.y.append(yclick)
            self.ax.scatter(xclick, yclick, color='darkorange')            
            # print("New amino acid added at location:",[int(xclick),int(yclick)])
        
        def add_h(self, xclick, yclick):
            self.hx.append(xclick)
            self.hy.append(yclick)
            self.ax.scatter(xclick, yclick, color='darkcyan')            
            # print("Amino acid added at location:",[int(xclick),int(yclick)],"noted as hydrophobic")
        
        def in_list(self, x_list, y_list, xclick, yclick):
            check_x = np.where(x_list == xclick)
            check_y = np.where(y_list == yclick)
            check = np.intersect1d(check_x, check_y)
            return not(len(check) == 0), check

        xclick = np.rint(event.xdata)
        yclick = np.rint(event.ydata)

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
                    inchain = in_list(self, self.x, self.y, xclick, yclick)[0]
                    
                    if not(inchain):
                        add_point(self, xclick, yclick)
                        self.ax.plot(self.x[-2:], self.y[-2:],color='black')
                        self.fig.canvas.draw()

                    else:
                        return

        elif self.key_count == 1:
            inchain, ind = in_list(self, self.x, self.y, xclick, yclick)
            if inchain:
                inhydro = in_list(self, self.hx, self.hy, xclick, yclick)[0]
                if not(inhydro):
                    self.h.append(ind[0])
                    add_h(self, xclick, yclick)
                    self.fig.canvas.draw()

        else:
            return



#####################################################



## define our class which defines the amino acid chain
class AminoAcidChain():
    
    def __init__(self, numAA, H, x = None, y = None):
        self.length = numAA
        if np.all(x) == None:
            x = np.arange(0,numAA)
        if np.all(y) == None:
            y = np.zeros(numAA, dtype ='int')

        self.x = x
        self.y = y
        self.temp_x = np.arange(0,numAA)
        self.temp_y = np.zeros(numAA, dtype ='int')
        self.H = np.array(H)
        self.prev_energy = self.measure_energy(self.x, self.y)

    def new_loc(self, x, y, rot):
        if rot == 0:
            new_x = x
            new_y = y + 1
        elif rot == 1:
            new_x = x + 1 
            new_y = y
        elif rot == 2:
            new_x = x
            new_y = y - 1
        elif rot == 3:
            new_x = x - 1 
            new_y = y
        return new_x, new_y
    
    def check_occupied(self, x, y, new_x, new_y):
        check_x = np.where(x == new_x)
        check_y = np.where(y == new_y)
        
        adj_aa = np.intersect1d(check_x, check_y)
        
        if len(adj_aa) == 0:
            return False, adj_aa ## is not occupied
        else:
            return True, adj_aa  ## is occupied
        
    def move(self):
        
        ## a temporary list of amino acid locations
        temp_x = np.array(self.x)
        temp_y = np.array(self.y)
        
        ## choose a random aa to start the fold at
        init_loc = randint(0, self.length-1)
                
        ## the amino acids which will move (from random point to end, shortest direction)
        if init_loc >= self.length//2:
            swap_locs = np.arange(init_loc, self.length, 1)
        elif init_loc < self.length//2:
            swap_locs = np.arange(init_loc, -1, -1)                
                
        for i_loc in range(len(swap_locs)-1):
            rot = randint(0, 3)
            
            curr_loc = swap_locs[i_loc]
            next_loc = swap_locs[i_loc + 1]
            
            ## try all four rotations that are possible, until one is found
            for j_rot in [0,1,2,3]:
                temp_rot = (rot + j_rot) % 4

                new_x, new_y = self.new_loc(temp_x[curr_loc], temp_y[curr_loc], temp_rot )
                occ, adj_aa = self.check_occupied(temp_x, temp_y, new_x, new_y)
                                             
                ## if the new location is not occupied, break loop
                if occ == False:
                    temp_x[next_loc] = new_x
                    temp_y[next_loc] = new_y
                    break
                
                ## failed, cannot move anywhere
                elif occ == True and j_rot == 3:
                    success = False
                    return success

        success = True
        self.temp_x = temp_x
        self.temp_y = temp_y
        return success
        
    def measure_energy(self, x, y):
        energy = 0
        for Hloc in self.H:
            ## check 
            for rot in [0,1,2,3]:
                new_x, new_y = self.new_loc(x[Hloc], y[Hloc], rot)
                occ, adj_aa = self.check_occupied(x, y, new_x, new_y)
                
                ## if there is no amino acid in that bond location, add energy quanta (is next to water)
                if occ == False:
                    energy += 1
                
                ## if there is an amino acid in adjacent location
                elif occ == True:
                    ## is it another hydrophobic AA?
                    if np.isin(adj_aa, self.H):
                        pass
                        # do nothing, as this is favourable
                        
                    elif adj_aa == Hloc + 1 or adj_aa == Hloc - 1:
                        pass
                        # do nothing, also favourable
                        
                    else:
                        energy += 1.0 # this value could change if desired
                    
        return energy