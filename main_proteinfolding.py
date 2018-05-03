import numpy as np
import matplotlib.pyplot as plt
from random import randint, random, sample
from matplotlib.ticker import MaxNLocator

plt.close('all')

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
    plt.title(("Energy: %i"%aa.prev_energy))
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

## define our class which defines the amino acid chain
class AminoAcidChain():
    
    def __init__(self, numAA, H, x = None, y = None):
        self.length = numAA
        if x == None:
            x = np.arange(0,numAA)
        if y == None:
            y = np.zeros(numAA, dtype ='int')

        self.x = x
        self.y = y
        self.temp_x = np.arange(0,numAA)
        self.temp_y = np.zeros(numAA, dtype ='int')
        self.H = np.array(H)
        self.measure_energy(self.x, self.y)

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
        if init_loc >= aa.length//2:
            swap_locs = np.arange(init_loc, aa.length, 1)
        elif init_loc < aa.length//2:
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
                        energy += 0.5
                    
        return energy

## initialize starting linear strand
n = 30
aa = AminoAcidChain( n , sample( range(0, n), n//3 ) )

im_step = 10
save_path = 'temp_images/'

aa.init_energy = aa.measure_energy(aa.x, aa.y)
aa.prev_energy = aa.init_energy

plot_aa(aa)
save_im(aa, 0, save_path)

num_steps = 2000
energy_change = np.zeros(num_steps)
step = 0
prob_equal = 0.0
T = 0.0

while step < num_steps:
    success = aa.move() 

    if success == True:
        aa.temp_energy = aa.measure_energy(aa.temp_x, aa.temp_y)
        
        if (aa.temp_energy - aa.prev_energy < 0) or ((aa.temp_energy == aa.prev_energy) and (random() < prob_equal)) :            
            if random() < T:
                pass
            else:
                print("Step:", step, "Update... Old energy: ", aa.prev_energy, "New energy:", aa.temp_energy)
                
                aa.x = np.array(aa.temp_x)
                aa.y = np.array(aa.temp_y)
                aa.prev_energy = aa.temp_energy
                
                ## Save certain configurations for visualization
                if im_step != None:
            #if step%im_step == 0:
                    save_im(aa, step+1, save_path)                
                
        else:
            pass
    
        energy_change[step] = aa.prev_energy
        step += 1
    
plot_aa(aa)
print("Started at:", aa.init_energy, ", Ended at:", aa.prev_energy)

fig, ax = plt.subplots(1,1)
ax.plot(energy_change)
ax.set_xlabel("Steps")
ax.set_ylabel("Number of Unfavourable Bonds")