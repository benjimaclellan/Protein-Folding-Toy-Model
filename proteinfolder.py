from proteinfolding_funcs import plot_aa, save_im
from proteinfolding_classes import AminoAcidChain, ChainBuilder

from random import random, randint, sample
import matplotlib.pyplot as plt
import numpy as np

plt.close('all')

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