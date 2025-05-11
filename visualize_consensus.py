import sys

import math as m
import numpy as np
import matplotlib.pyplot as plt

T= 160

def get_traj_array(filename, s_or_i = 'state'):
    """
    Parses output of consensus.cc script.

    Toggle s_or_i to 'state' or 'input' to get inputs or states
    """
    traj = []
    inputs = []
    with open(filename) as f:
        l = 0 
        for line in f:
            x = [i for i in line.rstrip(' \n').split(' ')]
            if 'State:' in x:
                x = [float(i) for i in x[1:]]
                traj.append(x)
            if 'Input:' in x:
                x = [float(i) for i in x[1:]]
                inputs.append(x)
            # cur_pre_iter = np.rint(x[0])
    traj = np.array(traj)
    inputs = np.array(inputs)
    if s_or_i == 'state':
        return traj
    elif s_or_i == 'input':
        return inputs
    
t= np.arange(T)

# Extract the trajectories and input
traj = get_traj_array("traj_active1.txt", 'state')[0:T, :]
div_traj = get_traj_array("traj_passive1.txt", 'state')[0:T, :]
inputs = get_traj_array("traj_active1.txt", 'input')[0:T, :]


# Visualize state trajectories 
statefig = plt.figure()
ax = statefig.add_subplot(111)
lw = 0.5
ax.step(t, traj[:,0],
        t, traj[:,1], 
        t, traj[:,2], 
        t, traj[:,3], 
        t, traj[:,4],
        t, traj[:,5],
        t, traj[:,6],
        t, traj[:,7],
        t, traj[:,8],
        t, traj[:,9], linewidth = 2*lw)
ax.step(t, div_traj[:,0], '--',
        t, div_traj[:,1], '--', 
        t, div_traj[:,2], '--', 
        t, div_traj[:,3], '--', 
        t, div_traj[:,4], '--', 
        t, div_traj[:,5], '--',
        t, div_traj[:,6], '--',
        t, div_traj[:,7], '--',
        t, div_traj[:,8], '--',
        t, div_traj[:,9], '--', linewidth = 2*lw, alpha = .4)
ax.plot(t, [10 for _ in t], ':', linewidth = 2*lw, alpha = .9)
ax.plot(t, [50 for _ in t], ':', linewidth = 2*lw, alpha = .9)
ax.set_xlim(0,T)
ax.set_title("State Trajectories")

# Visualize inputs1
inputfig = plt.figure()
ax2 = inputfig.add_subplot(111)
td =.03 # temporal offset for visualization
d = .04 # vertical offset for visualization
ax2.plot(t-4.5*td, inputs[:,0]-4.5*d,
         t-3.5*td, inputs[:,1]-3.5*d, 
         t-2.5*td, inputs[:,2]-2.5*d,
         t-1.5*td, inputs[:,3]-1.5*d, 
         t-.5*td, inputs[:,4]-.5*d,
         t+.5*td, inputs[:,5]+.5*d,
         t+1.5*td, inputs[:,6]+1.5*d, 
         t+2.5*td, inputs[:,7]+2.5*d,
         t+3.5*td, inputs[:,8]+3.5*d, 
         t+4.5*td, inputs[:,9]+4.5*d, linestyle = '-', linewidth = 2*lw, drawstyle='steps-post', marker='o') 
ax2.set_xlim(-td-.1,T-1)
ax2.set_xlim(0,T)
ax2.set_ylim(0,10)
ax2.set_title("Input Trajectories")

plt.show()
