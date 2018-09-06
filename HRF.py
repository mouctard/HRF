import numpy as np
import mdtraj as md
import sys
from core import corr as corr
from core import diffusion as diff

"""
This is an example of how you can start putting code together to make a nice little
software package. There are much better ways of doing this, but this is a good
way to get started. I recommend adding any code to this that requires looping 
through a trajectory. It makes less sense to add anything that doesn't, but in
principle you could. I added the correlation code as an example, but it's not useful
for anything so I don't recommend using it.

DMF 
8/2018

"""

# here we specify a list of calculations we want to run
calcs = ['diffusion'] #diffusion, correlation

# import things per usual
dcd = sys.argv[1] ; psf = sys.argv[2]
if dcd == "bulk_TRG_GOL_last_114ns.dcd": 
	print "not a good test case for diffusion"
	print "use dcd+psf in diffusion folder"
	if 'diffusion' in calcs: exit()

# load traj
traj = md.load_dcd(dcd, top=psf)[450:]

# pre-allocate memory
if 'diffusion' in calcs:
	R_list = np.array([0, 0.25, 0.5, 0.75, 1.0, 1.25,1.5])
	diff_data = np.zeros((len(traj)-1,len(R_list)))

# loop through traj
for i,fr in enumerate(traj):
	print 'FRAME', i+1, 'of', len(traj)
	# ONLY run the things specified in calcs
	if 'diffusion' in calcs:
		# only start on second frame
		if i>0:
			for ri in range(1,len(R_list)):
				# 'diff' is the name of our code, 'run' is the function that runs it
				# since we start recording at i=1, we can shift everything so diff_data[0] isn't just a bunch of zeros
				diff_data[i-1][ri] = diff.D_shells(fr,traj[i-1],R_list[ri-1],R_list[ri])
	if 'correlation' in calcs:
		c = corr.correlation_function(traj)
		corr.plot_corr(c)

# plot
if 'diffusion' in calcs:
	print 'PLOTTING: diffusion'
	D = np.mean(diff_data.T,axis=1)[1:]
	R = [(R_list[i]+R_list[i-1])/2. for i in range(1,len(R_list))]
	diff.plot_shells(R,D)
