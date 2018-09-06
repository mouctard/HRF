import matplotlib as mpl
font = {'family' : 'normal','weight' : 'normal','size'   : 15} ; mpl.rc('font', **font) ; mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import mdtraj as md

def pbc_dist(r_i,r_j,pbc):
        d = r_i - r_j
        for q in range(3):
                if abs(d[q]) > pbc[q]/2.0:
                        if d[q] > 0: d[q] -= pbc[q]
                        else: d[q] += pbc[q]
        return np.array(d)

def plot_corr(data):
	cmap = plt.hist2d(data.T[0], data.T[1],bins=40) ; plt.colorbar()
	plt.subplots_adjust(top=0.9,bottom=0.15,left=0.2,right=0.8,hspace=0.2,wspace=0.2)
	plt.xlabel("Distance (nm)") ; plt.ylabel("cos(theta)")# ; plt.show()
	plt.savefig("corr_bulk.png")
	return None

def correlation_function(traj, sel="resname GOL"):
	dipole = 'n' ; data = []
	if sel == "resname GOL":
		charges = np.array([0.050000,  0.090000, 0.090000, -0.650000, 0.420000, 0.140000,\
				0.090000, -0.650000, 0.420000, 0.050000, 0.090000, 0.090000,\
				-0.650000, 0.420000])
	elif sel == "resname TRG":
		charges = np.array([-0.760000, -0.760000, -0.312000,  0.157000,  0.056000,  0.033000,\
				-0.025000,  0.017000, -0.102000,  0.696000,  0.199000,  0.193000,\
				0.183000,  0.108000,  0.108000,  0.108000,  0.101000])
	else: print sel, "not a known selection"
	pbc = np.array(traj[0].unitcell_lengths)[0]
	atoms = traj[0].topology.select(sel) ; atom_slice = traj[0].atom_slice(atoms)
	res_ind = atoms.reshape(atom_slice.n_residues, atom_slice.n_atoms/atom_slice.n_residues)
	for fr in traj:
		if dipole == 'y':
			mat = np.array([[md.dipole_moments(fr.atom_slice(res),charges)[0],md.compute_center_of_mass(fr.atom_slice(res))[0]] for res in res_ind])
		elif sel == "resname GOL":
			mat = np.array([[pbc_dist(fr.atom_slice([res[0]]).xyz[0][0],fr.atom_slice([res[9]]).xyz[0][0],pbc),md.compute_center_of_mass(fr.atom_slice(res))[0]] for res in res_ind])
		else: print "these options dont make sense"
		for i,mi in enumerate(mat):
			for j,mj in enumerate(mat):
				if i != j:
					r = np.linalg.norm(pbc_dist(mi[1], mj[1], pbc))
					S = np.abs(np.dot(mi[0], mj[0])/(np.linalg.norm(mi[0])*np.linalg.norm(mj[0])))
					data.append([r,S])
	return np.array(data)

dcd = "bulk_TRG_GOL_last_114ns.dcd" ; psf = "2copies_autopsf.psf"
#dcd = "/home/dillion/data/HRF/ABF/simulations/GOL_GOL-TRG/7_TRG/solvate.eq310.0002.coor.dcd" ; psf = "/home/dillion/data/HRF/ABF/simulations/GOL_GOL-TRG/7_TRG/solvate.psf"

# if we are running the code directly, then python assigns a variable called "__name__ = __main__"
# if you import this code into another code, __name__ is something else, and this won't run
if __name__ == "__main__":
	traj = md.load_dcd(dcd,top=psf) 
	data = correlation_function(traj)
	plot_corr(data)

