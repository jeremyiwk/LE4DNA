# coding: utf-8

def index_P(TOP):
	idx = []
	with open(TOP, "r+") as f:
		for line in f:
			if line[:4] == 'ATOM':
				if line.split()[2] == 'P':
					idx.append(line.split()[1])

	f.close()

	with open('P.ndx', "w+") as f:
		f.write('[ P ]\n')
		for i in idx:
			f.write(str(i) + '\n')
	f.close()
					

def dump_traj(XTC, TOP, SYS):
	import subprocess
	
	#subprocess.call('module load gromacs', shell = True)
	subprocess.call('gmx_mpi trjconv -f ' + str(XTC) + ' -s ' + str(TOP) + ' -o ' + str(SYS) + '.g96 -n P.ndx', shell = True)

def gen_protinfo(PROTNAME,G96,TOP):
	import subprocess
	import numpy as np

	protname = str(PROTNAME)
	#The way I am using awk requires that TOP be a PDB file
	if TOP.split(".")[1] != 'pdb':
		print('Structure file must be in PDB format for gen_protinfo to work correctly.')
	N = int(subprocess.check_output("awk '{ print $3 }' " + str(TOP) + ' | grep -cx "P" ',shell=True))
	NFRS = int(subprocess.check_output('grep -c "TIMESTEP" '+str(G96),shell=True))
	NATOMS = int(subprocess.check_output('grep -c "ATOM" '+str(TOP),shell=True))

	return N, NFRS, NATOMS
	#array = np.array([protname,N,NFRS,NATOMS],dtype=str)
	#np.savetxt("protname.txt",array.T,fmt="%s")

def get_chain_info(PROTNAME,TOP):
	import subprocess
	import numpy as np
	from string import ascii_uppercase

	subprocess.call("grep '^ATOM' " + TOP +" > tmp",shell=True)

	data = []
	with open("tmp","r+") as f:
		for line in f:
			data.append(line)


	chain_list = []
	atype_list = []
	for i in range(len(data)):
		chain_list.append(data[i][21])
		atype_list.append(data[i].split()[2])
	chain_list = np.array(chain_list)

	sum_list = []
	for c in ascii_uppercase:
		sum_list.append((chain_list == str(c)).sum())
		
	sum_list = np.array(sum_list)

	for i in range(len(sum_list)):
		if sum_list[i] != 0:
				sum_list[i] = int(sum_list[i])
			
	nres_list = np.zeros(len(ascii_uppercase))         
	for i,chain in enumerate(chain_list):
		if atype_list[i] == "P":
			for num,c in enumerate(ascii_uppercase):
				if chain == str(c):
					nres_list[num] += 1
					continue
					
	for i in range(len(nres_list)):
		if sum_list[i] != 0:
				nres_list[i] = int(nres_list[i])

	#with open('nmol.dat','w+') as f:
	#	f.write(str((sum_list != 0).sum())+"\n")
	#f.close()

	subprocess.call("rm -rfv tmp",shell=True)

	return (sum_list != 0).sum(), nres_list[nres_list != 0], sum_list[sum_list != 0]

'''def gen_protinfo(PROTNAME,PDB,TOP):
		import subprocess
		import numpy as np
		protname = str(PROTNAME)
		N = int(subprocess.check_output('grep -c "CA" '+str(TOP),shell=True))
		NFRS = int(subprocess.check_output('grep -c "MODEL" '+str(PDB),shell=True)-1)
		NATOMS = int(subprocess.check_output('grep -c "ATOM" '+str(TOP),shell=True))
		array = np.array([protname,N,NFRS,NATOMS],dtype=str)
		np.savetxt("protname.txt",array.T,fmt="%s")'''

'''def convert_traj(traj):
	import numpy as np
	traj = np.loadtxt(str(traj))
	#Save the trajectory in units of nm
	np.save('unformatted_traj.npy',traj/10.0)'''

def convert_traj(G96):
	import numpy as np
	import subprocess
	import platform

	#Check which system the code is run on. If Linux, I can use 'sed'. For Darwin (macOS) I need 'gsed'.
	if platform.system() == 'Linux':
		status = subprocess.call("sed '/BOX/, +1 d' " + str(G96) + " | sed '/TITLE/, +1 d' | awk 'NF==3' > tmp",shell=True)
	elif platform.system() == 'Darwin':
		status = subprocess.call("gsed '/BOX/, +1 d' " + str(G96) + " | gsed '/TITLE/, +1 d' | awk 'NF==3' > tmp",shell=True)
	else:
		raise OSError("System platform not recognized.")
	if status == 0:
		traj = np.loadtxt('tmp')

		#np.save('unformatted_traj.npy',traj)
		subprocess.call('rm -rfv tmp', shell=True)
		return traj
	else:
		raise OSError('''Something has gone incorrectly and the unformatted trajectory was not generated.
		Please check where the .g96 file is located and make sure the correct PATH is specified 
		in the call to this function.''')

def format_traj(traj, N, NFRS):
	import numpy as np

	ftraj = np.zeros((3*N, NFRS))
	for numba, k in enumerate(range(0, N*NFRS, N)):
		ftraj[::3, numba] = traj[k:k+N, 0]
		ftraj[1::3, numba] = traj[k:k+N, 1]
		ftraj[2::3, numba] = traj[k:k+N, 2]

	return ftraj

def Umatrix(traj, protname, N, nfrs, nmol, nca):
	import numpy as np

	nca_cumsum = np.cumsum(nca)
	print(protname, N, nfrs, nmol)

	rx = np.zeros((N,nfrs))
	ry = np.zeros((N,nfrs))
	rz = np.zeros((N,nfrs))
	lx = np.zeros((N - nmol, nfrs))
	ly = np.zeros((N - nmol, nfrs))
	lz = np.zeros((N - nmol, nfrs))
	Rinv = np.zeros((N,N))
	for numba, i in enumerate(range(0, 3*N, 3)):
		rx[numba, :] = traj[i, :]
		ry[numba, :] = traj[i + 1, :]
		rz[numba, :] = traj[i + 2, :]

	#Define bond vectors
	counter = 0
	mol_counter = 0
	i = 1
	while i <= N:
		if i == nca_cumsum[mol_counter]:
			i += 1
			mol_counter += 1
			if mol_counter == nmol: break
		else:
			lx[counter,:] = rx[i,:] - rx[i-1,:] 
			ly[counter,:] = ry[i,:] - ry[i-1,:]
			lz[counter,:] = rz[i,:] - rz[i-1,:]
			counter += 1
			i += 1

	#Calculate average inverse distances for the
	#hydrodynamic interaction matrix
	for i in range(N):
			for j in range(i,N):
				if i == j:
					Rinv[i,i] = np.nan
				else:
					Rinv[i,j] = (1/np.sqrt((rx[i,:] - rx[j,:])**2 + (ry[i,:] - ry[j,:])**2 + (rz[i,:] - rz[j,:])**2)).sum()
				Rinv[j,i] = Rinv[i,j]
	Rinv = Rinv/nfrs
	lavm = np.zeros(N - nmol)
	lavmsq = np.zeros(N - nmol)
	avdot = np.zeros(N - nmol)
	avblsq = 0
	avbl = 0

	avgx = lx.mean(1)
	avgy = ly.mean(1)
	avgz = lz.mean(1)

	for i in range(N - nmol):
		for k in range(nfrs):
			dummy = lx[i,k]**2 + ly[i,k]**2 + lz[i,k]**2
			lavm[i] += np.sqrt(dummy)
			lavmsq[i] += dummy
			avblsq += dummy 
		avbl += lavm[i]
		avdot[i] = avgx[i]**2 + avgy[i]**2 + avgz[i]**2
	lavm = lavm/nfrs
	lavmsq = lavmsq/nfrs
	avblsq = (avblsq/(N - nmol))/nfrs
	avbl = (avbl/(N - nmol))/nfrs

	print(avbl,avblsq)

	Umat = np.zeros((N - nmol,N - nmol))
	for i in range(N - nmol):
		for j in range(N - nmol):
			Umat[i,j] = (np.dot(lx[i,:],lx[j,:]) + np.dot(ly[i,:],ly[j,:]) + np.dot(lz[i,:],lz[j,:]))/(lavm[i]*lavm[j]*nfrs)

	#np.save('Umatrix.npy',Umat)
	#np.savetxt('Umatrix',np.insert(np.ravel(Umat),N - nmol,0))
	#np.save('Rinv.npy',Rinv)
	#np.savetxt('Rij',np.ravel(Rinv))
	#np.savetxt('length',lavm)
	#np.savetxt('lengthsq',lavmsq)
	#np.savetxt('avldot.dat',avdot)
	#np.savetxt('avbl',np.array([avbl]))
	#np.savetxt('avblsq',np.array([avblsq]))

	return Umat, Rinv, avbl, avblsq

def fric_calc(TOP, protname, N, nfrs, natoms, avblsq, T, intv = 2.71828, viscosity = 1e-3, fd20 = 0.0, path_to_resarea = './'):
	import numpy as np
	import sys
	import os
	import subprocess

	#TOP = str(sys.argv[1])
	pi = np.pi

	print(protname, N, nfrs, natoms)

	#Calculate the Miller radius per bead
	mradlist = []
	with open(TOP) as f:
		for line in f:
			if line[0:4] != 'ATOM':
				#print(line)
				pass
			elif line[0:4] == 'ATOM' and line.split()[2] == "P":
				dummy = line.split()

				#Really horrendous and vestigial; probably smoother
				#to make a dictionary with the nucleotide names plus their
				#assoicated radii

				#Deoxynucleotide volumes from "Standard atomic volumes in double-stranded DNA andpacking in protein–DNA interfaces" by Nadassy et al. NAR, 2001.
				#Ribonucleotide volumes from "Calculation of Standard Atomic Volumes for RNA andComparison with Proteins: RNA is Packed More Tightly" by Voss and Gersteain, JMB, 2005.
				if dummy[3] == "DA5" or  dummy[3] == "DA3" or  dummy[3] == "DA" : mradlist.append(np.cbrt((3*302.0) / (4*pi)))
				elif dummy[3] == "DC5" or  dummy[3] == "DC3" or  dummy[3] == "DC" : mradlist.append(np.cbrt((3*283.7) / (4*pi)))
				elif dummy[3] == "DG5" or  dummy[3] == "DG3" or  dummy[3] == "DG" : mradlist.append(np.cbrt((3*314.7) / (4*pi)))
				elif dummy[3] == "DT5" or  dummy[3] == "DT3" or  dummy[3] == "DT" : mradlist.append(np.cbrt((3*301.7) / (4*pi)))
				elif dummy[3] == "A5" or  dummy[3] == "A3" or  dummy[3] == "A" : mradlist.append(np.cbrt((3*315.3) / (4*pi)))
				elif dummy[3] == "C5" or  dummy[3] == "C3" or  dummy[3] == "C" : mradlist.append(np.cbrt((3*291.1) / (4*pi)))
				elif dummy[3] == "G5" or  dummy[3] == "G3" or  dummy[3] == "G" : mradlist.append(np.cbrt((3*322.0) / (4*pi)))
				elif dummy[3] == "U5" or  dummy[3] == "U3" or  dummy[3] == "U": mradlist.append(np.cbrt((3*286.9) / (4*pi)))
	#array = np.array(mradlist,dtype=str)
	#np.savetxt('mrad.dat',np.array(mradlist).T,fmt='%s')

	#Calculate the average solvent-exposed surface area per bead   gmx_mpi sasa -f after_rot.xtc -s dsAXA_pro.tpr -or resarea.xvg -dt 1000
    
    #gmx_mpi sasa -f after_rot_sim0.xtc -s dsAXA.pdb -n index.ndx -or resarea_DNA_base2.xvg -dt 1000 -surface 1 -output 3
    
    
	if os.path.exists(path_to_resarea + "resarea.xvg"):
		pass
	else:
		raise FileNotFoundError('''I can't find the resarea.xvg file containing the solvent-exposed surface area of each residue.
								Please either run the process.sh file, if you have not already done so, and move the resarea.xvg 
								file into the current working directory.''')

	resarea = []
	with open(path_to_resarea + 'resarea.xvg') as f:
		for line in f:
			if line[0] == '#' or line[0] == '@':
				pass
			else:
				resarea.append(float(line.split()[1]))

	rad = []
	for area in resarea:
		rad.append(((area/(4*np.pi))**0.5)*10)

	#np.savetxt('avresrad',np.array(rad),fmt="%f")
	fratio = (np.array(rad).sum()/N)/10
	print('fratio: ', fratio)

	#np.savetxt('fratio',np.array([fratio]))    

	#Calculate the friction coefficients

	kB = 1.38066E-23
	print('Temperature (K): ',T)
	print('Internal viscosity factor: ',intv)

	#Use NIST formula for viscosity -- good NEAR room temperature and physiological.
	#Won't work higher than, say, 360 K.

	if viscosity == 0:
		print('No viscosity given. Using the NIST formula, which is only valid for physiological conditions,\n')
		print('i.e. between about 273 and 310 K.')
		viscosity = (.2131590-1.96290E-3*T+(.00246411*T)**2+(-.0018462*T)**3)

	print("Viscosity (Pa s): ", viscosity)
	print("fd20", fd20)

	rv = np.array(mradlist)
	rw = np.array(rad)
	rp = np.zeros(N)
	friw = np.zeros(N)
	fri = np.zeros(N)
	friwt = 0
	frit = 0
	for i in range(N):
		if rw[i] < rv[i]: 
			rp[i] = (rv[i]**2 - rw[i]**2)**0.5
		else:
			rp[i] = 0

		friw[i] = 6.0*pi*(rw[i]/10)*viscosity
		fri[i] = 6.0*pi*(rp[i]/10)*(intv*viscosity) + 6.0*pi*(rw[i]/10)*viscosity
		friwt += friw[i]
		frit += fri[i]

	avfr = frit/float(N)
	avfrw = friwt/float(N)
	#np.savetxt('avfr',np.array([avfr*1.0E-9]))

	#avblsq = float(np.loadtxt('avblsq'))
	sigma = (3*kB*T*1E15)/(avblsq*avfr)

	#with open('sigma','w') as f:
	#	f.write('sigma, 1/ps\n')
	#	f.write(str(sigma)+'\n')

	#with open('sigma.dat','w') as f:
	#	f.write(str(sigma)+'\n')

	fric = np.zeros((N+1, 2))

	fric[0,0] = avfrw
	fric[0,1] = avfr
	for i in range(N):
		fric[i+1,:] = np.column_stack([friw[i],fri[i]])

	#np.savetxt('fric',fric)

	return fratio, sigma, fric, avfr

def LUI_calc(protname, N, nfrs, nmol, nca, U, fratio, avblsq, sigma, fric, Rinv, T):
	import numpy as np

	nca_cumsum = np.cumsum(nca)
	print(protname, N, nfrs, nmol)

	avfr = fric[0,1]

	M = np.zeros((N,N))

	idx = 0
	counter = 1
	for i in range(N):
		if i == idx:
			for j in range(N):
				M[i,j] = 1/N
			idx = nca_cumsum[counter - 1]
			counter += 1
		else:
			for j in range(N):
				if i == j + 1:
					M[i,j] = -1
				elif i == j:
					M[i,j] = 1

	a = np.zeros((N - nmol,N))
	counter = 0
	for i in range(N):
		if ((i == 0) or (np.where(nca_cumsum == i)[0].size == 1)):
			#print('here')
			pass
		else:
			a[counter,:] = M[i,:]
			counter += 1

	H = np.zeros((N,N))

	for i in range(N):
		for j in range(i,N):
			if i == j:
				H[i,i] = avfr/fric[i+1,1]
			else:
				H[i,j] = fratio*Rinv[i,j]
				H[j,i] = H[i,j]

	L = np.matmul(a,np.matmul(H,a.T))

	LINV = np.linalg.inv(L)
	L = np.matmul(a,np.matmul(H,a.T))
	LINV = np.linalg.inv(L)
	UINV = U #np.load("Umatrix.npy")
	UILI = np.matmul(UINV,LINV)

	#Eigendecomposition of UILI to find eigenvalues and eigenvectors 
	eigval,Q = np.linalg.eig(UILI)
	eigval = 1/eigval
	QINV = np.linalg.inv(Q)

	#Implement fratio loop to avoid negative eigenvalues of LU matrix
	while (np.any(1/eigval <= -1e-9) and fratio >= 0.0):
		fratio -= 0.01
		#print(fratio,(1/eigval)[np.where(1/eigval <= -1e-9)])
		H = np.zeros((N,N))
		for i in range(N):
			for j in range(i,N):
				if i == j:
					H[i,i] = avfr/fric[i+1,1]
				else:
					H[i,j] = fratio*Rinv[i,j]
					H[j,i] = H[i,j]

		L = np.matmul(a,np.matmul(H,a.T))

		LINV = np.linalg.inv(L)
		L = np.matmul(a,np.matmul(H,a.T))
		LINV = np.linalg.inv(L)
		UINV = U #np.load("Umatrix.npy")
		UILI = np.matmul(UINV,LINV)

		#Eigendecomposition of UILI to find eigenvalues and eigenvectors 
		eigval,Q = np.linalg.eig(UILI)
		eigval = 1/eigval
		QINV = np.linalg.inv(Q)

	perm = np.argsort(np.abs(eigval))

	Q_sorted = np.copy(Q)[:,perm]
	QINV_sorted = np.linalg.inv(Q_sorted)

	eigval_sorted = abs(eigval)[perm]

	mu = 1/(np.diag(np.matmul(QINV_sorted,np.matmul(UINV,QINV_sorted.T))))

	#np.save("UILImatrix",UILI)
	#np.savetxt('UILImatrix',np.ravel(UILI))
	#np.save('Lmatrix',L)
	#np.save("Hmatrix",H)
	#np.save("Qmatrix.npy",Q_sorted)
	#np.savetxt("Qmatrix",np.ravel(Q_sorted))
	#np.save("QINVmatrix.npy",QINV_sorted)
	#np.savetxt("QINVmatrix",np.ravel(QINV_sorted))
	#np.save("lambda_eig.npy",eigval_sorted)
	#np.savetxt("lambda_eig",eigval_sorted)
	#np.save("mu_eig.npy",mu)
	#np.savetxt("mu_eig",mu)

	return UILI, H, Q_sorted, QINV_sorted, eigval_sorted, mu

def mode_mad(traj, protname, N, nfrs, nmol, nca, Q, QINV, T, nmodes = 10):
	import numpy as np
	import matplotlib.pyplot as plt
	import os
	import subprocess
	import physt

	pi = np.pi

	print(protname, N, nfrs, nmol)

	nca_cumsum = np.cumsum(nca)

	rx = np.zeros((N, nfrs))
	ry = np.zeros((N, nfrs))
	rz = np.zeros((N, nfrs))
	lx = np.zeros((N - nmol, nfrs))
	ly = np.zeros((N - nmol, nfrs))
	lz = np.zeros((N - nmol, nfrs))

	#Already saved the formatted trajectory
	for numba, i in enumerate(range(0, 3*N, 3)):
		rx[numba, :] = traj[i, :]
		ry[numba, :] = traj[i + 1, :]
		rz[numba, :] = traj[i + 2, :]

	#Define bond vectors
	counter = 0
	mol_counter = 0
	i = 1
	while i <= N:
		if i == nca_cumsum[mol_counter]:
			i += 1
			mol_counter += 1
			if mol_counter == nmol: break
		else:
			lx[counter,:] = rx[i,:] - rx[i-1,:] 
			ly[counter,:] = ry[i,:] - ry[i-1,:]
			lz[counter,:] = rz[i,:] - rz[i-1,:]
			counter += 1
			i += 1

	xix = np.matmul(QINV,lx)
	xiy = np.matmul(QINV,ly)
	xiz = np.matmul(QINV,lz)
	xim = np.sqrt(xix**2 + xiy**2 + xiz**2)
	xi = xix + xiy +xiz

	theta = np.arccos(xiz/xim)
	phi = np.arctan(xiy/xix)

	if nmodes > N - nmol:
		raise ValueError('''Specified more modes than are generated by the analysis.\n
			The isotropic LE4PD generates N - nmol modes, with N the number of residues in the protein.''')
	else:

		for a in range(N - nmol):
			for k in range(nfrs):
				if xix[a,k] <= 0.0:
					phi[a,k] += pi
				if phi[a,k] <= 0.0:
					phi[a,k] += 2*pi

		#theta = np.rad2deg(theta)
		#phi = np.rad2deg(phi)

		#Make histogram
		kT = 0.00198*T
		fmadlist = []
		feslist = []
		for a in range(N-nmol):
			x=xim[a,:]*np.sin(theta[a,:])*np.cos(phi[a,:])
			y=xim[a,:]*np.sin(theta[a,:])*np.sin(phi[a,:])
			z=xim[a,:]*np.cos(theta[a,:])
			if int(physt.__version__.split('.')[1]) <= 3:
				h=physt.special.spherical_histogram(np.column_stack([x,y,z]),theta_bins=50,phi_bins=50,radial_bins=1)
			else:
				h=physt.special_histograms.spherical(np.column_stack([x,y,z]),theta_bins=50,phi_bins=50,radial_bins=1)
			his=(h.densities[0]/h.densities[0].sum()).T
			fes = -kT*np.log(his)
			#fes -= fes.min()
			femax = -kT*np.log(1/nfrs)
			for j in range(fes.shape[0]):
				for k in range(fes.shape[1]):
					if np.isinf(fes[j,k]) == True:
						fes[j,k] = femax
						
			#xx,yy=np.meshgrid(np.linspace(0,180,his.shape[1]),np.linspace(0,360,his.shape[0]))
			#im=plt.contourf(xx,yy,fes,levels=np.linspace(fes.min(),fes.max(),25),cmap='gnuplot')
			#cbar=plt.colorbar(im)
			#cbar.set_label(r'Free Energy $(k_BT)$')
			#plt.contour(xx,yy,fes,levels=np.linspace(fes.min(),fes.max(),25),colors='k')
			#plt.xlabel(r'$\theta$ (deg)')
			#plt.ylabel(r'$\phi$ (deg)')
			#plt.savefig('fes'+str(a+1)+'.eps',dpi=300)
			#plt.show()
			#plt.close()
			dummy = list(np.ravel(fes))
			#Cut-off for outliers 
			cut = -kT*np.log(2/nfrs)
			for num,i in enumerate(dummy):
				if i >= cut: dummy.remove(i)
			fmad = np.median(abs(np.array(dummy) - fes.min()))

			fmadlist.append(fmad)
			feslist.append(fes)

		#Make a directory for this mode-dependent analysis
		if os.path.exists('mode_analysis'):
			pass
		else:
			os.mkdir('mode_analysis/')
		path = 'mode_analysis/'

		for a in range(nmodes):
			np.save(path + 'fes'+str(a+1)+'.npy', feslist[a])
			np.save(path + 'theta_phi_'+str(a+1)+'.npy', np.column_stack([np.rad2deg(theta[a,:]),np.rad2deg(phi[a,:]), xim[a,:]]))
			np.save(path + "xi_"+str(a+1)+'.npy', xi[a,:])
		np.savetxt(path + 'barriers_kcal.dat', np.array(fmadlist))
		np.savetxt(path + 'barriers.dat', np.array(fmadlist)/kT)

		return np.array(fmadlist)/kT, xi, np.column_stack([np.rad2deg(theta[a,:]),np.rad2deg(phi[a,:]), xim[a,:]])

def tau_convert(eigvals, sigma, bar, T):
	import numpy as np
	import os
	
	if os.path.exists('mode_analysis'):
		pass
	else:
		os.mkdir('mode_analysis/')
	path = 'mode_analysis/'

	tau = (eigvals*sigma)**-1
	tau_scaled = tau*np.exp(bar)
	np.savetxt(path + 'tau.dat', np.column_stack([np.arange(1,len(tau)+1),tau]))
	np.savetxt(path + 'tau_scaled.dat', np.column_stack([np.arange(1,len(tau)+1),tau_scaled]))

	return tau, tau_scaled

def LML(Qmatrix, avbl, mu):
	import numpy as np
	import os
	
	LML = np.zeros_like(Qmatrix)
	for a in range(Qmatrix.shape[0]):
		for i in range(Qmatrix.shape[1]):
			LML[i,a] = (Qmatrix[i,a]**2*avbl**2)/mu[a]
	LML = np.sqrt(LML)

	if os.path.exists('mode_analysis'):
		pass
	else:
		os.mkdir('mode_analysis/')
	path = 'mode_analysis/'

	np.save(path + "LML.npy", LML)

	return LML
