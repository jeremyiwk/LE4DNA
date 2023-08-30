#
# ****  
# ****  
# ****  
# 
# 
# 
# 
# 
# ****  
# ****  
# ****  
#

import numpy as np
import LE4DNA as le4pd
import matplotlib.pyplot as plt
import physt
import LE4DNA_extras as le4pd_extras

TOP = 'SP_dsAXA.pdb'
SYS = 'dsAXA'
P_G96 = 'P_traj.g96'
S_G96 = 'S_traj.g96'

T = 300
N_S = 624
N_P = 230
N = 94
NATOMS = 94

NFRS_S = 10001
NFRS_P = 10001
NFRS   = 10001

nmol = 2
nres_list = np.array([47.,47.])
sum_list = np.array([])

Smassarr = np.array([ 12.01000, 1.00800, 1.00800, 12.01000, 1.00800, 16.00000 , 12.01000 , 1.00800, 12.01000, 1.00800, 12.01000, 1.00800, 1.00800])
Pmassarr = np.array([16.00000, 30.97000, 16.00000, 16.00000, 16.00000])

P_traj = np.load('P_traj.npy')
S_traj = np.load('S_traj.npy')
# 
P_ftraj = le4pd.format_traj(P_traj, N_P, NFRS_P)
S_ftraj = le4pd.format_traj(S_traj, N_S, NFRS_S)
# 
np.save('P_ftraj.npy',P_ftraj)
np.save('S_ftraj.npy',S_ftraj)
# 
P_com_ftraj = le4pd_extras.com_traj(P_ftraj, NFRS_P, 46, Pmassarr)
S_com_ftraj = le4pd_extras.com_traj(S_ftraj, NFRS_S, 48, Smassarr)
# 
np.save('P_com_ftraj.npy',P_com_ftraj)
np.save('S_com_ftraj.npy',S_com_ftraj)
# 
SP_com_ftraj = np.zeros((3*(46+48), NFRS_S))
# 
H = int(3*(46+48)/2)
SH = int(3*(48)/2)
PH = int(3*(46)/2)
# 
SP_com_ftraj[0:H:6]=S_com_ftraj[0:SH:3]   #sugar x coords, strand 1
SP_com_ftraj[1:H:6]=S_com_ftraj[1:SH:3]   #sugar y coords, strand 1
SP_com_ftraj[2:H:6]=S_com_ftraj[2:SH:3]   #sugar z coords, strand 1
###################################################################
SP_com_ftraj[H::6]=S_com_ftraj[SH::3]     #sugar x coords, strand 2
SP_com_ftraj[H+1::6]=S_com_ftraj[SH+1::3] #sugar y coords, strand 2
SP_com_ftraj[H+2::6]=S_com_ftraj[SH+2::3] #sugar z coords, strand 2
###################################################################
SP_com_ftraj[3:H:6]=P_com_ftraj[0:PH:3]   #phos x coords, strand 1
SP_com_ftraj[4:H:6]=P_com_ftraj[1:PH:3]   #phos y coords, strand 1
SP_com_ftraj[5:H:6]=P_com_ftraj[2:PH:3]   #phos z coords, strand 1
###################################################################
SP_com_ftraj[H+3::6]=P_com_ftraj[PH::3]   #phos x coords, strand 2
SP_com_ftraj[H+4::6]=P_com_ftraj[PH+1::3] #phos y coords, strand 2
SP_com_ftraj[H+5::6]=P_com_ftraj[PH+2::3] #phos z coords, strand 2
# 
np.save('SP_com_ftraj.npy',SP_com_ftraj)

#SP_com_ftraj = np.load('SP_com_ftraj.npy')

print('Done creating, weaving, and saving SP model trajectory!')

le4pd_extras.save_single_particle_traj('tcfrun/', SP_com_ftraj, 94)

Umat, Rinv, avbl, avblsq =  le4pd.Umatrix(SP_com_ftraj, SYS, N, NFRS, nmol, nres_list)

np.save('Umat.npy',Umat)
np.save('Rinv.npy',Rinv)
np.save('avbl.npy',avbl)
np.save('avblsq.npy',avblsq)

#Umat, Rinv, avbl, avblsq = np.load('Umat.npy'), np.load('Rinv.npy'), np.load('avbl.npy'), np.load('avblsq.npy')

print('Done computing Umat!')
print('avblsq = ' + str(avblsq))

fratio, sigma, fric, avfr = le4pd_extras.fric_calc_SPmodel('SP_model.pdb', SYS, N, NFRS, avblsq, T=300, intv = 2.71828, viscosity = 3.1e-4, fd20 = 0.0, path_to_resarea = './', cgmodel = 'SP')
np.save('fratio.npy',fratio)
np.save('sigma.npy',sigma)
np.save('fric.npy',fric)
np.save('avfr.npy',avfr)

print('Done computing fric!')
print('sigma = ' + str(sigma))

UILI, H, Q_sorted, QINV_sorted, eigval_sorted, mu = le4pd.LUI_calc(SYS, N, NFRS, nmol, nres_list, Umat, fratio, avblsq, sigma, fric, Rinv, T)
np.save('UILI.npy',UILI)
np.save('H.npy',H)
np.save('Q_sorted.npy',Q_sorted)
np.save('QINV_sorted.npy',QINV_sorted)
np.save('eigval_sorted.npy',eigval_sorted)
np.save('mu.npy',mu)

print('Done diagonalizing!')

barlist, xi, dummy = le4pd.mode_mad(SP_com_ftraj, SYS, N, NFRS, nmol, nres_list, Q_sorted, QINV_sorted, T, nmodes = 10)

np.save('barlist.npy',barlist)
np.save('xi.npy',xi)
np.save('dummy.npy',dummy)

print('Done with mode mad!')

tau, tau_scaled = le4pd.tau_convert(eigval_sorted, sigma, barlist, T)

np.save('tau.npy',tau)
np.save('tau_scaled.npy',tau_scaled)

LML = le4pd.LML(Q_sorted, avbl, mu)

np.savetxt('Qmatrix', np.ravel(Q_sorted))
np.savetxt('QINVmatrix', np.ravel(QINV_sorted))
np.savetxt('lambda_eig', np.ravel(eigval_sorted))
np.savetxt('mu_eig', np.ravel(mu))
np.savetxt('avblsq.dat', np.array([avblsq]))
np.savetxt('sigma.dat', np.array([sigma]))
      
print('Done!')
