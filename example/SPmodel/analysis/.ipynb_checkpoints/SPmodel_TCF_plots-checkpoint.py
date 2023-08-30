import numpy as np
import matplotlib.pyplot as plt
import LE4DNA as le4dna
import LE4DNA_extras as le4dna_extras
import subprocess
import os

CWD = os.getcwd()
print(f'cwd: {CWD}')

path_to_data = CWD + '/tcfdata/'
path_to_npy  = CWD + '/numpy_arrays/'
path_to_imgs = CWD + '/plots/'

simula_path  = path_to_data + 'm1CAsimint_'
theory_path  = path_to_data + 'm1_'

protname_data = subprocess.check_output('head protname.txt',shell=True).decode('ascii').strip()
nbonds = int(protname_data.split('\n')[1]) - 2
nfrs   = int(protname_data.split('\n')[2])
print(f'nbonds = {nbonds}')
print(protname_data)

if os.path.exists(path_to_npy + 'tcf_simula.npy') and os.path.exists(path_to_npy + 'tcf_theory.npy'): 
    tcf_simula = np.load(path_to_npy + 'tcf_simula.npy')
    tcf_theory = np.load(path_to_npy + 'tcf_theory.npy')
else:
    #npy = le4dna_extras.autocorr(nbonds, save)
    tcf_simula = le4dna_extras.autocorr(nbonds, simula_path)
    tcf_theory = le4dna_extras.autocorr(nbonds, theory_path)
    np.save(path_to_npy + 'tcf_simula.npy', tcf_simula)
    np.save(path_to_npy + 'tcf_theory.npy', tcf_theory)
    print(f'TCFs of bond autocorrelation saved as .npy arrays in {path_to_npy}')
    
# The loop below will plot TCFs for each bond on a linear scale and log scale.
# The plots are then saved individually, and as a combined document in path_to_imgs
# plotting params:
bonds = 2 # nbonds
excld = 100
xmin  = 0
xmax  = (nfrs-excld)*0.2
ymin  = 0.75
ymax  = 1.0



fig, ax = plt.subplots(nrows=bonds, ncols=2, figsize=(12,4*bonds))
for bond in range(0,bonds):
    ax[bond,0].plot(tcf_simula[bond,0],tcf_simula[bond,1],'k',label='SP sim.',linewidth=1.5)
    ax[bond,0].plot(tcf_theory[bond,0],tcf_theory[bond,1],'b',label='SP theory',linewidth=2.0)
    ax[bond,0].set_xlim(xmin,xmax)
    #ax[bond,0].set_ylim(0.85,1.0)
    ax[bond,0].set_xlabel('t (ps)',size=12)
    ax[bond,0].set_ylabel('$M_{1,i}(t)$',size=14)
    ax[bond,0].set_title('Bond ' + str(bond+1),size=14)
    ax[bond,0].legend()
    ax[bond,0].set_facecolor('xkcd:pale blue')
    ax[bond,0].grid(color='k', linestyle='dashed',alpha=0.15)
    ax[bond,0].set_axisbelow(True)
    #
    ax[bond,1]..plot(tcf_simula[bond,0],tcf_simula[bond,1],'k',label='SP sim.',linewidth=1.5)
    ax[bond,1]..plot(tcf_theory[bond,0],tcf_theory[bond,1],'b',label='SP theory',linewidth=2.0)
    ax[bond,1].set_xlim(xmin + 1, xmax)
    #ax[bond,1].set_ylim(0.85,1.0)
    ax[bond,1].set_xscale('log')
    ax[bond,1].set_xlabel('t (ps)',size=12)
    ax[bond,1].set_ylabel('$M_{1,i}(t)$',size=14)
    ax[bond,1].set_title('Bond ' + str(bond+1),size=14)
    ax[bond,1].legend()
    ax[bond,1].set_facecolor('xkcd:pale blue')
    ax[bond,1].grid(color='k', linestyle='dashed',alpha=0.15)
    ax[bond,1].set_axisbelow(True)
    
plt.tight_layout()
plt.show()
plt.savefig(path_to_imgs + 'SP_model_TCF_plots.pdf')
plt.close()
