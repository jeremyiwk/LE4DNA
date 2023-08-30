# coding: utf-8

def com_traj(ftraj, nfrs, nmon, massarr):

    '''converts a formatted trajectory of monomers from the output of
    LE4DNA.format_traj() into a trajectory of the center of mass of each monomer.

    args:

    ftraj (numpy array) : output of LE4DNA.format_traj(), trajectory of monomers
    nfrs (int) : np.shape(ftraj)[1], number of frames in trajectory
    nmon (int) : np.shape(ftraj)[0]/(3*len(massarr)), number of monomers in molecule
    massarr (numpy array) : masses of atoms in each monomer. i.e. massarr[i] = mass
    of ith atom in the monomer.
    '''

    import numpy as np

    monsize = len(massarr)

    totalM = massarr.sum()

    com_traj = np.zeros( ( 3*nmon , nfrs ) )

    for i in range(0, 3*nmon, 3):

        for j in range(monsize):

            com_traj[i] += massarr[j]*ftraj[3*j + monsize*i] #xij

            com_traj[i+1] += massarr[j]*ftraj[3*j + monsize*i + 1] #yij

            com_traj[i+2] += massarr[j]*ftraj[3*j + monsize*i + 2] #zij

    com_traj *= 1/totalM

    return com_traj

def save_single_particle_traj( filename, ftraj, nparticles):

    import numpy as np

    if int(np.shape(ftraj)[0]/3) != nparticles:
        print("Wrong number of particles or wrong trajectory.")

    for i in range(nparticles):

        np.savetxt(filename + f'particle_{i+1}', ftraj[3*i:3*i+3].T)

def fric_calc_SPmodel(TOP, protname, N, nfrs, avblsq, T, intv = 2.71828, viscosity = 3.1e-4, fd20 = 0.0, path_to_resarea = './', cgmodel = 'SP'):
    
    # !!!! TOP must contain only coarse grained sites. !!!!
    
    import numpy as np
    import sys
    import os
    import subprocess
    
    pi = np.pi
    
    print(protname, cgmodel, N, nfrs)
    
    if cgmodel != 'SP':
        raise ValueError('''Available coarse grain models are: SP.''')
        
    # All volumes given are pulled from the papers below.   
    #Deoxynucleotide volumes from "Standard atomic volumes in double-stranded DNA andpacking in protein–DNA interfaces" by Nadassy et al. NAR, 2001.
    #Ribonucleotide volumes from "Calculation of Standard Atomic Volumes for RNA andComparison with Proteins: RNA is Packed More Tightly" by Voss and Gersteain, JMB, 2005.
    mrad_dict = { "DA5" : np.cbrt((3*310.9) / (4*pi)) ,
                  "DA3" : np.cbrt((3*310.9) / (4*pi)) ,
                  "DA"  : np.cbrt((3*310.9) / (4*pi)) ,
                  "DC5" : np.cbrt((3*288.0) / (4*pi)) ,
                  "DC3" : np.cbrt((3*288.0) / (4*pi)) ,
                  "DC"  : np.cbrt((3*288.0) / (4*pi)) ,
                  "DG5" : np.cbrt((3*318.6) / (4*pi)) ,
                  "DG3" : np.cbrt((3*318.6) / (4*pi)) ,
                  "DG"  : np.cbrt((3*318.6) / (4*pi)) ,
                  "DT5" : np.cbrt((3*307.4) / (4*pi)) ,
                  "DT3" : np.cbrt((3*307.4) / (4*pi)) ,
                  "DT"  : np.cbrt((3*307.4) / (4*pi)) ,
                  "A5"  : np.cbrt((3*315.3) / (4*pi)) ,
                  "A3"  : np.cbrt((3*315.3) / (4*pi)) ,
                  "A"   : np.cbrt((3*315.3) / (4*pi)) ,
                  "C5"  : np.cbrt((3*291.1) / (4*pi)) ,
                  "C3"  : np.cbrt((3*291.1) / (4*pi)) ,
                  "C"   : np.cbrt((3*291.1) / (4*pi)) ,
                  "G5"  : np.cbrt((3*322.0) / (4*pi)) ,
                  "G3"  : np.cbrt((3*322.0) / (4*pi)) ,
                  "G"   : np.cbrt((3*322.0) / (4*pi)) ,
                  "U5"  : np.cbrt((3*286.9) / (4*pi)) ,
                  "U3"  : np.cbrt((3*286.9) / (4*pi)) ,
                  "U"   : np.cbrt((3*286.9) / (4*pi))  }
    
    # sugar and base volumes, excluding volumes of atoms in phosphates (P, O1P, O2P, O3', O5')
    # A = 243.7, C = 220.8, G = 251.4, T = 240.2
    
    #volume of phosphate (P + O1P + O2P + O3' + O5')
    Prad_dict = { "P" : np.cbrt((3*67.2) / (4*pi)) }
    
    #volume of sugar+base
    Srad_dict = { "DA5" : np.cbrt((3*243.7) / (4*pi)) , 
                  "DA3" : np.cbrt((3*243.7) / (4*pi)) ,
                  "DA"  : np.cbrt((3*243.7) / (4*pi)) ,
                  "DC5" : np.cbrt((3*220.8) / (4*pi)) ,
                  "DC3" : np.cbrt((3*220.8) / (4*pi)) , 
                  "DC"  : np.cbrt((3*220.8) / (4*pi)) ,
                  "DG5" : np.cbrt((3*251.4) / (4*pi)) ,
                  "DG3" : np.cbrt((3*251.4) / (4*pi)) ,
                  "DG"  : np.cbrt((3*251.4) / (4*pi)) ,
                  "DT5" : np.cbrt((3*240.2) / (4*pi)) ,
                  "DT3" : np.cbrt((3*240.2) / (4*pi)) ,
                  "DT"  : np.cbrt((3*240.2) / (4*pi))  }
    
                 
                 
    #Calculate the Miller radius per bead
    mradlist = []
    #P_mradlist = []
    #S_mradlist = []
    with open(TOP) as f:
        for line in f:
            if line[0:4] != 'ATOM':
                pass
            elif line[0:4] == 'ATOM':   # and line.split()[2] == "P":
                dummy = line.split()
                if dummy[2] == "P":
                    mradlist.append( Prad_dict[ "P" ] )
                elif dummy[2] == "S":
                    mradlist.append( Srad_dict[ dummy[3] ])
    
    mradlist = np.array(mradlist)
            
    #file should look like: "./SP_resarea.xvg"       
    if os.path.exists(path_to_resarea + cgmodel + "_resarea.xvg"):
        pass
    else:
        raise FileNotFoundError('''I can't find the resarea.xvg file containing the solvent-exposed surface area of each residue.
                                Please either run the process.sh file, if you have not already done so, and move the resarea.xvg 
                                file into the current working directory.''')
    
    # To compute the solvent exposed surface area, run the following GROMACS command:
    # gmx_mpi sasa -f after_rot_sim0.xtc -s dsAXA.pdb -n SP_index.ndx -or SP_resarea.xvg -dt 1000 -surface 1 -output 2 3
    # This assumes that SP_index.ndx contains the atom indices for the sugar+base and the phosphates in selections
    # 2 and 3, respectively. 
    
    #determine solvent exposed surface areas of each group from the .xvg file
    # this assumes lines in .xvg have the following format:
    # res  resarea  std      Sarea    std      Parea    std
    # 1    2.211    0.179    1.118    0.102    0.000    0.000
    # 2    1.942    0.195    0.603    0.085    0.842    0.076 ...
    S_resarea = []
    P_resarea = []
    with open(path_to_resarea + cgmodel + "_resarea.xvg") as f:
        for line in f:
            if line[0] == '#' or line[0] == '@':
                pass
            else:
                S_resarea.append(float(line.split()[3]))
                P_resarea.append(float(line.split()[5]))
                       
    S_resarea = np.array(S_resarea)
    P_resarea = np.array(P_resarea)
    P_resarea = P_resarea[ P_resarea != 0.0] #each res has a sugar+base but not all have a P.
    
    #build the list of SASA for each group in the same order as mrad_list.
    #the list slicing below yields SP_resarea = [S,P,...,P,S,S,P,...,P,S]
    SP_resarea = np.zeros(N)
    SP_resarea[0:int(N/2):2] = S_resarea[0: int(len(S_resarea)/2)]
    SP_resarea[int(N/2)::2] = S_resarea[int(len(S_resarea)/2):]
    
    SP_resarea[1:int(N/2):2] = P_resarea[0: int(len(P_resarea)/2)]
    SP_resarea[int(N/2)+1:-1:2] = P_resarea[int(len(P_resarea)/2):]

    #create list of CG-site solvent exposed radii
    # S_rad = ((S_resarea/(4*pi))**0.5)*10
    # P_rad = ((P_resarea/(4*pi))**0.5)*10
    SP_rad = ((SP_resarea/(4*pi))**0.5)*10
            
    
    #np.savetxt('avresrad',np.array(rad),fmt="%f")
    fratio = (SP_rad.sum()/N)/10
    print('fratio: ', fratio)
    
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
    
    rv = mradlist    #np.array(mradlist)
    rw = SP_rad    #np.array(rad)
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
        

def autocorr(nbonds, path):
    import numpy as np
    autocorr = []
    for num in range(nbonds):
        col1 = []
        col2 = []
        with open( path + str(num + 1)) as file:
            for line in file:
                line = ' '.join(line.split())
                line = line.split()
                col1.append(float(line[0]))
                col2.append(float(line[1]))
        file.close()
        autocorr.append( [col1, col2] )
    autocorr = np.array(autocorr)
    print(autocorr.shape)
    return autocorr

def sequence_data(pdb_filename, path = './'):

    import numpy as np
    from string import ascii_uppercase

    if '.pdb' in pdb_filename:
        file = path + pdb_filename
    else:
        file = path + pdb_filename + '.pdb'
    linelist = []
    with open(file) as pdb:
        for line in pdb:
            if line.split()[0] == 'ATOM':
                linelist.append(line.split()[1:6])
    linelist = np.array(linelist)
    chains = len(set(linelist[:,3]))
    chainindex = ascii_uppercase[0:chains]
    bases_by_chain = np.zeros(chains)
    for i in range(chains):
        bases_by_chain[i] = max( linelist[:,4][ linelist[:,3] == chainindex[i] ].astype('int32'))
    chain_sequence = []
    for i in range(chains):
        enum_list = np.unique( linelist[:,2:5:2][ linelist[:,3] == chainindex[i] ] , axis = 0 )
        for j in range(1, int(bases_by_chain[i]) + 1):
            chain_sequence.append( enum_list[:,0][ enum_list[:,1 ] == str(j) ][0] )


    return  chain_sequence, chains, bases_by_chain, linelist

def ftraj2pdb( ftraj_frame , atomlist, chain_sequence, ref_pdb, filename, path = './'):

    import numpy as np
    import platform
    import subprocess
    import os
    from string import ascii_uppercase

    linelist = []
    with open(ref_pdb) as pdb:
        for line in pdb:
            if line.split()[0] == 'ATOM':
                linelist.append(line.split()[1:6])
    linelist = np.array(linelist)

    chains = len(set(linelist[:,3]))
    chainindex = ascii_uppercase[0:chains]

    sites = len(atomlist)

    with open( path + filename , 'w+') as newpdb:

        for i in range(sites):

            c1 = 'ATOM'
            c2 = i+1 #str(i+1).rjust(7)
            c3 = atomlist[i,0] #str(atomlist[i,0])
            c4 = ''
            c5 = atomlist[i,1] #str(atomlist[i,1])
            c6 = atomlist[i,2]
            c7 = i+1 #str(i+1)
            c8 = ''
            c9 = "%.3f" % np.around(ftraj_frame[3*i], 3)
            c10 = "%.3f" % np.around(ftraj_frame[3*i+1],3)
            c11 = "%.3f" % np.around(ftraj_frame[3*i+2],3)
            c12= 1.00
            c13= 0.00
            c14= atomlist[i,0] #str(atomlist[i,0])
            c15= ''

            #line=c1+c2+c3+c4+c5+c6+c7+c8+c9+c10+c11+c12+c13
            line= "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8s}{:8s}{:8s}{:6.2f}{:6.2f}          {:>2s}{:2s}".format(c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15) +'\n'

            newpdb.write(line)

    return print('Done!')


def fric_calc_Pmodel(TOP, protname, N, nfrs, avblsq, T, intv = 2.71828, viscosity = 3.1e-4, fd20 = 0.0, path_to_resarea = './', cgmodel = 'P'):
    
    # !!!! TOP must contain only coarse grained sites. !!!!
    
    import numpy as np
    import sys
    import os
    import subprocess
    
    pi = np.pi
    
    print(protname, cgmodel, N, nfrs)
    
    if cgmodel != 'P':
        raise ValueError('''Available coarse grain models are: P.''')
        
    # All volumes given are pulled from the papers below.   
    #Deoxynucleotide volumes from "Standard atomic volumes in double-stranded DNA andpacking in protein–DNA interfaces" by Nadassy et al. NAR, 2001.
    #Ribonucleotide volumes from "Calculation of Standard Atomic Volumes for RNA andComparison with Proteins: RNA is Packed More Tightly" by Voss and Gersteain, JMB, 2005.
    mrad_dict = { "DA5" : np.cbrt((3*310.9) / (4*pi)) ,
                  "DA3" : np.cbrt((3*310.9) / (4*pi)) ,
                  "DA"  : np.cbrt((3*310.9) / (4*pi)) ,
                  "DC5" : np.cbrt((3*288.0) / (4*pi)) ,
                  "DC3" : np.cbrt((3*288.0) / (4*pi)) ,
                  "DC"  : np.cbrt((3*288.0) / (4*pi)) ,
                  "DG5" : np.cbrt((3*318.6) / (4*pi)) ,
                  "DG3" : np.cbrt((3*318.6) / (4*pi)) ,
                  "DG"  : np.cbrt((3*318.6) / (4*pi)) ,
                  "DT5" : np.cbrt((3*307.4) / (4*pi)) ,
                  "DT3" : np.cbrt((3*307.4) / (4*pi)) ,
                  "DT"  : np.cbrt((3*307.4) / (4*pi)) ,
                  "A5"  : np.cbrt((3*315.3) / (4*pi)) ,
                  "A3"  : np.cbrt((3*315.3) / (4*pi)) ,
                  "A"   : np.cbrt((3*315.3) / (4*pi)) ,
                  "C5"  : np.cbrt((3*291.1) / (4*pi)) ,
                  "C3"  : np.cbrt((3*291.1) / (4*pi)) ,
                  "C"   : np.cbrt((3*291.1) / (4*pi)) ,
                  "G5"  : np.cbrt((3*322.0) / (4*pi)) ,
                  "G3"  : np.cbrt((3*322.0) / (4*pi)) ,
                  "G"   : np.cbrt((3*322.0) / (4*pi)) ,
                  "U5"  : np.cbrt((3*286.9) / (4*pi)) ,
                  "U3"  : np.cbrt((3*286.9) / (4*pi)) ,
                  "U"   : np.cbrt((3*286.9) / (4*pi))  }
                 
                 
    #Calculate the Miller radius per bead
    mradlist = []
    #P_mradlist = []
    #S_mradlist = []
    with open(TOP) as f:
        for line in f:
            if line[0:4]=='ATOM':
                dummy=line.split()
                if dummy[2] == "P":
                    mradlist.append( mrad_dict[dummy[3]]) 
                else:
                    raise Exception('The TOP file contains atoms which are not phosphates.')
            else:
                pass
    mradlist = np.array(mradlist)
            
    #file should look like: "./SP_resarea.xvg"       
    if os.path.exists(path_to_resarea + cgmodel + "_resarea.xvg"):
        pass
    else:
        raise FileNotFoundError('''I can't find the resarea.xvg file containing the solvent-exposed surface area of each residue.
                                Please either run the process.sh file, if you have not already done so, and move the resarea.xvg 
                                file into the current working directory.''')
    
    #gmx_mpi sasa -f after_rot_sim0.xtc -s SP_dsAXA.pdb -n SP_index.ndx -or SP_resarea.xvg -dt 1000 -surface 1 -output 2 3
    
    #determine solvent exposed surface areas of each group from the .xvg file
    # this assumes lines in .xvg have the following format:
    # res  resarea  std      Sarea    std      Parea    std
    # 1    2.211    0.179    1.118    0.102    0.000    0.000
    # 2    1.942    0.195    0.603    0.085    0.842    0.076 ...
    P_resarea = []
    with open('./' + 'P' + "_resarea.xvg") as f:
        for line in f:
            if line[0] == '#' or line[0] == '@':
                pass
            else:
                dummy=line.split()
                if float(dummy[3])==0.0:
                    pass
                else:
                    P_resarea.append(float(dummy[1]))
                       
    P_resarea = np.array(P_resarea)
    
    #create list of CG-site solvent exposed radii
    P_rad = ((P_resarea/(4*pi))**0.5)*10
    
    #np.savetxt('avresrad',np.array(rad),fmt="%f")
    fratio = (P_rad.sum()/N)/10
    print('fratio: ', fratio)
    
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
    
    rv = mradlist    #np.array(mradlist)
    rw = P_rad    #np.array(rad)
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
