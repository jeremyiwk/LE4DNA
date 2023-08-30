!******************************************************************************************
! This program computes the theoretical bond autocorrelation function from the LE4PD theory
! as described in the paper:
!
! Copperman, J., and M. G Guenza. "Coarse-Grained Langevin Equation for Protein Dynamics: 
! Global Anisotropy and a Mode Approach to Local !Complexity." The Journal of Physical 
! Chemistry. B 119.29 (1900): 9195-211. Web.
! 
! The rescaling of the barrier crossing timescales can be found in this paper:
!
! Beyerle, Eric R, and Marina G Guenza. "Kinetics Analysis of Ubiquitin Local 
! Fluctuations with Markov State Modeling of the LE4PD Normal Modes." 
! The Journal of Chemical Physics 151.16 (2019): 164119. Web.
!******************************************************************************************

        program inputreader
        IMPLICIT NONE
        INTEGER :: n,nfrs,nf
        open(unit=114,file='protname.txt',status='old')
        read (114,*)
        read(114,*)n
        read(114,*)nfrs
        close(114)
        nf=(nfrs-1)/4  !use same number of frames as in sim. bond tcf
        call m1(n,nf)
        end program inputreader



        subroutine m1(n,nf)
        
        IMPLICIT NONE
        INTEGER :: i,j,k,a
        INTEGER :: nm1,nmol
        INTEGER :: n,nfrs,nf
        REAL :: tstep
        DOUBLE PRECISION :: sigma,T,Rb,lsq
        DOUBLE PRECISION :: decay,lnorm
        DOUBLE PRECISION :: eiglam1,eiglam2,eiglam3
        DOUBLE PRECISION, dimension(n,n) :: qinvm(n,n),qm(n,n)
        DOUBLE PRECISION, dimension(n,n) :: cc(n,n)
        DOUBLE PRECISION, dimension(n) :: eigmu(n),eiglam(n)
        DOUBLE PRECISION, dimension(n) :: fricorr(n),tau(n)
        DOUBLE PRECISION, dimension(nf) :: tc(nf),tt(nf)
        DOUBLE PRECISION, dimension(nf) :: m1sum(nf)
        character(32)protname,ii

        open(unit=51, file='nmol.dat', status='old')
        open(unit=52, file='sigma.dat', status='old')
        open(unit=53, file='barriers_kcal.dat', status='old')
        open(unit=54, file='lambda_eig', status='old')
        open(unit=55, file='mu_eig', status='old')
        open(unit=56, file='Qmatrix', status='old')
        open(unit=57, file='QINVmatrix', status='old')
        open(unit=58, file='avblsq.dat',status='old')
        
        read(51,*)nmol
        read(52,*)sigma
        read(58,*)lsq
        
        nm1=n-nmol     ! number of internal modes
        T=300.0        ! Kelvin
        Rb=0.00198     ! boltzmanns constant in kcal/mol*K
        tstep=0.2      ! time-step in simulation in ps
        
        do i=1,nm1
            read(53,*)fricorr(i)
            fricorr(i)=exp(fricorr(i)/(Rb*T))
        end do
        
        do i=1,nm1
            read(54,*)eiglam(i)
        end do
        
        do i=1,nm1
            read(55,*)eigmu(i)
        end do
        
        do i=1,nm1
            do j=1,nm1
                read(56,*)qm(i,j)
            end do
        end do
        
        do i=1,nm1
            do j=1,nm1
                read(57,*)qinvm(i,j)
            end do
        end do
        
        close(51)
        close(52)
        close(53)
        close(54)
        close(55)
        close(56)
        close(57)
        close(58)
        
        do i=1,nm1
            do j=1,nm1
                cc(i,j)=lsq*(qm(i,j)**2)/eigmu(j)
            end do
        end do
        
        do a=1,nm1
            tau(a)=(fricorr(a)/(sigma*eiglam(a)))
        end do
        
        ! The lines below correct the eigenvalues for
        ! diffusion in the lab-frame
        ! eiglam1=.5*(eiglam(1)+eiglam(2))
        ! eiglam2=.5*(eiglam(1)+eiglam(3))
        ! eiglam3=.5*(eiglam(2)+eiglam(3))
        ! eiglam(1)=eiglam1
        ! eiglam(2)=eiglam2
        ! eiglam(3)=eiglam3
        
        do i=1,nm1  !loop over bonds to write each m1_i
        
            write(ii,*)i   !set the string "ii" equal to i
            ii=adjustl(ii) !remove leading spaces in "ii"
            open(unit=60,file='m1_'//trim(ii),status='unknown')
            write(60,*)dble(0.0),dble(1.0)
            do k=1,nf   !loop over frames  
            
                m1sum(k)=0.0
                tt(k)= tstep*k
                lnorm=0.0

                do a=1,nm1   !loop over modes 1-(N-2)
                
                    if (a.le.3) then                  !
                        decay=1.0                     !
                    else if (a.ge.4) then             !
                        decay=exp(-tt(k)/tau(a))      !
                    end if                            ! the first three modes do not decay
                    
                    lnorm=lnorm+cc(i,a)
                    m1sum(k)=m1sum(k)+cc(i,a)*decay
                
                end do

                write(60,*)tt(k),m1sum(k)/lnorm
                
            end do
                close(60)
        end do
        
        return
    
        end

