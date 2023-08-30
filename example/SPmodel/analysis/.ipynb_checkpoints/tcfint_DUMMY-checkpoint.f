	program inputreader
	integer nfrs,n
	open(unit=5,file="protname.txt",status='old')
	read(5,*)
	read(5,*)n
	read(5,*)nfrs
	close(5)
        nfrs = nfrs - 1
	write(*,*)n,nfrs
!	nfrs=500 !testing!
	call umatrix(n,nfrs)
	End program inputreader

	subroutine umatrix(n,nfrs)
	integer i,j,k,l,m,imin,jmin,a,nca,n,nfrs
	DOUBLE PRECISION :: lavmsq
        DOUBLE PRECISION :: rx1, ry1, rz1, rx2, ry2, rz2
	DOUBLE PRECISION, dimension(nfrs) :: lx,ly,lz
	DOUBLE PRECISION, dimension(0:nfrs) :: dipcorr
	character(32)protname
	character(16)aa,ii
	open(unit=5,file="protname.txt",status='old')
	read(5,'(A)')protname
	close(5)
	protname=adjustl(protname)
	lx = 0.0
	ly = 0.0
	lz = 0.0
	lavmsq = 0.0

	!read from trajectory
	open(unit=13,file='particle_spam',status='old')
        open(unit=14,file='particle_eggs',status='old')
        
        !Use the residue trajectories to calculate the bond
        !vector at each time frame
	do k=1,nfrs
	  rx1 = 0.0
	  ry1 = 0.0
	  rz1 = 0.0
          rx2 = 0.0
          ry2 = 0.0
          rz2 = 0.0
	  if(mod(k,5000).eq.0)write(*,*)"frame",k
	  read(13,*)rx1, ry1, rz1
          read(14,*)rx2, ry2, rz2

          !Calculate bond vector and mean-square bond length
	  lx(k) = rx2 - rx1
	  ly(k) = ry2 - ry1
	  lz(k) = rz2 - rz1
	  lavmsq = lavmsq + lx(k)**2 + ly(k)**2 + lz(k)**2


	end do
	close(13)
	close(14)

	lavmsq=lavmsq/real(nfrs)

	dipcorr=0.0
	!calculate tcf of bond autocorrelation
	l=dummy
	write(ii,*)l
	ii=adjustl(ii)
	open(unit=200,file='m1CAsimint_'//trim(ii))
	!open(unit=100,file='avm1CAsimint')
	do j=0,nfrs/4
	!bond loop
	if(mod(j,100).eq.0)write(*,*)"tcf frame:",j
	do k=1,nfrs-j
!and now the tcf
	dipcorr(j)=dipcorr(j)+((lx(k)*lx(j+k)
     &+ly(k)*ly(j+k)+lz(k)*lz(j+k)))
     &/(lavmsq)
	end do
	dipcorr(j)=dipcorr(j)/real(nfrs-j)
	write(200,*)j*.2, dipcorr(j)
	!avdipcorr(j)=avdipcorr(j)/real(n-1)
	!write(100,*)j*.2,avdipcorr(j)
	end do !end j time loop
	!close(100)
	close(200)

	end subroutine

