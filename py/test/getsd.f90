!pgf90 getsd.f90 -llapack

program analyze
  
  implicit none
  
  integer, parameter :: natoms = 105*3, nmodes = natoms
  real, parameter :: scfactor=1.00d0                     !scaling factor for frequencies
  
  ! unit conversion factors
  real(8), parameter :: &
       planck = 6.626068d-34,                                                         & ! Planck's constant in J*s
       avog = 6.0221415d23,                                                           & ! Avogadro's number
       light = 2.99792458d10,                                                         & ! Speed of light (SI) in cm/s 
       pi = 4.0d0*atan(1.0d0)   
  ! derived conversion factors
  real(8), parameter :: &! **************** ENERGY *************** !
       joules2wvnbr = 1.d0/planck/light,                                              & ! Joules (SI) to wavenumber
       au2joules = 4.35974434d-18,                                                    & ! Hartrees to Joules
       wvnbr2joules = 1.d0/joules2wvnbr,                                              & ! wavenumber to Joules
       joules2au = 1.d0/au2joules,                                                    & ! Joules to Hartrees
       ! **************** LENGTH *************** !
       au2ang = 0.52917725d0,                                                         & ! atomic units of length to Angstrom
       ang2au = 1.0d0/au2ang,                                                         & ! Angstrom to atomic units of length 
       ang2m = 1.d-10,                                                                & ! Angstrom to meter
       ! ***************** TIME **************** !
       amu2kg = 1.66053892d-27,                                                       & ! atomic mass unit to kg
       ! ************** FREQUENCY ************** !
       wvnbr2Hz = 3.0d10,                                                             & ! wavenumbers to Hertz
       Hz2wvnbr = 1.0d0/wvnbr2Hz,                                                     & ! Hertz to wavenumbers
       Hz2aufreq = 2.418884d-17 
  
  !for data workup
  integer :: i, j, l, npoint, itmp, imass
  real(8) :: disp=0.d0, reorg=0.d0, hrfactor=0.d0, sbc=0.d0, totallambda=0.d0, lambdanorm=0.d0, sum=0.d0
  real(8) :: fc00=0.d0, fc01=0.d0, fc02=0.d0, fc11=0.d0
  real(8) :: shrfactor=0.d0, sreorg=0.d0, sfc00=0.d0, sfc01=0.d0, sfc02=0.d0, sfc11=0.d0
  real(8), dimension(natoms) :: eforce=0.d0, gforce=0.d0, mass=0.d0, ntot=0.d0
  real(8), dimension(nmodes) :: presd = 0.0d0, freq=0.d0, tempfreq=0.d0, redmass=0.d0, tempredmass=0.d0, k=0.d0
  real(8), dimension(nmodes,natoms) :: nmvector=0.d0, tempnmvector=0.d0
  real(8), allocatable :: gsd(:), gfreq(:)
  character(100) :: infile, outfile, ifile
  
  !for gaussian broadening J(w)
  integer(8) :: nfreq, nbin
  real(8) :: sigma=20.d0, dw=0.53195053d0, norm=0.d0, dumfreq=0.d0 ! in cm-1
  
  !for dgetri and dgetrf
  integer(8), dimension(natoms) :: ipiv=0
  real(8), dimension(nmodes,natoms) :: inv=0.d0
  
  !stupid stuff
  integer(8) :: idum
  real(8) :: rdum
  character(100) :: cdum
  
  !for rotational modes
  character(2), dimension(natoms/3) :: identity
  real(8), dimension(natoms/3,3) :: coord=0.d0
  real(8), dimension(3) :: center=0.d0
  real(8), dimension(natoms/3) :: rotmass=0.d0
  real(8) :: totalmass=0.d0
  real(8), dimension(3,3) :: inertia=0.d0
  real(8), dimension(3) :: w=0.d0
  real(8), dimension(natoms/3) :: px=0.d0, py=0.d0, pz=0.d0

  !for adding orthogonal rotational and translational vectors
  integer :: ivec
  real(8) :: dot_nm, new_vec(natoms)
  ! for testing inversion
  real(8), dimension(2,2) :: a
  

  
!  !read forces in amu, hartrees/bohr
!  itmp=0 ; open(unit=10, file='eforces.dat')
!  do i=1, natoms/3
!     read(10,*) idum, imass, eforce(itmp+1), eforce(itmp+2), eforce(itmp+3)
!     if (imass.eq.15) mass(itmp+1)=2.d0*dble(imass)
!     if (imass.eq.1) mass(itmp+1)=dble(imass)
!     if (imass.eq.6) mass(itmp+1)=2.d0*dble(imass)
!     if (imass.eq.7) mass(itmp+1)=2.d0*dble(imass)
!     if (imass.eq.8) mass(itmp+1)=2.d0*dble(imass)
!     mass(itmp+2)=mass(itmp+1) ; mass(itmp+3)=mass(itmp+1)
!     itmp=i*3
!  end do
!  close(10)
!
!
!  itmp=0 ; open(unit=10, file='gforces.dat')
!  do i=1, natoms/3
!     read(10,*) idum, idum, gforce(itmp+1), gforce(itmp+2), gforce(itmp+3)
!     itmp=i*3
!  end do
!  close(10)
 
 
!  mass = mass * amu2kg
  eforce = eforce * au2joules / (au2ang*ang2m) / dsqrt(mass)
  gforce = gforce * au2joules / (au2ang*ang2m) / dsqrt(mass)
 
!  open(unit=10, file='modes.dat')
!  do i=1 ,nmodes-6
!     read(10,*) cdum, idum, cdum, freq(i), cdum, redmass(i)
!     if ( idum .ne. i ) exit
!     do j=1,natoms
!        read(10,*) idum, idum, idum, nmvector(i,j)
!        nmvector(i,j)=nmvector(i,j)*dsqrt(mass(j))
!        norm=norm+nmvector(i,j)*nmvector(i,j)
!     end do
!     nmvector(i,:)=nmvector(i,:)/dsqrt(norm)
!     norm=0.d0
!  end do
!  close(10)


!  tempfreq=freq ; tempredmass=redmass ; tempnmvector=nmvector
!  freq=0.d0 ; redmass=0.d0 ; nmvector=0.d0 !freq
  
!  do i=1, nmodes-6
!     freq(i+6)=tempfreq(i)
!     redmass(i+6)=tempredmass(i)
!     nmvector(i+6,:)=tempnmvector(i,:)
!  end do
  
!  do i=1, nmodes
!     if (redmass(i).eq.0.d0) then
!        redmass(i) = 1.d0
!     end if
!  end do
  
  !add in translational eigenvectors 
  itmp=1
  do i=1, natoms/3
     nmvector(1,itmp) = dsqrt(mass(i)) 
     nmvector(2,itmp+1) = dsqrt(mass(i))
     nmvector(3,itmp+2) = dsqrt(mass(i))
     itmp=itmp+3     
  end do

  !  read in coordinates of atoms for rotational eigenvector
  open(unit=60, file='coords.dat')
  do i=1,natoms/3
     read(60,*), identity(i), (coord(i,j), j=1,3)
  end do
  close(60)
  
  ! calculate center of mass
  do i=1, natoms/3
     if (identity(i) .eq. 'P') rotmass(i)=2.d0*dble(15)
     if (identity(i) .eq. 'C ') rotmass(i)=2.d0*dble(6)
     if (identity(i) .eq. 'N ') rotmass(i)=2.d0*dble(7)
     if (identity(i) .eq. 'O ') rotmass(i)=2.d0*dble(8)
     if (identity(i) .eq. 'H ') rotmass(i)=dble(1)
  end do
  
  
  do i=1, natoms/3
     coord(i,:)=coord(i,:)*rotmass(i)
     totalmass=totalmass+rotmass(i)
     center=center+coord(i,:)
  end do

  center=center/totalmass
  
  do i=1,natoms/3
     coord(i,:)=coord(i,:)-center
  end do

  ! add translational/rotational vectors that are orthogonal
!  ivec = 7
!333 continue
!  if (ivec .gt. 1) then
!    new_vec = 1.0
!    do i=ivec,nmodes,1
!      new_vec = new_vec-dot_product(new_vec,nmvector(i,:))*nmvector(i,:)
!    end do
!    ivec = ivec-1
!    nmvector(ivec,:) = new_vec/dsqrt(dot_product(new_vec,new_vec))
!go to 333 
! end if

!do i=1,nmodes,1
!  do j=0,natoms,1
!    write(*,*) nmvector(i,3*j+1), nmvector(i,3*j+2), nmvector(i,3*j+3)
!  end do
!end do

  ! calculate  moment of inertia tensor
  do i=1,natoms/3
     inertia(1,1) = inertia(1,1) + rotmass(i)*(coord(i,2)*coord(i,2) + coord(i,3)*(coord(i,3)))
     inertia(2,2) = inertia(2,2) + rotmass(i)*(coord(i,1)*coord(i,1) + coord(i,3)*(coord(i,3)))
     inertia(3,3) = inertia(3,3) + rotmass(i)*(coord(i,1)*coord(i,1) + coord(i,2)*(coord(i,2)))
     inertia(1,2) = inertia(1,2) - rotmass(i)*coord(i,1)*coord(i,2)
     inertia(1,3) = inertia(1,3) - rotmass(i)*coord(i,1)*coord(i,3)
     inertia(2,3) = inertia(2,3) - rotmass(i)*coord(i,2)*coord(i,3)
  end do
  inertia(2,1) = inertia(1,2) ; inertia(3,1) = inertia(1,3) ; inertia(3,2) = inertia(2,3)

  ! now diagonalized the moment of inertia tensor
  call get_zheev(inertia,3,w)
  
  ! build the rotational eigenvectors
  itmp=0
  do i=1,natoms/3
     px(i) = ( inertia(1,1)*coord(i,1) + inertia(2,1)*coord(i,2) + inertia(3,1)*coord(i,3) )
     py(i) = ( inertia(1,2)*coord(i,1) + inertia(2,2)*coord(i,2) + inertia(3,2)*coord(i,3) )
     pz(i) = ( inertia(1,3)*coord(i,1) + inertia(2,3)*coord(i,2) + inertia(3,3)*coord(i,3) )
     nmvector(4,itmp+1) = (py(i)*inertia(3,1) - pz(i)*inertia(2,1))/dsqrt(rotmass(i))
     nmvector(4,itmp+2) = (py(i)*inertia(3,2) - pz(i)*inertia(2,2))/dsqrt(rotmass(i))
     nmvector(4,itmp+3) = (py(i)*inertia(3,3) - pz(i)*inertia(2,3))/dsqrt(rotmass(i))
     nmvector(5,itmp+1) = (pz(i)*inertia(1,1) - px(i)*inertia(3,1))/dsqrt(rotmass(i))
     nmvector(5,itmp+2) = (pz(i)*inertia(1,2) - px(i)*inertia(3,2))/dsqrt(rotmass(i))
     nmvector(5,itmp+3) = (pz(i)*inertia(1,3) - px(i)*inertia(3,3))/dsqrt(rotmass(i))
     nmvector(6,itmp+1) = (px(i)*inertia(2,1) - py(i)*inertia(1,1))/dsqrt(rotmass(i))
     nmvector(6,itmp+2) = (px(i)*inertia(2,2) - py(i)*inertia(1,2))/dsqrt(rotmass(i))
     nmvector(6,itmp+3) = (px(i)*inertia(2,3) - py(i)*inertia(1,3))/dsqrt(rotmass(i))
     itmp=i*3
  end do
  
  !  normalize translational + rotational eigenvectors
  do i=1,6
     norm=0.d0
     do j=1,natoms
        norm = norm + nmvector(i,j)*nmvector(i,j)
     end do
     nmvector(i,:)=nmvector(i,:)/dsqrt(norm)
  end do

  !test  
  open(50,file='test.dat')
  do i=1,nmodes
     write(50,*), 'freq=',freq(i),'redmass=',redmass(i)
     do j=1,natoms
        write(50,*), 'nmvec=',nmvector(i,j)
     end do
  end do
  close(50)
  

     
     freq = freq * wvnbr2Hz * 2.d0 * pi * scfactor 
     redmass = redmass * amu2kg

          
     inv=nmvector

     call get_dgetrf(natoms,nmodes,inv,ipiv)
     call get_dgetri(natoms,nmodes,inv,ipiv)
     
     open(unit=30,file='results.dat')
     open(unit=40,file='fc.dat')
     open(unit=50,file='hr0d1_fc.dat')

     do i=7, nmodes  ! 7 to nmodes
        do j=1,natoms
           k(i) = k(i) + eforce(j)*inv(j,i) - gforce(j)*inv(j,i)
        end do


        disp=k(i)/freq(i)/freq(i)
        hrfactor=0.5d0*freq(i)*disp*disp/(planck/2.d0/pi)
        reorg=0.5d0*freq(i)*freq(i)*disp*disp*joules2wvnbr
        sbc = freq(i)*freq(i)*disp
        fc00=exp(-0.5d0*hrfactor)
        fc01=sqrt(hrfactor)*fc00
        fc02=hrfactor*fc00
        fc11=(1.d0-hrfactor)*fc00
        shrfactor=0.1d0+hrfactor
        !sreorg=0.5d0*freq(i)*freq(i)*disp*disp*joules2wvnbr+(shrfactor*freq(i)*(planck/2.d0/pi))
        sreorg=reorg
        sfc00=exp(-0.5d0*shrfactor)
        sfc01=sqrt(shrfactor)*fc00
        sfc02=shrfactor*sfc00
        sfc11=(1.d0-(shrfactor))*fc00
        totallambda=totallambda + reorg

        write(30,'(F30.20,F30.20,F30.20,F30.20,F30.20)') freq(i)*Hz2wvnbr/2.d0/pi, reorg, hrfactor, sbc, sqrt(sbc**2)
        write(40,'(F30.20,F30.20,F30.20,F30.20,F30.20,F30.20,F30.20)') freq(i)*Hz2wvnbr/2.d0/pi, reorg, hrfactor, fc00**2, fc01**2, fc02**2, fc11**2
        write(50,'(F30.20,F30.20,F30.20,F30.20,F30.20,F30.20,F30.20)') freq(i)*Hz2wvnbr/2.d0/pi, sreorg, shrfactor, sfc00**2, sfc01**2, sfc02**2, sfc11**2


        presd(i)=0.5d0*pi*k(i)*k(i)/freq(i)

     end do
     close(unit=30)
     close(unit=40)
     close(unit=50)
     
     
     
     !gaussianbroadening
     freq=freq/wvnbr2Hz/2.d0/pi
     nfreq=nint(3000.0/dble(dw))
     allocate(gsd(nfreq*nfreq),gfreq(nfreq*nfreq))
     
     do i=1, nmodes
        do j = 1, nfreq
           dumfreq = freq(i) + dble(j - dble(nfreq)/2.d0)*dw
           npoint = nint(dumfreq/dw)
           gsd(npoint) = gsd(npoint) + presd(i)*exp( -(freq(i) - dumfreq)**2/sigma/sigma)
        end do
     enddo
     
     !normalizing SD to give correct total reorganization energy
     do i=1,nfreq
        gfreq(i) = dw * dble(i)
        lambdanorm=lambdanorm + gsd(i)/gfreq(i)*dw
     end do
     lambdanorm=lambdanorm/pi ; lambdanorm=totallambda/lambdanorm
     
     open (unit=60, file='sd-gaussbroad.dat')
     do j=1, nfreq
        write(60,*) gfreq(j) , gsd(j) * lambdanorm
     enddo
     close(60)
     
   end program analyze
   
   subroutine get_dgetrf(natoms,nmodes,inv,ipiv)
    implicit none
     integer(8) :: nmodes,natoms,info   
     integer(8) :: ipiv(natoms)    
     real(8) :: inv(natoms,natoms)          
     
     call dgetrf(nmodes,natoms,inv,natoms,ipiv,info)
     !write(*,*) info
     !write(*,*) ipiv 
   end subroutine get_dgetrf
   subroutine get_dgetri(natoms,nmodes,inv,ipiv)
     implicit none
     integer :: natoms,nmodes,info
     integer :: ipiv(natoms)
     real(8) :: inv(nmodes,natoms)
     real(8) :: work(100*nmodes)
     integer :: lwork
     
     !ipiv = 0
     lwork = 100*nmodes
     
     !write(*,*) ipiv
     call dgetri(natoms,inv,nmodes,ipiv,work,lwork,info)

   end subroutine get_dgetri
   
   
   subroutine get_zheev(inertia,n,w)
     implicit none
     real(8),dimension(n,n) :: inertia
     complex*16 :: a(n,n)
     integer :: n
     real*8 :: w(n)
     character(1) :: jobz='V'
     character(1) :: uplo='U'
     integer::info
     complex*16,dimension(:),allocatable :: work
     real*8,dimension(:),allocatable :: rwork
     
     a = (0.d0,0.d0)
     a = inertia
     
     if(allocated(work) == .true.) deallocate(work)
     if(allocated(rwork) == .true.) deallocate(rwork)
     allocate(work(2*n),rwork(3*n-2))
     
     call zheev(jobz,uplo,n,a,n,w,work,2*n,rwork,info)
     
     !print*, inertia
     !!print*, a
     inertia = real(a)
     !print*, inertia
     
     if(allocated(work) == .true.) deallocate(work)
     if(allocated(rwork) == .true.) deallocate(rwork)
     
   end subroutine get_zheev
   
