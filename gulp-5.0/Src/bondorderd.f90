  subroutine bondorderd(ebondorder,xkv,ykv,zkv,lgrad1,lgrad2,ldophonon)
!
!  Calculates the energy and derivatives for the Bond Order potentials.
!
!  Distributed memory version. 
!
!  On entry : 
!
!  lgrad1          = if .true. calculate the first derivatives
!  lgrad2          = if .true. calculate the second derivatives
!  ldophonon       = if .true. then this is a phonon run and so
!                    second derivatives should be phased
!  xkv             = if ldophonon this is the x kvector component
!  ykv             = if ldophonon this is the y kvector component
!  zkv             = if ldophonon this is the z kvector component
!
!  On exit :
!
!  ebondorder      = the value of the energy contribution
!
!   1/17 Created from bondorder
!   4/17 Debug printing now only on ioproc
!
!  Conditions of use:
!
!  GULP is available free of charge to academic institutions
!  and non-commerical establishments only. Copies should be
!  obtained from the author only and should not be distributed
!  in any form by the user to a third party without the express
!  permission of the author. This notice applies to all parts
!  of the program, except any library routines which are
!  distributed with the code for completeness. All rights for
!  such routines remain with the original distributor.
!
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, January 2017
!
  use datatypes
  use bondorderdata
  use configurations, only : nregionno, nregions, lsliceatom, nregiontype, QMMMmode
  use control,        only : keyword, lseok
  use current
  use energies,       only : eattach, esregion12, esregion2
  use iochannels
  use neighbours
  use optimisation,   only : lfreeze, lopf
  use parallel
  use spatialbo,      only : lspatialok => lspatialBOok
  use times
  implicit none
!
!  Passed variables
!
  real(dp),    intent(out)                         :: ebondorder
  real(dp),    intent(in)                          :: xkv
  real(dp),    intent(in)                          :: ykv
  real(dp),    intent(in)                          :: zkv
  logical,     intent(in)                          :: lgrad1
  logical,     intent(in)                          :: lgrad2
  logical,     intent(in)                          :: ldophonon
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ii
  integer(i4)                                      :: ind
  integer(i4), dimension(:,:,:), allocatable, save :: ineigh         ! Cell indices for vector
  integer(i4)                                      :: itmp
  integer(i4)                                      :: j
  integer(i4)                                      :: k
  integer(i4)                                      :: kk
  integer(i4)                                      :: maxneigh2
  integer(i4)                                      :: maxneigh22
  integer(i4)                                      :: mA
  integer(i4)                                      :: mR
  integer(i4)                                      :: n
  integer(i4)                                      :: nati
  integer(i4)                                      :: natj
  integer(i4)                                      :: natk
  integer(i4)                                      :: nboij
  integer(i4)                                      :: nboAij
  integer(i4)                                      :: nboAik
  integer(i4)                                      :: nboRij
  integer(i4)                                      :: nboRik
  integer(i4)                                      :: nboZi
  integer(i4), dimension(:,:),   allocatable, save :: nbopotptr
  integer(i4)                                      :: nboatom
  integer(i4), dimension(:),     allocatable, save :: nboatomRptr
  integer(i4)                                      :: ni
  integer(i4)                                      :: nj
  integer(i4)                                      :: njk
  integer(i4)                                      :: nmin
  integer(i4)                                      :: nn
  integer(i4)                                      :: nn1
  integer(i4)                                      :: nn2
  integer(i4)                                      :: npki   
  integer(i4)                                      :: nptr
  integer(i4), dimension(:,:),   allocatable, save :: neighno
  integer(i4), dimension(:),     allocatable, save :: nfreeatom
  integer(i4), dimension(:),     allocatable, save :: nneigh
  integer(i4)                                      :: nneighi2
  integer(i4)                                      :: nneighj2
  integer(i4)                                      :: nneighi22
  integer(i4)                                      :: nneighj22
  integer(i4)                                      :: nregioni
  integer(i4)                                      :: nregionj
  integer(i4)                                      :: nregiontypi
  integer(i4)                                      :: nregiontypj
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: ntypj
  integer(i4)                                      :: ntypk
  integer(i4)                                      :: status
  logical,     dimension(:),     allocatable, save :: latomdone
  logical                                          :: lattach
  logical                                          :: lfound
  logical                                          :: lilocal
  logical                                          :: lmaxneighok
  logical                                          :: lneedBOcnderiv
  logical                                          :: lQMMMok
  logical                                          :: lreg2one
  logical                                          :: lreg2pair
  logical                                          :: lslicei
  logical                                          :: lslicej
  logical,     dimension(:),     allocatable, save :: lopanyneigh
  real(dp)                                         :: bijA
  real(dp)                                         :: bijR
  real(dp)                                         :: bijsumA
  real(dp)                                         :: bijsumA1
  real(dp)                                         :: bijsumAn1
  real(dp)                                         :: bijsumAn2
  real(dp)                                         :: bijsumR
  real(dp)                                         :: bijsumR1
  real(dp)                                         :: bijsumRn1
  real(dp)                                         :: bijsumRn2
  real(dp)                                         :: BOncoAij
  real(dp)                                         :: BOncoRij
  real(dp)                                         :: btotA
  real(dp)                                         :: btotR
  real(dp)                                         :: g_cpu_time
  real(dp),    dimension(:),     allocatable, save :: d1i
  real(dp),    dimension(:),     allocatable, save :: d2i
  real(dp),    dimension(:,:),   allocatable, save :: d1ZiCN
  real(dp),    dimension(:),     allocatable, save :: d1BtotiA
  real(dp),    dimension(:),     allocatable, save :: d1BtotiR
  real(dp),    dimension(:),     allocatable, save :: d2BtotiA
  real(dp),    dimension(:),     allocatable, save :: d2BtotiR
  real(dp)                                         :: dedZ
  real(dp)                                         :: d2edZ2
  real(dp)                                         :: dexpijkdr
  real(dp)                                         :: d2expijkdr2
  real(dp)                                         :: dfdr
  real(dp)                                         :: dfikdr
  real(dp)                                         :: dfzdz
  real(dp)                                         :: d2fdr2
  real(dp)                                         :: d2fikdr2
  real(dp)                                         :: d2fzdz2
  real(dp)                                         :: d3fdr3
  real(dp)                                         :: d3fikdr3
  real(dp)                                         :: d3fzdz3
  real(dp)                                         :: dGijkdr(3)
  real(dp)                                         :: d2Gijkdr2(6)
  real(dp)                                         :: d3Gijkdr3(10)
  real(dp)                                         :: dZ
  real(dp)                                         :: deltaZ
  real(dp)                                         :: eij
  real(dp)                                         :: expijk
  real(dp)                                         :: f
  real(dp)                                         :: fik
  real(dp)                                         :: fz
  real(dp)                                         :: Gijk
  real(dp)                                         :: RmA
  real(dp)                                         :: RmR
  real(dp)                                         :: rbijsumA1
  real(dp)                                         :: rbijsumR1
  real(dp),    dimension(:),     allocatable, save :: rBOcutmax
  real(dp)                                         :: rij
  real(dp)                                         :: rik
  real(dp)                                         :: rlambda
  real(dp)                                         :: drlambdadrij
  real(dp)                                         :: drlambdadrik
  real(dp)                                         :: rrij
  real(dp)                                         :: rrik
  real(dp)                                         :: rtmp
  real(dp)                                         :: t1
  real(dp)                                         :: t2
  real(dp),    dimension(:,:),   allocatable, save :: rneigh
  real(dp),    dimension(:,:),   allocatable, save :: xneigh
  real(dp),    dimension(:,:),   allocatable, save :: yneigh
  real(dp),    dimension(:,:),   allocatable, save :: zneigh
  real(dp)                                         :: Va
  real(dp)                                         :: Vr
  real(dp)                                         :: dVadr
  real(dp)                                         :: dVrdr
  real(dp)                                         :: d2Vadr2
  real(dp)                                         :: d2Vrdr2
  real(dp)                                         :: xdiff
  real(dp)                                         :: ydiff
  real(dp)                                         :: zdiff
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xki
  real(dp)                                         :: yki
  real(dp)                                         :: zki
  real(dp)                                         :: zAi3
  real(dp)                                         :: zRi3
  real(dp)                                         :: zsign
  real(dp),    dimension(:),     allocatable, save :: Zcn
!
  t1 = g_cpu_time()
  allocate(nboatomRptr(numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','nboatomRptr')
! 
!  Set up a dummy pointer for derivative array calls
! 
  nboatom = 0
  do i = 1,numat                                 
    nboatom = nboatom + 1                  
    nboatomRptr(i) = nboatom               
  enddo
!
!  Allocate memory that does not depend on maxneigh
!
  allocate(latomdone(numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','latomdone')
  allocate(lopanyneigh(numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','lopanyneigh')
  allocate(nfreeatom(numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','nfreeatom')
  allocate(nneigh(numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','nneigh')
  allocate(rBOcutmax(numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','rBOcutmax')
  allocate(Zcn(numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','Zcn')
!
!  Reinitialisation point should maxneigh be increased
!
100 continue
  lmaxneighok = .true.
  if (allocated(rneigh)) then
    deallocate(d2BtotiR,stat=status)
    if (status/=0) call deallocate_error('bondorderd','d2BtotiR')
    deallocate(d2BtotiA,stat=status)
    if (status/=0) call deallocate_error('bondorderd','d2BtotiA')
    deallocate(d2i,stat=status)
    if (status/=0) call deallocate_error('bondorderd','d2i')
    deallocate(d1BtotiR,stat=status)
    if (status/=0) call deallocate_error('bondorderd','d1BtotiR')
    deallocate(d1BtotiA,stat=status)
    if (status/=0) call deallocate_error('bondorderd','d1BtotiA')
    deallocate(d1i,stat=status)
    if (status/=0) call deallocate_error('bondorderd','d1i')
    deallocate(zneigh,stat=status)
    if (status/=0) call deallocate_error('bondorderd','zneigh')
    deallocate(yneigh,stat=status)
    if (status/=0) call deallocate_error('bondorderd','yneigh')
    deallocate(xneigh,stat=status)
    if (status/=0) call deallocate_error('bondorderd','xneigh')
    deallocate(rneigh,stat=status)
    if (status/=0) call deallocate_error('bondorderd','rneigh')
    deallocate(ineigh,stat=status)
    if (status/=0) call deallocate_error('bondorderd','ineigh')
    deallocate(nbopotptr,stat=status)
    if (status/=0) call deallocate_error('bondorderd','nbopotptr')
    deallocate(neighno,stat=status)
    if (status/=0) call deallocate_error('bondorderd','neighno')
  endif
!
!  Initialise Bond Order energy
!
  ebondorder = 0.0_dp
!
!  Set parameter for pairwise storage memory
!
  maxneigh2 = maxneigh + maxneigh*(maxneigh + 1)/2
  maxneigh22 = maxneigh2*(maxneigh2 + 1)/2
!
!  Allocate local memory
!
  allocate(neighno(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','neighno')
  allocate(nbopotptr(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','nbopotptr')
  allocate(ineigh(3_i4,maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','ineigh')
  allocate(rneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','rneigh')
  allocate(xneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','xneigh')
  allocate(yneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','yneigh')
  allocate(zneigh(maxneigh,numat),stat=status)
  if (status/=0) call outofmemory('bondorderd','zneigh')
  if (lgrad1) then
    allocate(d1i(maxneigh2),stat=status)
    if (status/=0) call outofmemory('bondorderd','d1i')
    allocate(d1BtotiA(maxneigh2),stat=status)
    if (status/=0) call outofmemory('bondorderd','d1BtotiA')
    allocate(d1BtotiR(maxneigh2),stat=status)
    if (status/=0) call outofmemory('bondorderd','d1BtotiR')
  else
    allocate(d1i(1),stat=status)
    if (status/=0) call outofmemory('bondorderd','d1i')
    allocate(d1BtotiA(1),stat=status)
    if (status/=0) call outofmemory('bondorderd','d1BtotiA')
    allocate(d1BtotiR(1),stat=status)
    if (status/=0) call outofmemory('bondorderd','d1BtotiR')
  endif
  if (lgrad2) then
    allocate(d2i(maxneigh22),stat=status)   
    if (status/=0) call outofmemory('bondorderd','d2i')
    allocate(d2BtotiA(maxneigh22),stat=status)
    if (status/=0) call outofmemory('bondorderd','d2BtotiA')
    allocate(d2BtotiR(maxneigh22),stat=status)
    if (status/=0) call outofmemory('bondorderd','d2BtotiR')
  else
    allocate(d2i(1),stat=status)   
    if (status/=0) call outofmemory('bondorderd','d2i')
    allocate(d2BtotiA(1),stat=status)
    if (status/=0) call outofmemory('bondorderd','d2BtotiA')
    allocate(d2BtotiR(1),stat=status)
    if (status/=0) call outofmemory('bondorderd','d2BtotiR')
  endif
!****************************
!  Find list of free atoms  *
!****************************
  if (lfreeze.and..not.ldophonon) then
    ii = 0
    do i = 1,numat
      if (lopf(nrelat(i))) then
        ii = ii + 1
        nfreeatom(i) = ii
      else
        nfreeatom(i) = 0
      endif
    enddo
  else
    do i = 1,numat
      nfreeatom(i) = i
    enddo
  endif
!*************************************
!  Find cut-off radii for all atoms  *
!*************************************
  do i = 1,numat
    nati = nat(i)
    ntypi = nftype(i)
    rBOcutmax(i) = 0.0_dp
!
!  Check twobody potentials
!
    do j = 1,nbopot
      if (nati.eq.nBOspec1(j).and.(ntypi.eq.nBOtyp1(j).or.nBOtyp1(j).eq.0)) then
        rBOcutmax(i) = max(rBOcutmax(i),rBOmax(j))
      endif
      if (nati.eq.nBOspec2(j).and.(ntypi.eq.nBOtyp2(j).or.nBOtyp2(j).eq.0)) then
        rBOcutmax(i) = max(rBOcutmax(i),rBOmax(j))
      endif
    enddo
  enddo
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
  call getBOneighbour(maxneigh,rBOcutmax,nBOpotptr,nneigh,neighno,rneigh, &
                      xneigh,yneigh,zneigh,ineigh,latomdone,lmaxneighok)
!**********************************************************************
!  If maxneigh has been exceeded, increase value and return to start  *
!**********************************************************************
  if (.not.lmaxneighok) then
    do i = 1,numat
      if (nneigh(i).gt.maxneigh) maxneigh = nneigh(i)
    enddo
    if (index(keyword,'verb').ne.0.and.ioproc) then
      write(ioout,'(/,''  Increasing maxneigh to '',i6)') maxneigh
    endif
    goto 100
  endif
!*******************************
!  Sort neighbours into order  *
!******************************* 
  if (lspatialok) then
    do i = 1,numat
!               
!  Build pointer
!               
      do nn = 1,nneigh(i)
        nmin = numat + 1 
        do nn2 = nn,nneigh(i) 
          if (neighno(nn2,i).lt.nmin) then
            nmin = neighno(nn2,i)
            nptr = nn2  
          endif
        enddo       
!         
!  Sort quantities
!
        if (nptr.ne.nn) then
          itmp = neighno(nptr,i)
          neighno(nptr,i) = neighno(nn,i)
          neighno(nn,i)  = itmp
          itmp = nbopotptr(nptr,i)
          nbopotptr(nptr,i) = nbopotptr(nn,i)
          nbopotptr(nn,i)  = itmp
          itmp = ineigh(1,nptr,i)
          ineigh(1,nptr,i) = ineigh(1,nn,i)
          ineigh(1,nn,i)  = itmp
          itmp = ineigh(2,nptr,i)
          ineigh(2,nptr,i) = ineigh(2,nn,i)
          ineigh(2,nn,i)  = itmp
          itmp = ineigh(3,nptr,i)
          ineigh(3,nptr,i) = ineigh(3,nn,i)
          ineigh(3,nn,i)  = itmp
          rtmp = rneigh(nptr,i)
          rneigh(nptr,i) = rneigh(nn,i)
          rneigh(nn,i)  = rtmp
          rtmp = xneigh(nptr,i)
          xneigh(nptr,i) = xneigh(nn,i)
          xneigh(nn,i)  = rtmp
          rtmp = yneigh(nptr,i)
          yneigh(nptr,i) = yneigh(nn,i)
          yneigh(nn,i)  = rtmp
          rtmp = zneigh(nptr,i)
          zneigh(nptr,i) = zneigh(nn,i)
          zneigh(nn,i)  = rtmp
        endif  
      enddo         
    enddo
  endif
!*********************************************
!  Calculate numbers of neighbours for each  *
!*********************************************
  do i = 1,numat
!
!  Set initial value for lopanyneigh - this
!  variable indicates whether an atom has 
!  any neighbours for which derivatives are
!  required
!
    if (.not.lfreeze.or.ldophonon) then
      lopanyneigh(i) = .true.
    else
      lopanyneigh(i) = lopf(nrelat(i))
    endif
    do n = 1,nneigh(i)
      j = neighno(n,i)
!
!  Check whether atom is free to optimise
!
      if (lopf(nrelat(j))) then
        lopanyneigh(i) = .true.
      endif
    enddo
  enddo
  if (index(keyword,'debu').ne.0.and.ioproc) then
    write(ioout,'(/,''  Number of neighbours for atoms :'',/)')
    write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    do i = 1,numat
      write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nneigh(i)
    enddo
    write(ioout,'(/,''  Neighbours of atoms :'',/)')
    do i = 1,numat
      write(ioout,'(i4,8(1x,i4))') i,(neighno(nn,i),nn=1,nneigh(i))
    enddo
  endif
!
!  Initialise coordination number
!
  Zcn(1:numat) = 0.0_dp
!*****************************
!  Loop over pairs of atoms  *
!*****************************
  do i = 1,numat
!
!  Set variables relating to i
!
    nati = nat(i)
    ntypi = nftype(i)
    nregioni = nregionno(nsft + nrelat(i))
    nregiontypi = nregiontype(nregioni,ncf)
    lslicei = lsliceatom(nsft + nrelat(i))
!
    lilocal = (atom2local(i).gt.0)
!
!  Set total number of distances for neighbours of i
!
    nneighi2 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 
    nneighi22 = nneighi2*(nneighi2 + 1)/2
!
!  Initialise derivative storage for neighbours of i
!
    if (lgrad1) then
      d1i(1:nneighi2) = 0.0_dp
      if (lgrad2) then
        d2i(1:nneighi22) = 0.0_dp
      endif
    endif
!
!  Check for self energy terms for atom i
!
    if (nboZ.gt.0.and.lilocal) then
      lfound = .false.
      nboZi = 0
      do while (.not.lfound.and.nboZi.lt.nboZ)
        nboZi = nboZi + 1
        if (nBOspecZ(nboZi).eq.nati) then
          if (nBOtypZ(nboZi).eq.ntypi.or.nBOtypZ(nboZi).eq.0) then
            lfound = .true.
          endif
        endif
      enddo
      if (lfound) then
        ebondorder = ebondorder + BOecoeffZ(nboZi)
      endif
    endif
!
!  Loop over neighbours of i (=> j)
!
    ni = 1
    niloop: do while (ni.le.nneigh(i))
      j = neighno(ni,i)
!
!  Do we need to do this pair of atoms
!
      if (lopanyneigh(i).or.lopanyneigh(j)) then
!
!  Set variables relating to j
!
      natj = nat(j)
      ntypj = nftype(j)
      nregionj = nregionno(nsft + nrelat(j))
      nregiontypj = nregiontype(nregionj,ncf)
      lslicej = lsliceatom(nsft + nrelat(j))
!
!  Set up i-j quantities
!
      rij = rneigh(ni,i)
      xji = xneigh(ni,i)
      yji = yneigh(ni,i)
      zji = zneigh(ni,i)
      rrij = 1.0_dp/rij
!
      lreg2one  = .false.
      lreg2pair = .false.
      if (lseok.and.nregions(ncf).gt.1) then
        lreg2pair = (nregioni.gt.1.and.nregionj.gt.1)
        if (.not.lreg2pair) lreg2one = (nregioni.gt.1.or.nregionj.gt.1)
      endif
      lattach = (lslicei.and..not.lslicej.or.lslicej.and..not.lslicei)
!
!  Find i in neighbour list for j
!
      nj = 1
      lfound = .false.
      do while (nj.lt.nneigh(j).and..not.lfound)
        if (neighno(nj,j).eq.i) then
          xdiff = xneigh(nj,j) + xji
          ydiff = yneigh(nj,j) + yji
          zdiff = zneigh(nj,j) + zji
          lfound = ((abs(xdiff)+abs(ydiff)+abs(zdiff)).lt.1.0d-6)
        endif
        if (.not.lfound) nj = nj + 1
      enddo
!
!  Set total number of distances for neighbours of j
!
      nneighj2 = nneigh(j) + nneigh(j)*(nneigh(j) + 1)/2
      nneighj22 = nneighj2*(nneighj2 + 1)/2
!
!  Find two-body bond order potential between i and j
!
      lfound = .false.
      nboij = 0
      do while (.not.lfound.and.nboij.lt.nbopot) 
        nboij = nboij + 1
        if (nBOspec1(nboij).eq.nati.and.nBOspec2(nboij).eq.natj) then
          if ((nBOtyp1(nboij).eq.ntypi.or.nBOtyp1(nboij).eq.0).and. &
              (nBOtyp2(nboij).eq.ntypj.or.nBOtyp2(nboij).eq.0)) then
            lfound = .true.
          endif
        elseif (nBOspec1(nboij).eq.natj.and.nBOspec2(nboij).eq.nati) then
          if ((nBOtyp1(nboij).eq.ntypj.or.nBOtyp1(nboij).eq.0).and. &
              (nBOtyp2(nboij).eq.ntypi.or.nBOtyp2(nboij).eq.0)) then
            lfound = .true.
          endif
        endif
      enddo
      if (.not.lfound) nboij = 0
!
!  Find repulsive bond order potential from j to i
!
      lfound = .false.
      nboRij = 0
      do while (.not.lfound.and.nboRij.lt.nboR) 
        nboRij = nboRij + 1
        if (nBOspecR1(nboRij).eq.nati.and.nBOspecR2(nboRij).eq.natj) then
          if ((nBOtypR1(nboRij).eq.ntypi.or.nBOtypR1(nboRij).eq.0).and. &
              (nBOtypR2(nboRij).eq.ntypj.or.nBOtypR2(nboRij).eq.0)) then
            lfound = .true.
          endif
        endif
      enddo
      if (.not.lfound) nboRij = 0
!
!  Find attractive bond order potential from j to i
!
      lfound = .false.
      nboAij = 0
      do while (.not.lfound.and.nboAij.lt.nboA)
        nboAij = nboAij + 1
        if (nBOspecA1(nboAij).eq.nati.and.nBOspecA2(nboAij).eq.natj) then
          if ((nBOtypA1(nboAij).eq.ntypi.or.nBOtypA1(nboAij).eq.0).and. &
              (nBOtypA2(nboAij).eq.ntypj.or.nBOtypA2(nboAij).eq.0)) then
            lfound = .true.
          endif
        endif
      enddo
      if (.not.lfound) nboAij = 0
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
      lQMMMok = .true.
      if (QMMMmode(ncf).gt.0) then
        if (nregiontypi.eq.1.and.nregiontypj.eq.1) lQMMMok = .false.
      endif
      if (nboij.gt.0.and.lQMMMok) then
!****************************************
!  Valid two-body bond order potential  *
!****************************************
!
!  Calculate fij
!
        if (nBOtapertype(nboij).eq.2) then
          call ctaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,.false.)
        elseif (nBOtapertype(nboij).eq.3) then
          call murtytaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,.false.)
        else
          call mdftaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,.false.)
        endif
!
!  Calculate Bij - loop over all other neighbours
!
        bijsumA = 0.0_dp
        bijsumR = 0.0_dp
        if (lgrad1) then
          d1BtotiA(1:nneighi2) = 0.0_dp
          d1BtotiR(1:nneighi2) = 0.0_dp
          if (lgrad2) then
            d2BtotiA(1:nneighi22) = 0.0_dp
            d2BtotiR(1:nneighi22) = 0.0_dp
          endif
        endif
!
!  Loop over neighbours of i .ne. j 
!
        do k = 1,nneigh(i)
          npki = nbopotptr(k,i)
          if (k.ne.ni) then
            rik = rneigh(k,i)
            xki = xneigh(k,i)
            yki = yneigh(k,i)
            zki = zneigh(k,i)
!
!  Repulsive component
!
            if (nboRij.gt.0) then
              if (rik.lt.rBOmax(npki)) then
!
!  Calculate fik
!
                if (nBOtapertype(npki).eq.2) then
                  call ctaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,lgrad2,.false.)
                elseif (nBOtapertype(npki).eq.3) then
                  call murtytaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,lgrad2,.false.)
                else
                  call mdftaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,lgrad2,.false.)
                endif
!
!  Calculate Gijk
!
                if (nBOtypeR(nboRij).ne.1) then
                  call GthetaBO(xji,yji,zji,xki,yki,zki,nBOtypeR(nboRij),BOccoeffR(1,nboRij), &
                                BOhcoeffR(nboRij),Gijk,dGijkdr,d2Gijkdr2,d3Gijkdr3,lgrad1,lgrad2,.false.)
                else
                  Gijk = 1.0_dp
                  dGijkdr = 0.0_dp
                  d2Gijkdr2 = 0.0_dp
                endif
!
!  Calculate exponential factor
!
                mR = nint(BOmcoeffR(nboRij))
                if (lBOzrlR(nboRij)) then
!
!  Find repulsive bond order potential from i to k
!
                  lfound = .false.
                  nboRik = 0
                  kk = neighno(k,i)
                  natk = nat(kk)
                  ntypk = nftype(kk)
                  do while (.not.lfound.and.nboRik.lt.nboR)
                    nboRik = nboRik + 1
                    if (nBOspecR1(nboRik).eq.nati.and.nBOspecR2(nboRik).eq.natk) then
                      if ((nBOtypR1(nboRik).eq.ntypi.or.nBOtypR1(nboRik).eq.0).and. &
                          (nBOtypR2(nboRik).eq.ntypk.or.nBOtypR2(nboRik).eq.0)) then
                        lfound = .true.
                      endif
                    endif
                  enddo
                  if (.not.lfound) then
                    call outerror('no bond order found for i-k in ZRL algorithm',0_i4)
                    call stopnow('bondorderd')
                  endif
                  rlambda = (BOlcoeffR(nboRij)*rij - BOlcoeffR(nboRik)*rik)
                else
                  rlambda = BOlcoeffR(nboRij)*(rij - rik)
                endif
                expijk = exp(rlambda**mR)
!
!  Combine terms
!
                bijsumR = bijsumR + Gijk*fik*expijk
                if (lgrad1) then
!
!  Derivatives
!
!  Find index for j-k 
!
                  if (ni.ge.k) then
                    njk = nneigh(i) + ni*(ni-1)/2 + k
                  else
                    njk = nneigh(i) + k*(k-1)/2 + ni
                  endif
!
                  rrik = 1.0_dp/rik
                  dfikdr = rrik*dfikdr
                  RmR = dble(mR)
                  if (mR.ge.1) then
                    if (lBOzrlR(nboRij)) then
                      drlambdadrij = BOlcoeffR(nboRij)
                      drlambdadrik = BOlcoeffR(nboRik)
                    else
                      drlambdadrij = BOlcoeffR(nboRij)
                      drlambdadrik = BOlcoeffR(nboRij)
                    endif
                    dexpijkdr = RmR*(rlambda**(mR-1))*expijk
                  else
                    dexpijkdr = 0.0_dp
                  endif
!
                  d1BtotiR(k)   = d1BtotiR(k)   + Gijk*dfikdr*expijk
!
                  d1BtotiR(ni)  = d1BtotiR(ni)  + dGijkdr(1)*fik*expijk
                  d1BtotiR(k)   = d1BtotiR(k)   + dGijkdr(2)*fik*expijk
                  d1BtotiR(njk) = d1BtotiR(njk) + dGijkdr(3)*fik*expijk
!
                  d1BtotiR(ni) = d1BtotiR(ni) + Gijk*fik*dexpijkdr*rrij*drlambdadrij
                  d1BtotiR(k)  = d1BtotiR(k)  - Gijk*fik*dexpijkdr*rrik*drlambdadrik
!
                  if (lgrad2) then
                    d2fikdr2 = rrik*rrik*(d2fikdr2 - dfikdr)
                    if (mR.ge.2) then
                      d2expijkdr2 = expijk*(RmR*(RmR - 1.0_dp)*rlambda**(mR-2) + (RmR*rlambda**(mR-1))**2)
                    elseif (mR.eq.1) then
                      d2expijkdr2 = expijk*(RmR*rlambda**(mR-1))**2
                    else
                      d2expijkdr2 = 0.0_dp
                    endif
!
!  i-j / i-j
!
                    nn = ni*(ni + 1)/2
                    d2BtotiR(nn) = d2BtotiR(nn) + Gijk*fik*d2expijkdr2*rrij*rrij*drlambdadrij*drlambdadrij
                    d2BtotiR(nn) = d2BtotiR(nn) - Gijk*fik*dexpijkdr*rrij*rrij*rrij*drlambdadrij
                    d2BtotiR(nn) = d2BtotiR(nn) + d2Gijkdr2(1)*fik*expijk
                    d2BtotiR(nn) = d2BtotiR(nn) + 2.0_dp*dGijkdr(1)*fik*dexpijkdr*rrij*drlambdadrij
!
!  i-j / i-k
!
                    if (ni.ge.k) then
                      nn = ni*(ni - 1)/2 + k
                    else
                      nn = k*(k - 1)/2 + ni
                    endif
                    d2BtotiR(nn) = d2BtotiR(nn) + Gijk*dfikdr*dexpijkdr*rrij*drlambdadrij
                    d2BtotiR(nn) = d2BtotiR(nn) - Gijk*fik*d2expijkdr2*rrij*rrik*drlambdadrij*drlambdadrik
                    d2BtotiR(nn) = d2BtotiR(nn) + d2Gijkdr2(2)*fik*expijk
                    d2BtotiR(nn) = d2BtotiR(nn) - dGijkdr(1)*fik*dexpijkdr*rrik*drlambdadrik
                    d2BtotiR(nn) = d2BtotiR(nn) + dGijkdr(1)*dfikdr*expijk
                    d2BtotiR(nn) = d2BtotiR(nn) + dGijkdr(2)*fik*dexpijkdr*rrij*drlambdadrij
!
!  i-j / j-k
!
                    if (ni.ge.njk) then
                      nn = ni*(ni - 1)/2 + njk
                    else
                      nn = njk*(njk - 1)/2 + ni
                    endif
                    d2BtotiR(nn) = d2BtotiR(nn) + d2Gijkdr2(3)*fik*expijk
                    d2BtotiR(nn) = d2BtotiR(nn) + dGijkdr(3)*fik*dexpijkdr*rrij*drlambdadrij
!
!  i-k / i-k
!
                    nn = k*(k + 1)/2
                    d2BtotiR(nn) = d2BtotiR(nn) + Gijk*d2fikdr2*expijk
                    d2BtotiR(nn) = d2BtotiR(nn) - 2.0_dp*Gijk*dfikdr*dexpijkdr*rrik*drlambdadrik
                    d2BtotiR(nn) = d2BtotiR(nn) + Gijk*fik*d2expijkdr2*rrik*rrik*drlambdadrik*drlambdadrik
                    d2BtotiR(nn) = d2BtotiR(nn) + Gijk*fik*dexpijkdr*rrik*rrik*rrik*drlambdadrik
                    d2BtotiR(nn) = d2BtotiR(nn) + 2.0_dp*dGijkdr(2)*dfikdr*expijk
                    d2BtotiR(nn) = d2BtotiR(nn) + d2Gijkdr2(4)*fik*expijk
                    d2BtotiR(nn) = d2BtotiR(nn) - 2.0_dp*dGijkdr(2)*fik*dexpijkdr*rrik*drlambdadrik
!
!  i-k / j-k
!
                    if (k.ge.njk) then
                      nn = k*(k - 1)/2 + njk
                    else
                      nn = njk*(njk - 1)/2 + k
                    endif
                    d2BtotiR(nn) = d2BtotiR(nn) + d2Gijkdr2(5)*fik*expijk
                    d2BtotiR(nn) = d2BtotiR(nn) + dGijkdr(3)*dfikdr*expijk
                    d2BtotiR(nn) = d2BtotiR(nn) - dGijkdr(3)*fik*dexpijkdr*rrik*drlambdadrik
!
!  j-k / j-k
!
                    nn = njk*(njk + 1)/2
                    d2BtotiR(nn) = d2BtotiR(nn) + d2Gijkdr2(6)*fik*expijk
                  endif
                endif
              endif
            endif
!
!  Attractive component
!
            if (nboAij.gt.0) then
              if (rik.lt.rBOmax(npki)) then
!
!  Calculate fik
!
                if (nBOtapertype(npki).eq.2) then
                  call ctaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,lgrad2,.false.)
                elseif (nBOtapertype(npki).eq.3) then
                  call murtytaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,lgrad2,.false.)
                else
                  call mdftaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,lgrad2,.false.)
                endif
!
!  Calculate Gijk
!
                if (nBOtypeA(nboAij).ne.1) then
                  call GthetaBO(xji,yji,zji,xki,yki,zki,nBOtypeA(nboAij),BOccoeffA(1,nboAij), &
                                BOhcoeffA(nboAij),Gijk,dGijkdr,d2Gijkdr2,d3Gijkdr3,lgrad1,lgrad2,.false.)
                else
                  Gijk = 1.0_dp
                  dGijkdr = 0.0_dp
                  d2Gijkdr2 = 0.0_dp
                endif
!
!  Calculate exponential factor
!
                mA = nint(BOmcoeffA(nboAij))
                if (lBOzrlA(nboAij)) then
!
!  Find attractive bond order potential from i to k
!
                  lfound = .false.
                  nboAik = 0
                  kk = neighno(k,i)
                  natk = nat(kk)
                  ntypk = nftype(kk)
                  do while (.not.lfound.and.nboAik.lt.nboA)
                    nboAik = nboAik + 1
                    if (nBOspecA1(nboAik).eq.nati.and.nBOspecA2(nboAik).eq.natk) then
                      if ((nBOtypA1(nboAik).eq.ntypi.or.nBOtypA1(nboAik).eq.0).and. &
                          (nBOtypA2(nboAik).eq.ntypk.or.nBOtypA2(nboAik).eq.0)) then
                        lfound = .true.
                      endif
                    endif
                  enddo
                  if (.not.lfound) then
                    call outerror('no bond order found for i-k in ZRL algorithm',0_i4)
                    call stopnow('bondorderd')
                  endif
                  rlambda = (BOlcoeffA(nboAij)*rij - BOlcoeffA(nboAik)*rik)
                else
                  rlambda = BOlcoeffA(nboAij)*(rij - rik)
                endif
                expijk = exp(rlambda**mA)
!
!  Combine terms
!
                bijsumA = bijsumA + Gijk*fik*expijk
                if (lgrad1) then
!
!  Derivatives
!
!  Find index for j-k 
!
                  if (ni.ge.k) then
                    njk = nneigh(i) + ni*(ni-1)/2 + k
                  else
                    njk = nneigh(i) + k*(k-1)/2 + ni
                  endif
!
                  rrik = 1.0_dp/rik
                  dfikdr = rrik*dfikdr
                  RmA = dble(mA)
                  if (mA.ge.1) then
                    if (lBOzrlA(nboAij)) then
                      drlambdadrij = BOlcoeffA(nboAij)
                      drlambdadrik = BOlcoeffA(nboAik)
                    else
                      drlambdadrij = BOlcoeffA(nboAij)
                      drlambdadrik = BOlcoeffA(nboAij)
                    endif
                    dexpijkdr = RmA*(rlambda**(mA-1))*expijk
                  else
                    dexpijkdr = 0.0_dp
                  endif
!
                  d1BtotiA(k)   = d1BtotiA(k)   + Gijk*dfikdr*expijk
!
                  d1BtotiA(ni)  = d1BtotiA(ni)  + dGijkdr(1)*fik*expijk
                  d1BtotiA(k)   = d1BtotiA(k)   + dGijkdr(2)*fik*expijk
                  d1BtotiA(njk) = d1BtotiA(njk) + dGijkdr(3)*fik*expijk
!
                  d1BtotiA(ni) = d1BtotiA(ni) + Gijk*fik*dexpijkdr*rrij*drlambdadrij
                  d1BtotiA(k)  = d1BtotiA(k)  - Gijk*fik*dexpijkdr*rrik*drlambdadrik
!
                  if (lgrad2) then
                    d2fikdr2 = rrik*rrik*(d2fikdr2 - dfikdr)
                    if (mA.ge.2) then
                      d2expijkdr2 = expijk*(RmA*(RmA - 1.0_dp)*rlambda**(mA-2) + (RmA*rlambda**(mA-1))**2)
                    elseif (mA.eq.1) then
                      d2expijkdr2 = expijk*(RmA*rlambda**(mA-1))**2
                    else
                      d2expijkdr2 = 0.0_dp
                    endif
!
!  i-j / i-j
!
                    nn = ni*(ni + 1)/2
                    d2BtotiA(nn) = d2BtotiA(nn) + Gijk*fik*d2expijkdr2*rrij*rrij*drlambdadrij*drlambdadrij
                    d2BtotiA(nn) = d2BtotiA(nn) - Gijk*fik*dexpijkdr*rrij*rrij*rrij*drlambdadrij
                    d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(1)*fik*expijk
                    d2BtotiA(nn) = d2BtotiA(nn) + 2.0_dp*dGijkdr(1)*fik*dexpijkdr*rrij*drlambdadrij
!
!  i-j / i-k
!
                    if (ni.ge.k) then
                      nn = ni*(ni - 1)/2 + k
                    else
                      nn = k*(k - 1)/2 + ni
                    endif
                    d2BtotiA(nn) = d2BtotiA(nn) + Gijk*dfikdr*dexpijkdr*rrij*drlambdadrij
                    d2BtotiA(nn) = d2BtotiA(nn) - Gijk*fik*d2expijkdr2*rrij*rrik*drlambdadrij*drlambdadrik
                    d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(2)*fik*expijk
                    d2BtotiA(nn) = d2BtotiA(nn) - dGijkdr(1)*fik*dexpijkdr*rrik*drlambdadrik
                    d2BtotiA(nn) = d2BtotiA(nn) + dGijkdr(1)*dfikdr*expijk
                    d2BtotiA(nn) = d2BtotiA(nn) + dGijkdr(2)*fik*dexpijkdr*rrij*drlambdadrij
!
!  i-j / j-k
!
                    if (ni.ge.njk) then
                      nn = ni*(ni - 1)/2 + njk
                    else
                      nn = njk*(njk - 1)/2 + ni
                    endif
                    d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(3)*fik*expijk
                    d2BtotiA(nn) = d2BtotiA(nn) + dGijkdr(3)*fik*dexpijkdr*rrij*drlambdadrij
!
!  i-k / i-k
!
                    nn = k*(k + 1)/2
                    d2BtotiA(nn) = d2BtotiA(nn) + Gijk*d2fikdr2*expijk
                    d2BtotiA(nn) = d2BtotiA(nn) - 2.0_dp*Gijk*dfikdr*dexpijkdr*rrik*drlambdadrik
                    d2BtotiA(nn) = d2BtotiA(nn) + Gijk*fik*d2expijkdr2*rrik*rrik*drlambdadrik*drlambdadrik
                    d2BtotiA(nn) = d2BtotiA(nn) + Gijk*fik*dexpijkdr*rrik*rrik*rrik*drlambdadrik
                    d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(4)*fik*expijk
                    d2BtotiA(nn) = d2BtotiA(nn) + 2.0_dp*dGijkdr(2)*dfikdr*expijk
                    d2BtotiA(nn) = d2BtotiA(nn) - 2.0_dp*dGijkdr(2)*fik*dexpijkdr*rrik*drlambdadrik
!
!  i-k / j-k
!
                    if (k.ge.njk) then
                      nn = k*(k - 1)/2 + njk
                    else
                      nn = njk*(njk - 1)/2 + k
                    endif
                    d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(5)*fik*expijk
                    d2BtotiA(nn) = d2BtotiA(nn) + dGijkdr(3)*dfikdr*expijk
                    d2BtotiA(nn) = d2BtotiA(nn) - dGijkdr(3)*fik*dexpijkdr*rrik*drlambdadrik
!
!  j-k / j-k
!
                    nn = njk*(njk + 1)/2
                    d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(6)*fik*expijk
                  endif
                endif
              endif
            endif
          endif
        enddo
!
        if (nboAij.gt.0) then
          BOncoAij = BOncoeffA(nboAij)
          zAi3 = BOecoeffA(nboAij)**BOncoAij
        else
          BOncoAij = 1.0_dp
          zAi3 = 0.0_dp
        endif
        if (nboRij.gt.0) then
          BOncoRij = BOncoeffR(nboRij)
          zRi3 = BOecoeffR(nboRij)**BOncoRij
        else
          BOncoRij = 1.0_dp
          zRi3 = 0.0_dp
        endif
!
        if (abs(bijsumA).gt.1.0d-12) then
          bijsumAn2 = bijsumA**(BOncoAij - 2.0_dp)
        else
          bijsumAn2 = 0.0_dp
        endif
        if (abs(bijsumR).gt.1.0d-12) then
          bijsumRn2 = bijsumR**(BOncoRij - 2.0_dp)
        else
          bijsumRn2 = 0.0_dp
        endif
!
        bijsumAn1 = bijsumAn2*bijsumA
        bijsumRn1 = bijsumRn2*bijsumR
!
        bijsumA1 = 1.0_dp + zAi3*bijsumA*bijsumAn1
        bijsumR1 = 1.0_dp + zRi3*bijsumR*bijsumRn1
!
        rbijsumA1 = 1.0_dp/bijsumA1
        rbijsumR1 = 1.0_dp/bijsumR1
!
        bijA = bijsumA1**(-0.5_dp/BOncoAij)
        bijR = bijsumR1**(-0.5_dp/BOncoRij)
!
!  Scale derivatives by bijsum factors
!
        if (lgrad1) then
          if (lgrad2) then
!
!  Second derivatives
!
            if (bijsumA.gt.0.0_dp) then
              rtmp  = 0.25_dp*zAi3*bijA*rbijsumA1*bijsumAn1
              do nn = 1,nneighi22
                d2BtotiA(nn) = - rtmp*d2BtotiA(nn)
              enddo
              rtmp = 0.25_dp*zAi3*bijA*rbijsumA1*(zAi3*(1.0_dp + 0.5_dp/BOncoAij)*(bijsumAn1**2.0_dp)* &
                BOncoAij*rbijsumA1 - ((BOncoAij-1)*bijsumAn2))
              nn = 0 
              do nn1 = 1,nneighi2
                do nn2 = 1,nn1
                  nn = nn + 1
                  d2BtotiA(nn) = d2BtotiA(nn) + rtmp*d1BtotiA(nn2)*d1BtotiA(nn1)
                enddo
              enddo
            endif
            if (bijsumR.gt.0.0_dp) then
              rtmp = 0.25_dp*zRi3*bijR*rbijsumR1*bijsumRn1
              do nn = 1,nneighi22
                d2BtotiR(nn) = - rtmp*d2BtotiR(nn)
              enddo
              rtmp = 0.25_dp*zRi3*bijR*rbijsumR1*(zRi3*(1.0_dp + 0.5_dp/BOncoRij)*(bijsumRn1**2.0_dp)* &
                BOncoRij*rbijsumR1 - ((BOncoRij-1)*bijsumRn2))
              nn = 0 
              do nn1 = 1,nneighi2
                do nn2 = 1,nn1
                  nn = nn + 1
                  d2BtotiR(nn) = d2BtotiR(nn) + rtmp*d1BtotiR(nn2)*d1BtotiR(nn1)
                enddo
              enddo
            endif
          endif
!
!  First derivatives
!
          if (bijsumA.gt.0.0_dp) then
            rtmp = 0.25_dp*zAi3*bijsumAn1*bijA*rbijsumA1
            do nn = 1,nneighi2
              d1BtotiA(nn) = - rtmp*d1BtotiA(nn)
            enddo
          endif
          if (bijsumR.gt.0.0_dp) then
            rtmp = 0.25_dp*zRi3*bijsumRn1*bijR*rbijsumR1
            do nn = 1,nneighi2
              d1BtotiR(nn) = - rtmp*d1BtotiR(nn)
            enddo
          endif
        endif
!
!  Calculate two-body component of potential
!
        Vr = BOacoeff(nboij)*exp(-BOzacoeff(nboij)*rij)
        Va = BObcoeff(nboij)*exp(-BOzbcoeff(nboij)*rij)
        if (lgrad1) then
          dVrdr = - BOzacoeff(nboij)*Vr
          dVadr = - BOzbcoeff(nboij)*Va
!
          dVrdr = rrij*dVrdr
          dVadr = rrij*dVadr
          if (lgrad2) then
            d2Vrdr2 = BOzacoeff(nboij)*BOzacoeff(nboij)*Vr
            d2Vadr2 = BOzbcoeff(nboij)*BOzbcoeff(nboij)*Va
            d2Vrdr2 = rrij*rrij*(d2Vrdr2 - dVrdr)
            d2Vadr2 = rrij*rrij*(d2Vadr2 - dVadr)
          endif
        endif
!
!  Calculate total i-j potential
!
        BtotA = 0.5_dp*bijA
        BtotR = 0.5_dp*bijR
        if (lilocal) then
          eij = f*(BtotR*Vr - BtotA*Va)
!
!  Add to surface energy totals if appropriate
!
          if (lseok) then
            if (lreg2one) then
              esregion12 = esregion12 + eij
            elseif (lreg2pair) then
              esregion2 = esregion2 + eij
            else
              ebondorder = ebondorder + eij
            endif
          else
            ebondorder = ebondorder + eij
          endif
          if (lattach) eattach = eattach + eij
        endif
!
!  Add contribution to coordination number
!
        Zcn(i) = Zcn(i) + f*bijA
!
!  Derivatives of Bond Order potential energy
!
        if (lgrad1) then
          dfdr = rrij*dfdr
          d1i(ni) = d1i(ni) + dfdr*(BtotR*Vr - BtotA*Va)
          d1i(ni) = d1i(ni) + f*(BtotR*dVrdr - BtotA*dVadr)
          do nn = 1,nneighi2
            d1i(nn) = d1i(nn) + f*(Vr*d1BtotiR(nn) - Va*d1BtotiA(nn))
          enddo
          if (lgrad2) then
            ind = ni*(ni + 1)/2
            d2fdr2 = rrij*rrij*(d2fdr2 - dfdr)
            d2i(ind) = d2i(ind) + d2fdr2*(BtotR*Vr - BtotA*Va)
            d2i(ind) = d2i(ind) + 2.0_dp*dfdr*(BtotR*dVrdr - BtotA*dVadr)
            d2i(ind) = d2i(ind) + f*(BtotR*d2Vrdr2 - BtotA*d2Vadr2)
            do nn = 1,nneighi2
              if (ni.ge.nn) then
                ind = ni*(ni - 1)/2 + nn
              else
                ind = nn*(nn - 1)/2 + ni
              endif
              d2i(ind) = d2i(ind) + ((dfdr*Vr + f*dVrdr)*d1BtotiR(nn) - (dfdr*Va + f*dVadr)*d1BtotiA(nn))
              if (ni.eq.nn) then
                d2i(ind) = d2i(ind) + ((dfdr*Vr + f*dVrdr)*d1BtotiR(nn) - (dfdr*Va + f*dVadr)*d1BtotiA(nn))
              endif
            enddo
            do nn = 1,nneighi22
              d2i(nn) = d2i(nn) + f*(Vr*d2BtotiR(nn) - Va*d2BtotiA(nn))
            enddo
          endif
        endif
      endif
!
!  End condition section on i or j being associated with moving atom
!
      endif
!
!  End of loop over neighbours of i
!
      ni = ni + 1
    enddo niloop
!
!  Add derivatives due to neighbours of i
!
    if (lgrad1) then
      if (lilocal) then
        call d1add(i,maxneigh,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nboatomRptr,d1i,.false.,.true.)
      endif
      if (lgrad2) then
        if (ldophonon) then
          call d2addpdm(i,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nboatomRptr,d1i,d2i,xkv,ykv,zkv,.false.)
        else
          call d2adddm(i,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nboatomRptr,d1i,d2i,.false.)
        endif
      endif
    endif
  enddo
  if (nboZ.gt.0) then
!
!  Allocate memory that is required for second derivative term
!
    if (lgrad2) then
      allocate(d1ZiCN(maxneigh2,maxneigh),stat=status)
      if (status/=0) call outofmemory('bondorderd','d1ZiCN')
    endif
!**************************************
!  Compute coordination contribution  *
!**************************************
    do i = 1,numat
!
!  Set variables relating to i
!
      nati = nat(i)
      ntypi = nftype(i)
!
      lilocal = (atom2local(i).gt.0)
!
!  Check for coordination terms for atom i
!
      lfound = .false.
      nboZi = 0
      do while (.not.lfound.and.nboZi.lt.nboZ)
        nboZi = nboZi + 1
        if (nBOspecZ(nboZi).eq.nati) then
          if (nBOtypZ(nboZi).eq.ntypi.or.nBOtypZ(nboZi).eq.0) then
            lfound = .true.
          endif
        endif
      enddo
      if (lfound) then
        dZ = Zcn(i) - BOzcoeffZ(nboZi)
        zsign = sign(1.0_dp,dZ)
        call botaper(dZ,fz,dfzdz,d2fzdz2,d3fzdz3,lgrad1,lgrad2,.false.)
        deltaZ = zsign*fz
        lneedBOcnderiv = (abs(dble(nint(deltaZ)) - deltaZ).gt.BOcntol)
        if (lilocal) then
          ebondorder = ebondorder + deltaZ*(BOccoeffZ(1,nboZi) + BOccoeffZ(2,nboZi)*deltaZ)
        endif
!
!  If deltaZ is not an integer then compute derivatives
!
        if (lneedBOcnderiv.and.lgrad1) then
          dedZ = (BOccoeffZ(1,nboZi) + 2.0_dp*BOccoeffZ(2,nboZi)*deltaZ)*zsign*dfzdz
          if (lgrad2) then
            d2edZ2 = (BOccoeffZ(1,nboZi) + 2.0_dp*BOccoeffZ(2,nboZi)*deltaZ)*zsign*d2fzdz2 + &
                     2.0_dp*BOccoeffZ(2,nboZi)*(zsign*dfzdz)**2
          endif
!
          nregioni = nregionno(nsft + nrelat(i))
          nregiontypi = nregiontype(nregioni,ncf)
!
!  Set total number of distances for neighbours of i
!
          nneighi2 = nneigh(i) + nneigh(i)*(nneigh(i) + 1)/2 
          nneighi22 = nneighi2*(nneighi2 + 1)/2
!
!  Initialise derivative storage for neighbours of i
!
          d1i(1:nneighi2) = 0.0_dp
          if (lgrad2) then
            d2i(1:nneighi22) = 0.0_dp
            d1ZiCN(1:nneighi2,1:nneigh(i)) = 0.0_dp
          endif
!
!  Loop over neighbours of i 
!
          ni = 1
          niloopnc: do while (ni.le.nneigh(i))
            j = neighno(ni,i)
!
!  Do we need to do this pair of atoms
!
            if (lopanyneigh(i).or.lopanyneigh(j)) then
!
!  Set variables relating to j
!
              natj = nat(j)
              ntypj = nftype(j)
              nregionj = nregionno(nsft + nrelat(j))
              nregiontypj = nregiontype(nregionj,ncf)
!
!  Set up i-j quantities
!
              rij = rneigh(ni,i)
              xji = xneigh(ni,i)
              yji = yneigh(ni,i)
              zji = zneigh(ni,i)
              rrij = 1.0_dp/rij
!
!  Find i in neighbour list for j
!
              nj = 1
              lfound = .false.
              do while (nj.lt.nneigh(j).and..not.lfound)
                if (neighno(nj,j).eq.i) then
                  xdiff = xneigh(nj,j) + xji
                  ydiff = yneigh(nj,j) + yji
                  zdiff = zneigh(nj,j) + zji
                  lfound = ((abs(xdiff)+abs(ydiff)+abs(zdiff)).lt.1.0d-6)
                endif
                if (.not.lfound) nj = nj + 1
              enddo
!
!  Find two-body bond order potential between i and j
!
              lfound = .false.
              nboij = 0
              do while (.not.lfound.and.nboij.lt.nbopot) 
                nboij = nboij + 1
                if (nBOspec1(nboij).eq.nati.and.nBOspec2(nboij).eq.natj) then
                  if ((nBOtyp1(nboij).eq.ntypi.or.nBOtyp1(nboij).eq.0).and. &
                      (nBOtyp2(nboij).eq.ntypj.or.nBOtyp2(nboij).eq.0)) then
                    lfound = .true.
                  endif
                elseif (nBOspec1(nboij).eq.natj.and.nBOspec2(nboij).eq.nati) then
                  if ((nBOtyp1(nboij).eq.ntypj.or.nBOtyp1(nboij).eq.0).and. &
                      (nBOtyp2(nboij).eq.ntypi.or.nBOtyp2(nboij).eq.0)) then
                    lfound = .true.
                  endif
                endif
              enddo
              if (.not.lfound) nboij = 0
!
!  Find attractive bond order potential from j to i
!
              lfound = .false.
              nboAij = 0
              do while (.not.lfound.and.nboAij.lt.nboA)
                nboAij = nboAij + 1
                if (nBOspecA1(nboAij).eq.nati.and.nBOspecA2(nboAij).eq.natj) then
                  if ((nBOtypA1(nboAij).eq.ntypi.or.nBOtypA1(nboAij).eq.0).and. &
                      (nBOtypA2(nboAij).eq.ntypj.or.nBOtypA2(nboAij).eq.0)) then
                    lfound = .true.
                  endif
                endif
              enddo
              if (.not.lfound) nboAij = 0
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
              lQMMMok = .true.
              if (QMMMmode(ncf).gt.0) then
                if (nregiontypi.eq.1.and.nregiontypj.eq.1) lQMMMok = .false.
              endif
              if (nboij.gt.0.and.lQMMMok) then
!****************************************
!  Valid two-body bond order potential  *
!****************************************
!
!  Calculate fij
!
                if (nBOtapertype(nboij).eq.2) then
                  call ctaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,.false.)
                elseif (nBOtapertype(nboij).eq.3) then
                  call murtytaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,.false.)
                else
                  call mdftaper(rij,rBOmin(nboij),rBOmax(nboij),f,dfdr,d2fdr2,d3fdr3,lgrad1,lgrad2,.false.)
                endif
!
!  Calculate Bij and Bji - loop over all other neighbours
!
                bijsumA = 0.0_dp
                d1BtotiA(1:nneighi2) = 0.0_dp
                if (lgrad2) then
                  d2BtotiA(1:nneighi22) = 0.0_dp
                endif
!
!  Loop over neighbours of i .ne. j 
!
                do k = 1,nneigh(i)
                  npki = nbopotptr(k,i)
                  if (k.ne.ni) then
                    rik = rneigh(k,i)
                    xki = xneigh(k,i)
                    yki = yneigh(k,i)
                    zki = zneigh(k,i)
!
!  Attractive component
!
                    if (nboAij.gt.0) then
                      if (rik.lt.rBOmax(npki)) then
!
!  Calculate fik
!
                        if (nBOtapertype(npki).eq.2) then
                          call ctaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,lgrad2,.false.)
                        elseif (nBOtapertype(npki).eq.3) then
                          call murtytaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,lgrad2,.false.)
                        else
                          call mdftaper(rik,rBOmin(npki),rBOmax(npki),fik,dfikdr,d2fikdr2,d3fikdr3,lgrad1,lgrad2,.false.)
                        endif
!
!  Calculate Gijk
!
                        if (nBOtypeA(nboAij).ne.1) then
                          call GthetaBO(xji,yji,zji,xki,yki,zki,nBOtypeA(nboAij),BOccoeffA(1,nboAij), &
                                        BOhcoeffA(nboAij),Gijk,dGijkdr,d2Gijkdr2,d3Gijkdr3,lgrad1,lgrad2,.false.)
                        else
                          Gijk = 1.0_dp
                          dGijkdr = 0.0_dp
                          d2Gijkdr2 = 0.0_dp
                        endif
!
!  Calculate exponential factor
!
                        mA = nint(BOmcoeffA(nboAij))
                        if (lBOzrlA(nboAij)) then
!
!  Find attractive bond order potential from i to k
!
                          lfound = .false.
                          nboAik = 0
                          kk = neighno(k,i)
                          natk = nat(kk)
                          ntypk = nftype(kk)
                          do while (.not.lfound.and.nboAik.lt.nboA)
                            nboAik = nboAik + 1
                            if (nBOspecA1(nboAik).eq.nati.and.nBOspecA2(nboAik).eq.natk) then
                              if ((nBOtypA1(nboAik).eq.ntypi.or.nBOtypA1(nboAik).eq.0).and. &
                                  (nBOtypA2(nboAik).eq.ntypk.or.nBOtypA2(nboAik).eq.0)) then
                                lfound = .true.
                              endif
                            endif
                          enddo
                          if (.not.lfound) then
                            call outerror('no bond order found for i-k in ZRL algorithm',0_i4)
                            call stopnow('bondorderd')
                          endif
                          rlambda = (BOlcoeffA(nboAij)*rij - BOlcoeffA(nboAik)*rik)
                        else
                          rlambda = BOlcoeffA(nboAij)*(rij - rik)
                        endif
                        expijk = exp(rlambda**mA)
!
!  Combine terms
!
                        bijsumA = bijsumA + Gijk*fik*expijk
!
!  Derivatives
!
!  Find index for j-k 
!
                        if (ni.ge.k) then
                          njk = nneigh(i) + ni*(ni-1)/2 + k
                        else
                          njk = nneigh(i) + k*(k-1)/2 + ni
                        endif
!
                        rrik = 1.0_dp/rik
                        dfikdr = rrik*dfikdr
                        RmA = dble(mA)
                        if (mA.ge.1) then
                          if (lBOzrlA(nboAij)) then
                            drlambdadrij = BOlcoeffA(nboAij)
                            drlambdadrik = BOlcoeffA(nboAik)
                          else
                            drlambdadrij = BOlcoeffA(nboAij)
                            drlambdadrik = BOlcoeffA(nboAij)
                          endif
                          dexpijkdr = RmA*(rlambda**(mA-1))*expijk
                        else
                          dexpijkdr = 0.0_dp
                        endif
!
                        d1BtotiA(k)   = d1BtotiA(k)   + Gijk*dfikdr*expijk
!
                        d1BtotiA(ni)  = d1BtotiA(ni)  + dGijkdr(1)*fik*expijk
                        d1BtotiA(k)   = d1BtotiA(k)   + dGijkdr(2)*fik*expijk
                        d1BtotiA(njk) = d1BtotiA(njk) + dGijkdr(3)*fik*expijk
!
                        d1BtotiA(ni) = d1BtotiA(ni) + Gijk*fik*dexpijkdr*rrij*drlambdadrij
                        d1BtotiA(k)  = d1BtotiA(k)  - Gijk*fik*dexpijkdr*rrik*drlambdadrik
!
                        if (lgrad2) then
                          d2fikdr2 = rrik*rrik*(d2fikdr2 - dfikdr)
                          if (mA.ge.2) then
                            d2expijkdr2 = expijk*(RmA*(RmA - 1.0_dp)*rlambda**(mA-2) + (RmA*rlambda**(mA-1))**2)
                          elseif (mA.eq.1) then
                            d2expijkdr2 = expijk*(RmA*rlambda**(mA-1))**2
                          else
                            d2expijkdr2 = 0.0_dp
                          endif
!
!  i-j / i-j
!
                          nn = ni*(ni + 1)/2
                          d2BtotiA(nn) = d2BtotiA(nn) + Gijk*fik*d2expijkdr2*rrij*rrij*drlambdadrij*drlambdadrij
                          d2BtotiA(nn) = d2BtotiA(nn) - Gijk*fik*dexpijkdr*rrij*rrij*rrij*drlambdadrij
                          d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(1)*fik*expijk
                          d2BtotiA(nn) = d2BtotiA(nn) + 2.0_dp*dGijkdr(1)*fik*dexpijkdr*rrij*drlambdadrij
!
!  i-j / i-k
!
                          if (ni.ge.k) then
                            nn = ni*(ni - 1)/2 + k
                          else
                            nn = k*(k - 1)/2 + ni
                          endif
                          d2BtotiA(nn) = d2BtotiA(nn) + Gijk*dfikdr*dexpijkdr*rrij*drlambdadrij
                          d2BtotiA(nn) = d2BtotiA(nn) - Gijk*fik*d2expijkdr2*rrij*rrik*drlambdadrij*drlambdadrik
                          d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(2)*fik*expijk
                          d2BtotiA(nn) = d2BtotiA(nn) - dGijkdr(1)*fik*dexpijkdr*rrik*drlambdadrik
                          d2BtotiA(nn) = d2BtotiA(nn) + dGijkdr(1)*dfikdr*expijk
                          d2BtotiA(nn) = d2BtotiA(nn) + dGijkdr(2)*fik*dexpijkdr*rrij*drlambdadrij
!
!  i-j / j-k
!
                          if (ni.ge.njk) then
                            nn = ni*(ni - 1)/2 + njk
                          else
                            nn = njk*(njk - 1)/2 + ni
                          endif
                          d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(3)*fik*expijk
                          d2BtotiA(nn) = d2BtotiA(nn) + dGijkdr(3)*fik*dexpijkdr*rrij*drlambdadrij
!
!  i-k / i-k
!
                          nn = k*(k + 1)/2
                          d2BtotiA(nn) = d2BtotiA(nn) + Gijk*d2fikdr2*expijk
                          d2BtotiA(nn) = d2BtotiA(nn) - 2.0_dp*Gijk*dfikdr*dexpijkdr*rrik*drlambdadrik
                          d2BtotiA(nn) = d2BtotiA(nn) + Gijk*fik*d2expijkdr2*rrik*rrik*drlambdadrik*drlambdadrik
                          d2BtotiA(nn) = d2BtotiA(nn) + Gijk*fik*dexpijkdr*rrik*rrik*rrik*drlambdadrik
                          d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(4)*fik*expijk
                          d2BtotiA(nn) = d2BtotiA(nn) + 2.0_dp*dGijkdr(2)*dfikdr*expijk
                          d2BtotiA(nn) = d2BtotiA(nn) - 2.0_dp*dGijkdr(2)*fik*dexpijkdr*rrik*drlambdadrik
!
!  i-k / j-k
!
                          if (k.ge.njk) then
                            nn = k*(k - 1)/2 + njk
                          else
                            nn = njk*(njk - 1)/2 + k
                          endif
                          d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(5)*fik*expijk
                          d2BtotiA(nn) = d2BtotiA(nn) + dGijkdr(3)*dfikdr*expijk
                          d2BtotiA(nn) = d2BtotiA(nn) - dGijkdr(3)*fik*dexpijkdr*rrik*drlambdadrik
!
!  j-k / j-k
!
                          nn = njk*(njk + 1)/2
                          d2BtotiA(nn) = d2BtotiA(nn) + d2Gijkdr2(6)*fik*expijk
                        endif
                      endif
                    endif
                  endif
                enddo
!
                if (nboAij.gt.0) then
                  BOncoAij = BOncoeffA(nboAij)
                  zAi3 = BOecoeffA(nboAij)**BOncoAij
                else
                  BOncoAij = 1.0_dp
                  zAi3 = 0.0_dp
                endif
!
                if (abs(bijsumA).gt.1.0d-12) then
                  bijsumAn2 = bijsumA**(BOncoAij - 2.0_dp)
                else
                  bijsumAn2 = 0.0_dp
                endif
!
                bijsumAn1 = bijsumAn2*bijsumA
!
                bijsumA1 = 1.0_dp + zAi3*bijsumA*bijsumAn1
!
                rbijsumA1 = 1.0_dp/bijsumA1
!
                bijA = bijsumA1**(-0.5_dp/BOncoAij)
!
!  Scale derivatives by bijsum factors
!
                if (lgrad2) then
!
!  Second derivatives
!
                  if (bijsumA.gt.0.0_dp) then
                    rtmp  = 0.25_dp*zAi3*bijA*rbijsumA1*bijsumAn1
                    do nn = 1,nneighi22
                      d2BtotiA(nn) = - rtmp*d2BtotiA(nn)
                    enddo
                    rtmp = 0.25_dp*zAi3*bijA*rbijsumA1*(zAi3*(1.0_dp + 0.5_dp/BOncoAij)*(bijsumAn1**2.0_dp)* &
                      BOncoAij*rbijsumA1 - ((BOncoAij-1)*bijsumAn2))
                    nn = 0 
                    do nn1 = 1,nneighi2
                      do nn2 = 1,nn1
                        nn = nn + 1
                        d2BtotiA(nn) = d2BtotiA(nn) + rtmp*d1BtotiA(nn2)*d1BtotiA(nn1)
                      enddo
                    enddo
                  endif
                endif
!
!  First derivatives
!
                if (bijsumA.gt.0.0_dp) then
                  rtmp = 0.25_dp*zAi3*bijsumAn1*bijA*rbijsumA1
                  do nn = 1,nneighi2
                    d1BtotiA(nn) = - rtmp*d1BtotiA(nn)
                  enddo
                endif
!
!  Derivatives of Bond Order coordination energy
!
                dfdr = rrij*dfdr
                d1i(ni) = d1i(ni) + dedZ*dfdr*bijA
                do nn = 1,nneighi2
                  d1i(nn) = d1i(nn) + 2.0_dp*dedZ*f*d1BtotiA(nn)
                enddo
                if (lgrad2) then
!
!  For second derivatives it is necessary to store first derivatives of Zi for all neighbours
!
                  do nn = 1,nneighi2
                    d1ZiCN(nn,ni) = 2.0_dp*f*d1BtotiA(nn)
                  enddo
                  d1ZiCN(ni,ni) = d1ZiCN(ni,ni) + dfdr*bijA
!
!  Do second derivatives that only depend on i-j pair
!
                  ind = ni*(ni + 1)/2
                  d2fdr2 = rrij*rrij*(d2fdr2 - dfdr)
                  d2i(ind) = d2i(ind) + dedZ*d2fdr2*bijA
                  do nn = 1,nneighi2
                    if (ni.ge.nn) then
                      ind = ni*(ni - 1)/2 + nn
                    else
                      ind = nn*(nn - 1)/2 + ni
                    endif
                    d2i(ind) = d2i(ind) + 2.0_dp*dedZ*dfdr*d1BtotiA(nn)
                    if (ni.eq.nn) then
                      d2i(ind) = d2i(ind) + 2.0_dp*dedZ*dfdr*d1BtotiA(nn)
                    endif
                  enddo
                  do nn = 1,nneighi22
                    d2i(nn) = d2i(nn) + 2.0_dp*dedZ*f*d2BtotiA(nn)
                  enddo
                endif
              endif
!
!  End condition section on i or j being associated with moving atom
!
            endif
!
!  End of loop over neighbours of i
!
            ni = ni + 1
          enddo niloopnc
!
          if (lgrad2) then
!
!  For second derivatives the product of first derivatives of different neighbours are needed
!
            do ni = 1,nneigh(i)
              do nj = 1,nneigh(i)
                ind = 0
                do nn1 = 1,nneighi2
                  do nn2 = 1,nn1
                    ind = ind + 1
                    d2i(ind) = d2i(ind) + d2edZ2*d1ZiCN(nn1,ni)*d1ZiCN(nn2,nj)
                  enddo
                enddo
              enddo
            enddo
          endif
!
!  Add derivatives due to neighbours of i
!
          if (lilocal) then
            call d1add(i,maxneigh,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nboatomRptr,d1i,.false.,.true.)
          endif
          if (lgrad2) then
            if (ldophonon) then
              call d2addpdm(i,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nboatomRptr,d1i,d2i,xkv,ykv,zkv,.false.)
            else
              call d2adddm(i,maxneigh,nneigh,neighno,xneigh,yneigh,zneigh,nfreeatom,nboatomRptr,d1i,d2i,.false.)
            endif
          endif
!
!  End test over whether derivatives are needed
!
        endif
!
!  End test over whether a coordination potential was found
!
      endif
    enddo
!
!  Free memory that was required for second derivative term
!
    if (lgrad2) then
      deallocate(d1ZiCN,stat=status)
      if (status/=0) call deallocate_error('bondorderd','d1ZiCN')
    endif
  endif
!
!  Free local memory
!
  deallocate(d2BtotiR,stat=status)
  if (status/=0) call deallocate_error('bondorderd','d2BtotiR')
  deallocate(d2BtotiA,stat=status)
  if (status/=0) call deallocate_error('bondorderd','d2BtotiA')
  deallocate(d2i,stat=status)
  if (status/=0) call deallocate_error('bondorderd','d2i')
  deallocate(d1BtotiR,stat=status)
  if (status/=0) call deallocate_error('bondorderd','d1BtotiR')
  deallocate(d1BtotiA,stat=status)
  if (status/=0) call deallocate_error('bondorderd','d1BtotiA')
  deallocate(d1i,stat=status)
  if (status/=0) call deallocate_error('bondorderd','d1i')
  deallocate(zneigh,stat=status)
  if (status/=0) call deallocate_error('bondorderd','zneigh')
  deallocate(yneigh,stat=status)
  if (status/=0) call deallocate_error('bondorderd','yneigh')
  deallocate(xneigh,stat=status)
  if (status/=0) call deallocate_error('bondorderd','xneigh')
  deallocate(rneigh,stat=status)
  if (status/=0) call deallocate_error('bondorderd','rneigh')
  deallocate(ineigh,stat=status)
  if (status/=0) call deallocate_error('bondorderd','ineigh')
  deallocate(nbopotptr,stat=status)
  if (status/=0) call deallocate_error('bondorderd','nbopotptr')
  deallocate(neighno,stat=status)
  if (status/=0) call deallocate_error('bondorderd','neighno')
  deallocate(Zcn,stat=status)
  if (status/=0) call deallocate_error('bondorderd','Zcn')
  deallocate(rBOcutmax,stat=status)
  if (status/=0) call deallocate_error('bondorderd','rBOcutmax')
  deallocate(nneigh,stat=status)
  if (status/=0) call deallocate_error('bondorderd','nneigh')
  deallocate(nfreeatom,stat=status)
  if (status/=0) call deallocate_error('bondorderd','nfreeatom')
  deallocate(lopanyneigh,stat=status)
  if (status/=0) call deallocate_error('bondorderd','lopanyneigh')
  deallocate(latomdone,stat=status)
  if (status/=0) call deallocate_error('bondorderd','latomdone')
  deallocate(nboatomRptr,stat=status)
  if (status/=0) call deallocate_error('bondorderd','nboatomRptr')
!
  t2 = g_cpu_time()
  tbondorder = tbondorder + t2 - t1
!
  return
  end