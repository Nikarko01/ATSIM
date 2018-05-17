  subroutine eemd(lmain)
!
!  Subroutine for performing electronegativity equilisation calcns
!  according to work of Mortier.
!
!  Distributed memory version for parallel case - uses iterative
!  algorithm.
!
!  lmain = .true.  => call from gulp, one off calc with output
!        = .false. => call from energy, needed for charges
!
!   2/14 Created from eem
!   3/14 Unused lgrad1 argument removed
!   4/14 Fixed parameters for iterative charge solve replaced by input values
!   3/15 lnoqeem added
!   7/17 External potential added
!   7/17 Trap added for nloc = 0
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
!  Julian Gale, CIC, Curtin University, July 2017
!
  use control
  use configurations
  use g_constants,     only : angstoev
  use current
  use derivatives
  use element
  use energies
  use field,           only : lfieldcfg, ntdfieldcfg
  use iochannels
  use parallel
  use partial
  use symmetry
  use times
  implicit none
!
!  Passed variables
!
  logical, intent(in)                          :: lmain
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ia
  integer(i4)                                  :: ii
  integer(i4)                                  :: il
  integer(i4)                                  :: iopr
  integer(i4)                                  :: iter
  integer(i4)                                  :: j
  integer(i4)                                  :: jj
  integer(i4)                                  :: k
  integer(i4)                                  :: kk
  integer(i4)                                  :: maxp
  integer(i4),                            save :: ncfold = 0
  integer(i4)                                  :: neemfoc
  integer(i4)                                  :: ni
  integer(i4)                                  :: niter
  integer(i4)                                  :: nitereem
  integer(i4)                                  :: nloc
  integer(i4), dimension(:), allocatable       :: nlocptr
  integer(i4)                                  :: nr
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: status
  logical,     dimension(:), allocatable       :: leemfoc
  logical                                      :: lconverged
  logical                                      :: literate
  real(dp)                                     :: chii
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: enega
  real(dp)                                     :: err
  real(dp)                                     :: eself_before
  real(dp),    dimension(:,:), allocatable     :: potl
  real(dp)                                     :: q0i
  real(dp)                                     :: qdiff
  real(dp)                                     :: qd
  real(dp)                                     :: qguesstot
  real(dp)                                     :: qi
  real(dp),    dimension(:), allocatable       :: qnmr
  real(dp)                                     :: qsum
  real(dp)                                     :: qtot
  real(dp)                                     :: reqv
  real(dp)                                     :: rjfac
  real(dp)                                     :: rmui
  real(dp)                                     :: rnguess
  real(dp)                                     :: time1
  real(dp)                                     :: time2
  real(dp)                                     :: zetah0
  real(dp),    dimension(:), allocatable       :: oldqa
  real(dp),    dimension(:), allocatable       :: vfield
  real(dp),    dimension(:), allocatable       :: z
  real(dp),    dimension(:), allocatable       :: z2
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(leemfoc(numat),stat=status)
  if (status/=0) call outofmemory('eemd','leemfoc')
!
!  Assign parameters
!  Note that there are parameter sets for hydrogen:
!     nat = 1 => H+
!     nat = 2 => H-
!
  if (.not.lqeq.and..not.lSandM.and.index(keyword,'oldeem').ne.0) then
    chi(14) = 3.478_dp
    rmu(14) = 6.408_dp
    rmu(8) = 9.466_dp
    chi(1) = 3.398_dp
    rmu(1) = 20.818_dp
    chi(2) = 4.706_dp
    rmu(2) = 8.899_dp
  endif
!
!  Set up chi/mu according to method
!
  if (lqeq) then
    do i = 1,maxele
      if (abs(qeqchi(i)).gt.1.0d-6.or.abs(qeqmu(i)).gt.1.0d-6) then
        lelementOK(i) = .true.
      else
        lelementOK(i) = .false.
      endif
    enddo
  elseif (lSandM) then
    do i = 1,maxele
      if (abs(smchi(i)).gt.1.0d-6.or.abs(smmu(i)).gt.1.0d-6) then
        lelementOK(i) = .true.
      else
        lelementOK(i) = .false.
      endif
    enddo
  elseif (lpacha) then
    do i = 1,maxele
      if (abs(chi_pacha(i)).gt.1.0d-6.or.abs(rad_pacha(i)).gt.1.0d-6) then
        lelementOK(i) = .true.
        mu_pacha(i) = 0.5_dp*angstoev/rad_pacha(i)
      else
        lelementOK(i) = .false.
      endif
    enddo
  else
    do i = 1,maxele
      if (abs(chi(i)).gt.1.0d-6.or.abs(rmu(i)).gt.1.0d-6) then
        lelementOK(i) = .true.
      else
        lelementOK(i) = .false.
      endif
    enddo
  endif
!
  qsum = 0.0_dp
  qtot = 0.0_dp
  neem = 0
!
!  Check elements
!
  if (lsymopt) then
    do i = 1,nasym
      ia = iatn(i)
      if (lelementOK(ia).and.nregionno(nsft+i).eq.1) then
        neem = neem + 1
        neemptr(neem) = i
        qsum = qsum - neqv(i)*qa(i)*occua(i)
      elseif (ia.gt.maxele) then
        call outerror('cannot use EEM with shells present',0_i4)
        call stopnow('eemd')
      else
        qtot = qtot + neqv(i)*qa(i)*occua(i)
      endif
    enddo
  else
    do i = 1,numat
      ia = nat(i)
      if (lelementOK(ia).and.nregionno(nsft+nrelat(i)).eq.1) then
        neem = neem + 1
        neemptr(neem) = i
        qsum = qsum - qf(i)*occuf(i)
      elseif (ia.gt.maxele) then
        call outerror('cannot use EEM with shells present',0_i4)
        call stopnow('eemd')
      else
        qtot = qtot + qf(i)*occuf(i)
      endif
    enddo
  endif
!
!  Now find the number of fully occupied sites for EEM/QEq
!
  leemfoc(1:numat) = .false.
  do i = 1,neem
    ii = iocptr(neemptr(i))
    leemfoc(ii) = .true.
  enddo
  neemfoc = 0
  do i = 1,ncfoc
    if (leemfoc(i)) neemfoc = neemfoc + 1
  enddo
!
!  Check the memory for the linear arrays
!
  if (numat+1.gt.maxat) then
    maxat = numat + 1
    call changemaxat
  endif
!
!  Set the pointer to where the electronegativity should be as well
!
  neemptr(neemfoc+1) = nasym + 1
!
!  Decide on total EEM fragment charge - for cluster
!  there is no constraint so use sum of initial charges.
!  For periodic system, charge must be equal to the
!  negative sum of the non-EEM ion charges.
!
  if (ndim.eq.0) then
    qtot = qsum
  endif
!*****************************************************************
!  Is hydrogen present in QEq? If so then solution is iterative  *
!*****************************************************************
  literate = .false.
  if (lqeq) then
    i = 0
    do while (i.lt.nasym.and..not.literate)
      i = i + 1
      literate = (iatn(i).eq.1)
    enddo
  endif
  if (literate) then
    nitereem = nqeqitermax
    zetah0 = 0.529177_dp*0.75_dp/qeqrad(1)
  else
    nitereem = 1
  endif
!
!  Allocate local memory that depends on neem
!
  allocate(z(max(neem+1,numat)),stat=status)
  if (status/=0) call outofmemory('eemd','z')
  allocate(z2(numat),stat=status)
  if (status/=0) call outofmemory('eemd','z2')
  allocate(nlocptr(neem+1),stat=status)
  if (status/=0) call outofmemory('eemd','nlocptr')
!
!  Set up pointers to local elements for parallelisation
!
  if (nprocs.eq.1) then
    nloc = neem
    do i = 1,neem
      nlocptr(i) = i
    enddo
  else
    nloc = 0
    do i = procid+1,neem,nprocs
      nloc = nloc + 1
      nlocptr(nloc) = i
    enddo
  endif
!
!  Allocate the memory for the potential array
!
!  procid = 0 carries the extra row for the constraints
!
  maxp = numat + 1
  if (procid.eq.0) then
    allocate(potl(maxp,nloc+1),stat=status)
    if (status/=0) call outofmemory('eemd','potl')
  else
    allocate(potl(maxp,nloc),stat=status)
    if (status/=0) call outofmemory('eemd','potl')
  endif
!
  if (literate) then
    allocate(oldqa(nasym),stat=status)
    if (status/=0) call outofmemory('eemd','oldqa')
  endif
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    allocate(vfield(numat),stat=status)
    if (status/=0) call outofmemory('eemd','vfield')
  endif
!
!  Calculate reciprocal lattice vectors
!
  if (lewald.and.ndim.gt.1) call kindex
!
!  For 1-D case set guess at charges based EEM parameters
!  and scaled up by 1.5 to allow for increase in ionicity.
!
  if (ndim.eq.1.and.ncf.ne.ncfold) then
    ncfold = ncf
    qguesstot = 0.0_dp
    rnguess = 0.0_dp
    do i = 1,neem
      ii = neemptr(i)
      if (lqeq) then   
        chii = qeqchi(iatn(ii))
        rmui = qeqmu(iatn(ii))
        q0i  = qeqq0(iatn(ii))
      elseif (lSandM) then   
        chii = smchi(iatn(ii))
        rmui = smmu(iatn(ii))
        q0i  = smq0(iatn(ii))
      elseif (lpacha) then   
        chii = chi_pacha(iatn(ii))
        rmui = mu_pacha(iatn(ii))
        q0i  = q0_pacha(iatn(ii))
      else
        chii = chi(iatn(ii))
        rmui = rmu(iatn(ii))
        q0i  = q0(iatn(ii))
      endif
      qa(ii) = q0i - chii/rmui
      qguesstot = qguesstot + qa(ii)*dble(neqv(ii))*occua(ii)
      rnguess = rnguess + dble(neqv(ii))*occua(ii)
    enddo
    qguesstot = (qguesstot + qtot)/rnguess
    do i = 1,neem
      ii = neemptr(i)
      qa(ii) = qa(ii) - qguesstot
      if (abs(qtot).lt.1.0d-12) qa(ii) = 1.5_dp*qa(ii)
      do j = 1,numat
        if (nrelat(j).eq.ii) then
          qf(j) = qa(ii)
        endif
      enddo
    enddo
    if (index(keyword,'debu').ne.0) then
      if (ioproc) then
        write(ioout,'(/,''  Initial guess for 1-D variable charges :'',/)')
        write(ioout,'('' Atom        Q'')')
      endif
      iopr = 0
      do i = 1,neem
        if (iopr.eq.procid) then
          ii = neemptr(i)
          write(ioout,'(i5,1x,f12.6)') ii,qa(ii)
        endif
        iopr = iopr + 1
        iopr = mod(iopr,nprocs)
      enddo
      write(ioout,'(/)')
    endif
  endif
!
!  Setup coordinates
!
  if (lsymopt) then
    do i = 1,nasym
      nr = nrel2(i)
      xalat(i) = xclat(nr)
      yalat(i) = yclat(nr)
      zalat(i) = zclat(nr)
    enddo
  else
    do i = 1,numat
      xalat(i) = xclat(i)
      yalat(i) = yclat(i)
      zalat(i) = zclat(i)
    enddo
  endif
!
!  Store charges for convergence check
!
  if (literate) then
    do i = 1,nasym
      oldqa(i) = qa(i)
    enddo
  endif
!
!  Generate electric field potential
!
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    vfield(1:numat) = 0.0_dp
    call electricfieldpotl(vfield)
  endif
!****************************
!  Start of iterative loop  *
!****************************
  lconverged = .false.
  niter = 0
  if (literate.and.lmain.and.ioproc) then
    write(ioout,'(''  Iterative solution of QEq :'',/)')
  endif
  do while (niter.lt.nitereem.and..not.lconverged)
    niter = niter + 1
!
!  Zero right hand vector
!
    z(1:neem) = 0.0_dp
!*******************************************
!  Generate electrostatic potential terms  *
!*******************************************
    if (lnoqeem) then
      potl(1:numat,1:nloc) = 0.0_dp
    else
      if (lsymopt) then
        call sympotd(nloc,nlocptr,potl,maxp,z)
      else
        call genpotd(nloc,nlocptr,potl,maxp,z)
      endif
    endif
!
!  Globalise z
!
    call sumall(z,z2,neem,"eemd","z")
    z(1:neem) = z2(1:neem)
!
!  From S & M, where z has been set without reference to neem, reduce elements to those that are needed
!
    if (lSandM) then
      do i = 1,neem
        z2(i) = z(neemptr(i))
      enddo
      z(1:neem) = z2(1:neem)
    endif
!
!  Reduce to nasym x nasym form
!
    do il = 1,nloc
      i = nlocptr(il)
!
!  Zero storage vector for potl array
!
      do j = 1,numat
        z2(j) = 0.0_dp
      enddo
!
!  Place i-j potential terms into potl
!
      do j = 1,numat
        k = nrelat(j)
        jj = 1
        kk = neemptr(jj)
        do while (jj.lt.neem.and.kk.ne.k)
          jj = jj + 1
          kk = neemptr(jj)
        enddo
!
!  Variable j charge case
!
        if (kk.eq.k) then
          z2(jj) = z2(jj) + potl(j,il)*occuf(j)
        else
          z(i) = z(i) - qf(j)*potl(j,il)*occuf(j)
        endif
      enddo
!
!  Copy temporary storage vector back into potl array
!
      do j = 1,neem
        potl(j,il) = z2(j)
      enddo
    enddo
!********************************
!  Form matrix of coefficients  *
!********************************
    if (nloc.gt.0) then
      if (lqeq) then
        do il = 1,nloc
          i = nlocptr(il)
          ii = neemptr(i)
          if (iatn(ii).ne.1) then
            potl(i,il) = potl(i,il) + 2.0_dp*qeqmu(iatn(ii))*occua(ii)
          else
!
!  For hydrogen charge dependant factor must be introduced
!
            rjfac = 1.0_dp + (qa(ii)/zetah0)
            potl(i,il) = potl(i,il) + 2.0_dp*qeqmu(1)*occua(ii)*rjfac
          endif
        enddo
      elseif (lSandM) then
        do il = 1,nloc
          i = nlocptr(il)
          ii = neemptr(i)
          potl(i,il) = potl(i,il) + 2.0_dp*smmu(iatn(ii))*occua(ii)
        enddo
      elseif (lpacha) then
        do il = 1,nloc
          i = nlocptr(il)
          ii = neemptr(i)
          potl(i,il) = potl(i,il) + 2.0_dp*mu_pacha(iatn(ii))*occua(ii)
        enddo
      else
        do il = 1,nloc
          i = nlocptr(il)
          ii = neemptr(i)
          potl(i,il) = potl(i,il) + 2.0_dp*rmu(iatn(ii))*occua(ii)
        enddo
      endif
      potl(nasym+1,nloc) = 0.0_dp
      do il = 1,nloc
        i = nlocptr(il)
        potl(neem+1,il) = 1.0_dp
      enddo
      if (procid.eq.0) then
        do i = 1,neem
          ii = neemptr(i)
          potl(i,nloc+1) = dble(neqv(ii))*occua(ii)
        enddo
        potl(neem+1,nloc+1) = 0.0_dp
      endif
    endif
!
!  Add external potential
!
    do i = 1,neem
      ii = neemptr(i)
      z(i) = z(i) - extpotcfg(nsft+ii)
    enddo
    if (lqeq) then
      do i = 1,neem
        ii = neemptr(i)
        ni = iatn(ii)
        z(i) = z(i) - qeqchi(ni) + 2.0_dp*qeqmu(ni)*qeqq0(ni)
      enddo
    elseif (lSandM) then
      do i = 1,neem
        ii = neemptr(i)
        ni = iatn(ii)
        z(i) = z(i) - smchi(ni) + 2.0_dp*smmu(ni)*smq0(ni)
      enddo
    elseif (lpacha) then
      do i = 1,neem
        ii = neemptr(i)
        ni = iatn(ii)
        z(i) = z(i) - chi_pacha(ni) + 2.0_dp*mu_pacha(ni)*q0_pacha(ni)
      enddo
    else
      do i = 1,neem
        ii = neemptr(i)
        ni = iatn(ii)
        z(i) = z(i) - chi(ni) + 2.0_dp*rmu(ni)*q0(ni)
      enddo
    endif
    if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
      do i = 1,neem
        ii = neemptr(i)
        z(i) = z(i) - vfield(ii)
      enddo
    endif
    z(neem+1) = - qtot
!
!  Debugging output
!
    if (index(keyword,'debu').ne.0) then
      if (ioproc) then
        write(ioout,'(/,''  EEM/QEq Matrix :'',/)')
      endif
      iopr = 0
      il = 0
      do i = 1,neem
        if (procid.eq.iopr) then
          il = il + 1
          write(ioout,'(i5,10(1x,f9.5))') i,(potl(j,il),j=1,neem+1),z(i)
        endif
        iopr = iopr + 1
        iopr = mod(iopr,nprocs)
      enddo
      if (procid.eq.0) then
        write(ioout,'(i5,10(1x,f9.5))') i,(potl(j,nloc+1),j=1,neem+1),z(neem+1)
      endif
    endif
!***************************
!  Iterative charge solve  *
!***************************
    do i = 1,neem
      ii = neemptr(i)
      qf(i) = qa(ii)
    enddo
    qf(neem+1) = 0.0_dp
!
!  Solve using iterative route
!
    call dbcgsolve(neem+1_i4,nloc,nlocptr,potl,maxp,z,qf,qitertol,nqitermax,iter,err)
!
    if (ioproc) then
      if (index(keyword,'verb').ne.0) then
        write(ioout,'('' Number of iterations / error in dbcgsolve = '',i4,1x,f16.14)') iter,err
      endif
    endif
!
!  Was iterative solution successful?
!
    if (iter.ge.nqitermax) then
      call outerror('iterative charge solution failed in eemd',0_i4)
      call stopnow('eemd')
    endif
    do i = 1,neem
      ii = neemptr(i)
      qa(ii) = qf(i)
    enddo
    enega = - qf(neem+1)
!
    if (literate) then
!
!  Check for convergence
!
      qdiff = 0.0_dp
      do i = 1,nasym
        qd = qa(i) - oldqa(i)
        qdiff = qdiff + abs(qd)
      enddo
      qdiff = qdiff/dble(nasym)
      lconverged = (qdiff.lt.qeqscfcrit)
      if (lmain.and.ioproc) then
        write(ioout,'(''  ** Cycle : '',i4,'' Qdiff : '',f10.8)') niter,qdiff
      endif
      if (.not.lconverged) then
!
!  Damp change to improve convergence
!
        do i = 1,neem
          ii = neemptr(i)
          qd = qa(ii) - oldqa(ii)
          qa(ii) = qa(ii) - 0.25_dp*qd
          oldqa(ii) = qa(ii)
        enddo
      endif
    endif
!
!  Transfer charges to qf
!
    do i = 1,numat
      nr = nrelat(i)
      qf(i) = qa(nr)
    enddo
!*****************************
!  End loop over iterations  *
!*****************************
  enddo
!
!  Store charges in configurational array
!
  do i = 1,nasym
    qlcfg(nsft+i) = qa(i)
  enddo
!**************************
!  Calculate self energy  *
!**************************
  eself = 0.0_dp
  do i = 1,neem
    ii = neemptr(i)
    qi = qa(ii)
    ni = iatn(ii)
    reqv = dble(neqv(ii))*occua(ii)
    eself_before = eself
    if (lqeq) then
      q0i = qeqq0(ni)
      if (ni.ne.1) then
        eself = eself + (qi-q0i)*reqv*(qeqchi(ni)+(qi-q0i)*qeqmu(ni))
      else
        eself = eself + (qi-q0i)*reqv*(qeqchi(ni)+(qi-q0i)*qeqmu(ni)*(1.0_dp+(2.0_dp*(qi-q0i)/(3.0_dp*zetah0))))
      endif
    elseif (lSandM) then
      q0i = smq0(ni)
      eself = eself + (qi-q0i)*reqv*(smchi(ni)+(qi-q0i)*smmu(ni))
    elseif (lpacha) then
      q0i = q0_pacha(ni)
      eself = eself + (qi-q0i)*reqv*(chi_pacha(ni)+(qi-q0i)*mu_pacha(ni))
    else
      q0i = q0(ni)
      eself = eself + (qi-q0i)*reqv*(chi(ni)+(qi-q0i)*rmu(ni))
    endif
!
!  Add external potential for site
!
    eself = eself + qi*reqv*extpotcfg(nsft+ii)
!
    nregioni = nregionno(nsft+ii)
    eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + eself - eself_before
!
    siteenergy(i) = siteenergy(i) + eself - eself_before 
  enddo
!
!  For Pacha we can also compute approximate NMR shifts for relevant nuclei
!
  if (lpacha) then
    allocate(qnmr(numat),stat=status)
    if (status/=0) call outofmemory('eemd','qnmr')
    call getnmr(nasym,iatn,qa,qnmr)
  endif
!*******************
!  Output results  *
!*******************
  if ((lmain.or.index(keyword,'debu').ne.0).and.ioproc) then
    if (lqeq) then
      write(ioout,'(//,''  Final charges from QEq :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,iatn(i),qa(i)
      enddo
    elseif (lpacha) then
      write(ioout,'(//,''  Final charges from PACHA-EEM :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge       Chemical shift'')')
      write(ioout,'(''                                                 (e)            (ppm)     '')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7,2x,f11.4)') i,iatn(i),qa(i),qnmr(i)
      enddo
    else
      write(ioout,'(//,''  Final charges from EEM :'',/)')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      write(ioout,'(''    Atom no.            Atomic No.             Charge'')')
      write(ioout,'(''--------------------------------------------------------------------------------'')')
      do i = 1,nasym
        write(ioout,'(4x,i6,17x,i3,16x,f10.7)') i,iatn(i),qa(i)
      enddo
    endif
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    write(ioout,'(''  Electronegativity = '',f16.6,'' eV'')') enega
    write(ioout,'(''  Self energy       = '',f16.6,'' eV'')') eself
    write(ioout,'(''--------------------------------------------------------------------------------'')')
    if (lqeq) then
      if (literate) then
        if (lconverged) then
          write(ioout,'(/,''  Charges converged in '',i3,'' iterations'',/)') niter
        else
          write(ioout,'(/,''  Failed to converged after '',i3,'' iterations'',/)') nitereem
        endif
      else
        write(ioout,'(/,''  No hydrogens present - no iteration needed'',/)')
      endif
    endif
  endif
!
!  Free local memory 
!
  if (lpacha) then
    deallocate(qnmr,stat=status)
    if (status/=0) call deallocate_error('eemd','qnmr')
  endif
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    deallocate(vfield,stat=status)
    if (status/=0) call deallocate_error('eemd','vfield')
  endif
  if (literate) then
    deallocate(oldqa,stat=status)
    if (status/=0) call deallocate_error('eemd','oldqa')
  endif
  deallocate(potl,stat=status)
  if (status/=0) call deallocate_error('eemd','potl')
  deallocate(nlocptr,stat=status)
  if (status/=0) call deallocate_error('eemd','nlocptr')
  deallocate(z2,stat=status)
  if (status/=0) call deallocate_error('eemd','z2')
  deallocate(z,stat=status)
  if (status/=0) call deallocate_error('eemd','z')
  deallocate(leemfoc,stat=status)
  if (status/=0) call deallocate_error('eemd','leemfoc')
!
!  Timing
!
  time2 = g_cpu_time()
  teem = teem + time2 - time1
!
  return
  end
