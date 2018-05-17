  subroutine eemdm(lmain,lgrad1,lgrad2)
!
!  Subroutine for performing electronegativity equilisation calcns
!  according to work of Mortier.
!
!  Full distributed memory parallel version. Symmetry not supported.
!  This version should only be called when second derivatives are
!  needed and so iterative algorithm is not necessary.
!
!  lmain = .true.  => call from gulp, one off calc with output
!        = .false. => call from energy, needed for charges
!
!   1/17 Created from eem
!   2/17 Parallelisation implemented and debugged
!   7/17 Correction to parallel calculation of charges made
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
  logical, intent(in)                          :: lgrad1
  logical, intent(in)                          :: lgrad2
  logical, intent(in)                          :: lmain
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: ia
  integer(i4)                                  :: ifail
  integer(i4)                                  :: ii
  integer(i4)                                  :: ilhs
  integer(i4)                                  :: iloc
  integer(i4)                                  :: j
  integer(i4)                                  :: n
  integer(i4),                            save :: ncfold = 0
  integer(i4)                                  :: neemfoc
  integer(i4)                                  :: ni
  integer(i4)                                  :: niter
  integer(i4)                                  :: nitereem
  integer(i4)                                  :: nodeeem
  integer(i4)                                  :: nr
  integer(i4)                                  :: nregioni
  integer(i4)                                  :: status
  logical,     dimension(:), allocatable       :: leemfoc
  logical                                      :: lconverged
  logical                                      :: lDoChargeDeriv
  logical                                      :: literate
  real(dp)                                     :: chii
  real(dp)                                     :: chiloc
  real(dp)                                     :: g_cpu_time
  real(dp)                                     :: enega
  real(dp)                                     :: eself_before
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
  real(dp),    dimension(:), allocatable       :: oldqf
  real(dp),    dimension(:), allocatable       :: vfield
  real(dp),    dimension(:), allocatable       :: z
  real(dp),    dimension(:), allocatable       :: z2
!
  time1 = g_cpu_time()
!
!  Allocate local memory
!
  allocate(leemfoc(numat),stat=status)
  if (status/=0) call outofmemory('eemdm','leemfoc')
!
  lDoChargeDeriv = (index(keyword,'dcha').ne.0)
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
  neemloc = 0
!
!  Check elements
!
  do i = 1,numat
    ia = nat(i)
    if (lelementOK(ia).and.nregionno(nsft+nrelat(i)).eq.1) then
      neem = neem + 1
      neemptr(neem) = i
      neemrptr(i) = neem
      qsum = qsum - qf(i)*occuf(i)
    elseif (ia.gt.maxele) then
      call outerror('cannot use EEM with shells present',0_i4)
      call stopnow('eemdm')
    else
      qtot = qtot + qf(i)*occuf(i)
    endif
  enddo
!
  do iloc = 1,natomsonnode
    i = node2atom(iloc)
    ia = nat(i)
    if (lelementOK(ia).and.nregionno(nsft+nrelat(i)).eq.1) then
      neemloc = neemloc + 1
      neemlocptr(neemloc) = i
      neemlocrptr(i) = neemloc
    endif
  enddo
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
!  Check the memory for the square arrays
!
  if (natomsonnode+1.gt.maxd2u) then
    maxd2u = natomsonnode + 1
    call changemaxd2
  endif
  if (numat+1.gt.maxd2) then
    maxd2 = numat + 1
    call changemaxd2
  endif
!
!  Find node that has electronegativity term
!
  i = numat / (nblocksize*nprocs)
  nodeeem = numat - i*nblocksize*nprocs
!
!  Set the pointer to where the electronegativity should be as well
!
  neemptr(neemfoc+1) = numat + 1
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
  if (status/=0) call outofmemory('eemdm','z')
  allocate(z2(numat),stat=status)
  if (status/=0) call outofmemory('eemdm','z2')
!
!  If iterative then set up pointers to local elements
!
  if (literate) then
    allocate(oldqf(numat),stat=status)
    if (status/=0) call outofmemory('eemdm','oldqf')
  endif
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    allocate(vfield(numat),stat=status)
    if (status/=0) call outofmemory('eemdm','vfield')
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
        chii = qeqchi(nat(ii))
        rmui = qeqmu(nat(ii))
        q0i  = qeqq0(nat(ii))
      elseif (lSandM) then   
        chii = smchi(nat(ii))
        rmui = smmu(nat(ii))
        q0i  = smq0(nat(ii))
      elseif (lpacha) then   
        chii = chi_pacha(nat(ii))
        rmui = mu_pacha(nat(ii))
        q0i  = q0_pacha(nat(ii))
      else
        chii = chi(nat(ii))
        rmui = rmu(nat(ii))
        q0i  = q0(nat(ii))
      endif
      qf(ii) = q0i - chii/rmui
      qguesstot = qguesstot + qf(ii)*occuf(ii)
      rnguess = rnguess + occuf(ii)
    enddo
    qguesstot = (qguesstot + qtot)/rnguess
    do i = 1,neem
      ii = neemptr(i)
      qf(ii) = qf(ii) - qguesstot
      if (abs(qtot).lt.1.0d-12) qf(ii) = 1.5_dp*qf(ii)
    enddo
!
!  Transfer to qa to ensure right values are set
!
    do i = 1,nasym
      qa(i) = qf(nrel2(i))
    enddo
!
    if (index(keyword,'debu').ne.0.and.ioproc) then
      write(ioout,'(/,''  Initial guess for 1-D variable charges :'',/)')
      write(ioout,'('' Atom        Q'')')
      do i = 1,neem
        ii = neemptr(i)
        write(ioout,'(i5,1x,f12.6)') ii,qf(ii)
      enddo
      write(ioout,'(/)')
    endif
  endif
!
!  Setup coordinates
!
  do i = 1,numat
    xalat(i) = xclat(i)
    yalat(i) = yclat(i)
    zalat(i) = zclat(i)
  enddo
!
!  Store charges for convergence check
!
  if (literate) then
    do i = 1,numat
      oldqf(i) = qf(i)
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
      derv2(1:numat,1:natomsonnode) = 0.0_dp
    else
      call genpotdm(derv2,maxd2,z,1_i4)
    endif
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
!  Reduce to neem x neemloc form
!
    do i = 1,neemloc
!
!  Zero storage vector for derv2 array
!
      do j = 1,numat
        z2(j) = 0.0_dp
      enddo
!
!  Place i-j potential terms into derv2
!
      do j = 1,numat
        if (neemrptr(j).gt.0) then
!
!  Variable charge atom
!
          derv2(neemrptr(j),i) = derv2(j,atom2local(neemlocptr(i)))*occuf(j)
        else
!
!  Fixed charge atom
!
          z(i) = z(i) - qf(j)*derv2(j,atom2local(neemlocptr(i)))*occuf(j)
        endif
      enddo
    enddo
!********************************
!  Form matrix of coefficients  *
!********************************
    if (lqeq) then
      do i = 1,neemloc
        ii = neemlocptr(i)
        ilhs = neemrptr(ii)
        if (nat(ii).ne.1) then
          derv2(ilhs,i) = derv2(ilhs,i) + 2.0_dp*qeqmu(nat(ii))*occuf(ii)
        else
!
!  For hydrogen charge dependant factor must be introduced
!
          rjfac = 1.0_dp+(qf(ii)/zetah0)
          derv2(ilhs,i) = derv2(ilhs,i) + 2.0_dp*qeqmu(1)*occuf(ii)*rjfac
        endif
      enddo
    elseif (lSandM) then
      do i = 1,neemloc
        ii = neemlocptr(i)
        ilhs = neemrptr(ii)
        derv2(ilhs,i) = derv2(ilhs,i) + 2.0_dp*smmu(nat(ii))*occuf(ii)
      enddo
    elseif (lpacha) then
      do i = 1,neemloc
        ii = neemlocptr(i)
        ilhs = neemrptr(ii)
        derv2(ilhs,i) = derv2(ilhs,i) + 2.0_dp*mu_pacha(nat(ii))*occuf(ii)
      enddo
    else
      do i = 1,neemloc
        ii = neemlocptr(i)
        ilhs = neemrptr(ii)
        derv2(ilhs,i) = derv2(ilhs,i) + 2.0_dp*rmu(nat(ii))*occuf(ii)
      enddo
    endif
    do i = 1,neemloc
      derv2(neem+1,i) = 1.0_dp
    enddo
    do i = 1,neem
      ii = neemptr(i)
      derv2(i,neemloc+1) = occuf(ii)
    enddo
    derv2(neem+1,neemloc+1) = 0.0_dp
!
!  Add external potential
!
    do i = 1,neem
      ii = neemptr(i)
      z(i) = z(i) - extpotcfg(nsft+nrelat(ii))
    enddo
    if (lqeq) then
      do i = 1,neem
        ii = neemptr(i)
        ni = nat(ii)
        z(i) = z(i) - qeqchi(ni) + 2.0_dp*qeqmu(ni)*qeqq0(ni)
      enddo
    elseif (lSandM) then
      do i = 1,neem
        ii = neemptr(i)
        ni = nat(ii)
        z(i) = z(i) - smchi(ni) + 2.0_dp*smmu(ni)*smq0(ni)
      enddo
    elseif (lpacha) then
      do i = 1,neem
        ii = neemptr(i)
        ni = nat(ii)
        z(i) = z(i) - chi_pacha(ni) + 2.0_dp*mu_pacha(ni)*q0_pacha(ni)
      enddo
    else
      do i = 1,neem
        ii = neemptr(i)
        ni = nat(ii)
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
      call mpbarrier
      do i = 1,neem + 1
        ii = neemlocrptr(i)
        if (ii.gt.0) then
          write(ioout,'(i5,1x,10(1x,f9.5))') i,(derv2(j,ii),j=1,neem+1),z(i)
        endif
        call mpbarrier
      enddo
    endif
!******************
!  Invert matrix  *
!******************
    ifail = 0
    n = neem + 1
!************************
!  Symmetric inversion  *
!************************
    call matrix_inversion_library(n,1_i4,maxd2,nblocksize,derv2,0_i4,ifail)
!
!  Was inversion successful?
!
    if (ifail.ne.0) then
      call outerror('matrix inversion failed in EEM/QEq',0_i4)
      call stopnow('eemdm')
    endif
!  
!  Multiply inverse matrix and chi matrix to get charges
!
    qf(1:numat+1) = 0.0_dp
    do i = 1,neemloc
      ii = neemlocptr(i)
      qf(ii) = 0.0_dp
      do j = 1,neem + 1
        qf(ii) = qf(ii) + z(j)*derv2(j,i)
      enddo
    enddo
    if (nprocs.gt.1) then
      if (procid.eq.nodeeem) then
        i = neemloc + 1
        chiloc = 0.0_dp
        do j = 1,neem + 1
          chiloc = chiloc + z(j)*derv2(j,i)
        enddo
        qf(numat+1) = chiloc
      endif
!
!  Global sum of charges
!
      call sumall(qf,qa,numat+1_i4,"eemdm","qf")
      qf(1:numat+1) = qa(1:numat+1)
    endif
    do i = 1,neem
      ii = neemptr(i)
      qa(ii) = qf(ii)
    enddo
    enega = - qf(numat+1)
    if (literate) then
!
!  Check for convergence
!
      qdiff = 0.0_dp
      do i = 1,numat
        qd = qf(i) - oldqf(i)
        qdiff = qdiff + abs(qd)
      enddo
      qdiff = qdiff/dble(numat)
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
          qd = qf(ii) - oldqf(ii)
          qf(ii) = qf(ii) - 0.25_dp*qd
          oldqf(ii) = qf(ii)
        enddo
      endif
    endif
!
!  Transfer charges to qa
!
    do i = 1,nasym
      nr = nrel2(i)
      qa(i) = qf(nr)
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
    qi = qf(ii)
    ni = nat(ii)
    reqv = occuf(ii)
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
    eself = eself + qi*reqv*extpotcfg(nsft+nrelat(ii))
!
!  Only add on one node to avoid duplication in parallel
!
    if (ioproc) then
      nregioni = nregionno(nsft+nrelat(ii))
      eregion2region(nregioni,nregioni) = eregion2region(nregioni,nregioni) + eself - eself_before
      siteenergy(i) = siteenergy(i) + eself - eself_before 
    endif
  enddo
!*********************************
!  Calculate charge derivatives  *
!*********************************
  if (lgrad2.or.(lgrad1.and.(lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0))).or.lDoChargeDeriv) then
    call dcharged(lmain,.false.,.true.)
  endif
!
!  For Pacha we can also compute approximate NMR shifts for relevant nuclei
!
  if (lpacha) then
    allocate(qnmr(numat),stat=status)
    if (status/=0) call outofmemory('eemdm','qnmr')
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
    if (status/=0) call deallocate_error('eemdm','qnmr')
  endif
  if (lfieldcfg(ncf).or.(ntdfieldcfg(ncf).gt.0)) then
    deallocate(vfield,stat=status)
    if (status/=0) call deallocate_error('eemdm','vfield')
  endif
  if (literate) then
    deallocate(oldqf,stat=status)
    if (status/=0) call deallocate_error('eemdm','oldqf')
  endif
  deallocate(z2,stat=status)
  if (status/=0) call deallocate_error('eemdm','z2')
  deallocate(z,stat=status)
  if (status/=0) call deallocate_error('eemdm','z')
  deallocate(leemfoc,stat=status)
  if (status/=0) call deallocate_error('eemdm','leemfoc')
!
!  Timing
!
  time2 = g_cpu_time()
  teem = teem + time2 - time1
!
  return
  end
