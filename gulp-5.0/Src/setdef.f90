  subroutine setdef(nldef,ndefst,lrestore)
!
!  Set up for defect calculation for current configuration
!
!  Scalars for defect calculation:
!
!  nreg1    = no. of atoms in defective region 1
!  nreg1old = no. of atoms in perfect region 1
!  qreg1    = charge on defective region 1
!  qreg1old = charge on perfect region 1
!  qdef     = net charge on defect for ML polarisation energy
!  rfind    = distance tolerance squared for locating vacancies
!
!  Arrays for defect calculation:
!
!  xdefe      = defective region cartesian coordinates
!  xperf      = perfect region 1 cartesian coordinates
!  qdefe      = charges for defective region 1
!  qp         = charges for perfect region 1
!  radefe     = radii for defective region 1
!  natdefe    = atomic nos for defective region 1
!  natp       = atomic nos for perfect region 1
!  ntypdefe   = atom types for defective region 1
!  ntypep     = atom types for perfect region 1
!  ndefmol    = molecule no. for defective region 1
!  ndefmolp   = molecule no. for perfect region 1
!  ndefind    = molecule cell index for defective region 1
!  ndefindp   = molecule cell index for perfect region 1
!  npsite     = perfect lattice site for perfect region 1
!  nbondsdef  = number of bonds for defective region 1
!  nbondeddef = bonding list for defective region 1
!  nreldef    = maps defect atoms to equivalent atom in
!               perfect lattice (0=> interstitial)
!               NB: This is a change from previous versions
!               where nreldef mapped to the perfect region 1
!               atom and then used in npsite look up.
!
!  Logicals:
!
!  lr1created = logical array, if true the atom has been created as
!             part of a defect - this is used to stop subsequent
!             defect commands from removing it
!  lrestore = if true then region 1 is to be set from configuration 
!             arrays rather than generated
!  lexpand  = if true then reg1last is greater than zero so exclude
!             all ions within this radius when creating region 1
!  ldefbsmat= logical for each atom in region 1 for whether it has
!             a breathing shell or not
!  ldbsm    = if .true. then BSM is present in region1
!  ldqmatom = logical for each atom in region 1 for whether it is
!             also a qunatum mechanically treated atom in an
!             embedded cluster calculation
!  ldeffix  = logical array, if .true. then the defect is to be
!             held fixed during optimisation
!  ldcellr  = if this is .true. then region 1 is constructed based
!             on unit cells rather than radius (special option for
!             Sasha!).
!
!  11/95 ldqmatom pointer added for embedded cluster calculations
!  11/95 sorting stopped if old region 1 and not being expanded
!        for benefit of embedded cluster stuff
!  11/96 fixing of interstitials and impurities added
!   7/00 lr1created now passed to setmoldef
!   3/03 Setting of nadd maded safer
!   7/05 Memory deallocation cleaned up
!   8/05 Bug in setting of ldfix corrected
!  12/07 Unused variables removed
!   5/08 Defect bonding array structure changed
!   4/10 Search for region 1 ions now uses iterative algorithm to 
!        avoid problems with choosing a suitable nadd value.
!   7/10 References to ijx wrapped so that they are only accessed
!        when lmol is true
!   3/13 Restore option modified to use data in configuration arrays
!        rather than scratch files
!   3/13 Meaning of nreldef value changed
!   3/13 ndptr replaced by ndintptr and ndvacptr
!   3/13 Sort moved to separate subroutine
!   9/15 cuts2 renamed to cut2s for consistency 
!   5/17 Forcing of lnewdefalg added here
!   5/17 Calling of changemaxr1at changed to ensure that there is 
!        nreg1 + 1 words in arrays.
!   5/17 Calculation of nreg1onnode for cores and shells added
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
!  Julian Gale, CIC, Curtin University, May 2017
!
  use configurations, only : lbsmat,radcfg
  use control
  use current
  use defects
  use element,        only : maxele
  use energies
  use general
  use molecule
  use parallel
  use reallocate
  use shells
  use species
  use sutton,         only : lsuttonc
  use symmetry
  implicit none
!
!  Passed variables
!
  integer(i4)                                  :: ndefst
  integer(i4)                                  :: nldef
  logical                                      :: lrestore
!
!  Local variables
!
  integer(i4)                                  :: i
  integer(i4)                                  :: icount
  integer(i4)                                  :: idir
  integer(i4)                                  :: ierror
  integer(i4)                                  :: ii
  integer(i4)                                  :: ijx
  integer(i4)                                  :: ijy
  integer(i4)                                  :: ijz
  integer(i4)                                  :: imid
  integer(i4)                                  :: inat
  integer(i4)                                  :: ind
  integer(i4)                                  :: indm
  integer(i4)                                  :: indmj
  integer(i4)                                  :: iptr
  integer(i4), dimension(:), allocatable       :: itmp
  integer(i4)                                  :: ityp
  integer(i4)                                  :: itype
  integer(i4)                                  :: ix
  integer(i4)                                  :: ixj
  integer(i4)                                  :: ixji
  integer(i4)                                  :: iy
  integer(i4)                                  :: iyj
  integer(i4)                                  :: iyji
  integer(i4)                                  :: iz
  integer(i4)                                  :: izj
  integer(i4)                                  :: izji
  integer(i4)                                  :: j
  integer(i4)                                  :: jdir
  integer(i4)                                  :: jj
  integer(i4)                                  :: jmid
  integer(i4)                                  :: jtyp
  integer(i4)                                  :: k
  integer(i4)                                  :: kdir
  integer(i4)                                  :: kk
  integer(i4)                                  :: kmid
  integer(i4)                                  :: mvar
  integer(i4)                                  :: nati
  integer(i4)                                  :: ncores
  integer(i4)                                  :: ndind
  integer(i4)                                  :: ndt
  integer(i4)                                  :: newind
  integer(i4)                                  :: nfix
  integer(i4)                                  :: ni
  integer(i4)                                  :: niind
  integer(i4)                                  :: ni2
  integer(i4)                                  :: nj
  integer(i4)                                  :: nmi
  integer(i4)                                  :: nmj
  integer(i4)                                  :: nmolat
  integer(i4)                                  :: node
  integer(i4)                                  :: nri
  integer(i4)                                  :: ntypi
  integer(i4)                                  :: nv
  integer(i4)                                  :: nvec
  integer(i4)                                  :: nvecj0
  integer(i4)                                  :: nveck0
  integer(i4)                                  :: nvind
  integer(i4)                                  :: status
  logical                                      :: lallfound1
  logical                                      :: lallfound2
  logical                                      :: lallfound3
  logical                                      :: lbre
  logical                                      :: lcore
  logical                                      :: lexpand
  logical                                      :: lfind
  logical                                      :: lfound
  logical                                      :: lnosort
  logical                                      :: lshel
  real(dp)                                     :: cut2s
  real(dp)                                     :: proj1
  real(dp)                                     :: proj2
  real(dp)                                     :: proj3
  real(dp)                                     :: q
  real(dp)                                     :: qreg1
  real(dp)                                     :: qreg1old
  real(dp)                                     :: r
  real(dp)                                     :: r1
  real(dp)                                     :: r12
  real(dp)                                     :: r12s
  real(dp)                                     :: r1ex
  real(dp)                                     :: r1xf
  real(dp)                                     :: r1yf
  real(dp)                                     :: r1zf
  real(dp)                                     :: r2i
  real(dp)                                     :: r2j
  real(dp)                                     :: r2k
  real(dp)                                     :: r2xf
  real(dp)                                     :: r2yf
  real(dp)                                     :: r2zf
  real(dp)                                     :: r3xf
  real(dp)                                     :: r3yf
  real(dp)                                     :: r3zf
  real(dp)                                     :: r2
  real(dp)                                     :: ra
  real(dp)                                     :: rcs
  real(dp)                                     :: rcx1
  real(dp)                                     :: rcy1
  real(dp)                                     :: rcz1
  real(dp)                                     :: rcx2
  real(dp)                                     :: rcy2
  real(dp)                                     :: rcz2
  real(dp)                                     :: rcx3
  real(dp)                                     :: rcy3
  real(dp)                                     :: rcz3
  real(dp)                                     :: rfind
  real(dp)                                     :: rmax
  real(dp)                                     :: rnm
  real(dp)                                     :: rnorm
  real(dp)                                     :: rtest
  real(dp)                                     :: rx
  real(dp)                                     :: ry
  real(dp)                                     :: rz
  real(dp)                                     :: rxi
  real(dp)                                     :: ryi
  real(dp)                                     :: rzi
  real(dp)                                     :: rxj
  real(dp)                                     :: ryj
  real(dp)                                     :: rzj
  real(dp)                                     :: rvt(3,3)
  real(dp)                                     :: xci
  real(dp)                                     :: yci
  real(dp)                                     :: zci
  real(dp)                                     :: xcrd
  real(dp)                                     :: ycrd
  real(dp)                                     :: zcrd
  real(dp)                                     :: xcs
  real(dp)                                     :: ycs
  real(dp)                                     :: zcs
  real(dp)                                     :: xd
  real(dp)                                     :: yd
  real(dp)                                     :: zd
  real(dp)                                     :: xdiff
  real(dp)                                     :: ydiff
  real(dp)                                     :: zdiff
  real(dp)                                     :: xi
  real(dp)                                     :: yi
  real(dp)                                     :: zi
  real(dp)                                     :: xij
  real(dp)                                     :: yij
  real(dp)                                     :: zij
  real(dp)                                     :: xm
  real(dp)                                     :: ym
  real(dp)                                     :: zm
  real(dp)                                     :: xv
  real(dp)                                     :: yv
  real(dp)                                     :: zv
!
  r1ex = reg1last(ncf)
  lexpand = (r1ex.gt.0.0_dp)
  r1ex = r1ex*r1ex
  lr2f = .true.
  e2b = 0.0_dp
  ndcon = 0
!
!  Set sorting flag
!
  lnosort = (lrestore.and..not.lexpand)
!
!  Set new defect algorithm flag - must use this algorithm if EAM model is being used.
!
  if (lsuttonc.or.ldcellr) lnewdefalg = .true.
!
!  Set up local variables
!
  r1 = reg1(ncf)
  cut2s = cuts*cuts
!
!  Check region 1 radius
!
  if (abs(r1).lt.1.0d-6) then
    call outerror('region 1 radius is zero',0_i4)
    call stopnow('setdef')
  endif
  r12 = r1*r1
  r12s = r1 + cuts
  r12s = r12s*r12s
  rfind = 0.04_dp
  nreg1 = 0
  nreg1old = 0
!
!  Need full cell vectors now for fractional coordinate
!  manipulations and image testing.
!
  do i = 1,3
    rvt(i,1) = rv(i,1)
    rvt(i,2) = rv(i,2)
    rvt(i,3) = rv(i,3)
  enddo
  if (ifhr(ncf).ne.1) call uncentre(rvt)
  r1xf = rvt(1,1)
  r1yf = rvt(2,1)
  r1zf = rvt(3,1)
  r2xf = rvt(1,2)
  r2yf = rvt(2,2)
  r2zf = rvt(3,2)
  r3xf = rvt(1,3)
  r3yf = rvt(2,3)
  r3zf = rvt(3,3)
!*************************
!  Locate defect centre  *
!*************************
  if (ndcentyp(ncf).eq.1.or.ndcentyp(ncf).eq.2) then
!
!  Atomic position
!
    if (ndcentyp(ncf).eq.2) then
      ni = nint(xdcent(ncf)) - nsft
    else
      ni = nint(xdcent(ncf))
    endif
    nri = nrel2(ni)
    xdc = xclat(nri)
    ydc = yclat(nri)
    zdc = zclat(nri)
  elseif (ndcentyp(ncf).eq.3) then
!
!  Fractional coordinates
!
    xi = xdcent(ncf)
    yi = ydcent(ncf)
    zi = zdcent(ncf)
    xdc = xi*r1xf + yi*r2xf + zi*r3xf
    ydc = xi*r1yf + yi*r2yf + zi*r3yf
    zdc = xi*r1zf + yi*r2zf + zi*r3zf
  elseif (ndcentyp(ncf).eq.4) then
!
!  Cartesian coordinates
!
    xdc = xdcent(ncf)
    ydc = ydcent(ncf)
    zdc = zdcent(ncf)
  else
!
!  Molecule centroid
!
    nmi = nint(xdcent(ncf))
    xdc = 0.0_dp
    ydc = 0.0_dp
    zdc = 0.0_dp
    nmolat = 0
    do i = 1,numat
      if (natmol(i).eq.nmi) then
        nmolat = nmolat + 1
        xm = xclat(i)
        ym = yclat(i)
        zm = zclat(i)
        indm = nmolind(i)
        iz = (indm/100) - 5
        ind = indm - 100*(iz + 5)
        iy = (ind/10) - 5
        ind = ind - 10*(iy + 5)
        ix = ind - 5
        xm = xm + ix*r1x + iy*r2x + iz*r3x
        ym = ym + ix*r1y + iy*r2y + iz*r3y
        zm = zm + ix*r1z + iy*r2z + iz*r3z
        xdc = xdc + xm
        ydc = ydc + ym
        zdc = zdc + zm
      endif
    enddo
    rnm = 1.0_dp/float(nmolat)
    xdc = rnm*xdc
    ydc = rnm*ydc
    zdc = rnm*zdc
  endif
!******************************
!  Generate perfect region 1  *
!******************************
  qreg1old = 0.0_dp
  qreg1 = 0.0_dp
!*********************
!  Find cores first  *
!*********************
  do j = 1,numat
    if (nat(j).le.maxele) then
      xcrd = xclat(j) - xdc
      ycrd = yclat(j) - ydc
      zcrd = zclat(j) - zdc
      if (lmol) then
        nmj = natmol(j)
        indmj = nmolind(j)
        izj = (indmj/100) - 5
        ind = indmj - 100*(izj + 5)
        iyj = (ind/10) - 5
        ind = ind - 10*(iyj + 5)
        ixj = ind - 5
      endif
!
!  Search unit cells iteratively for all images
!
      xij = xcrd
      yij = ycrd
      zij = zcrd
!
!  Find projection of cell vector 3 on to i - j vector
!
      rnorm = xij*xij + yij*yij + zij*zij
      if (rnorm.gt.1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
      proj3 = rnorm*recipc*(xij*r3x + yij*r3y + zij*r3z)
      kmid = nint(proj3)
      xij = xij - kmid*r3x
      yij = yij - kmid*r3y
      zij = zij - kmid*r3z
!
!  Find projection of cell vector 2 on to i - j vector
!
      rnorm = xij*xij + yij*yij + zij*zij
      if (rnorm.gt.1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
      proj2 = rnorm*recipb*(xij*r2x + yij*r2y + zij*r2z)
      jmid = nint(proj2)
      xij = xij - jmid*r2x
      yij = yij - jmid*r2y
      zij = zij - jmid*r2z
!
!  Find projection of cell vector 1 on to i - j vector
!
      rnorm = xij*xij + yij*yij + zij*zij
      if (rnorm.gt.1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
      proj1 = rnorm*recipa*(xij*r1x + yij*r1y + zij*r1z)
      imid = nint(proj1)
      xij = xij - imid*r1x
      yij = yij - imid*r1y
      zij = zij - imid*r1z
!
!  Adjust molecule indices for coordinate shift
!
      if (lmol) then
        ijx = ixj + imid
        ijy = iyj + jmid
        ijz = izj + kmid
      endif
!
!  Initialise number of distances to zero
!
      nvec = 0
!
!  Outer loop over first cell vector direction
!
      do idir = 1,-1,-2
!
!  Reinitialise distance squared
!
        r2i = 10000.0_dp*r12
!
!  Loop over first cell vector
!
        lallfound1 = .false.
        if (idir.eq.1) then
          ii = 0
        else
          ii = - 1
        endif
!
!  Set initial coordinate vector
!
        rxi = xij + dble(ii)*r1x
        ryi = yij + dble(ii)*r1y
        rzi = zij + dble(ii)*r1z
!
!  Set increment vector
!
        rcx1 = dble(idir)*r1x
        rcy1 = dble(idir)*r1y
        rcz1 = dble(idir)*r1z
!
        do while (.not.lallfound1)
!
!  Save number of vectors before search over second direction
!
          nvecj0 = nvec
!
!  Outer loop over second cell vector direction
!
          do jdir = 1,-1,-2
!
!  Reinitialise saved distance squared
!
            r2j = 10000.0_dp*r12
!
!  Loop over second cell vector
!
            lallfound2 = .false.
            if (jdir.eq.1) then
              jj = 0
            else
              jj = - 1
            endif
!
!  Set initial coordinate vector
!
            rxj = rxi + dble(jj)*r2x
            ryj = ryi + dble(jj)*r2y
            rzj = rzi + dble(jj)*r2z
!
!  Set increment vector
!
            rcx2 = dble(jdir)*r2x
            rcy2 = dble(jdir)*r2y
            rcz2 = dble(jdir)*r2z
!
            do while (.not.lallfound2)
!
!  Save number of vectors before search over third direction
!
              nveck0 = nvec
!
!  Outer loop over third cell vector direction
!
              do kdir = 1,-1,-2
!
!  Reinitialise saved distance squared
!
                r2k = 10000.0_dp*r12
!
!  Loop over third cell vector
!
                lallfound3 = .false.
                if (kdir.eq.1) then
                  kk = 0
                else
                  kk = - 1
                endif
!
!  Set initial coordinate vector
!
                rx = rxj + dble(kk)*r3x
                ry = ryj + dble(kk)*r3y
                rz = rzj + dble(kk)*r3z
!
!  Set increment vector
!
                rcx3 = dble(kdir)*r3x
                rcy3 = dble(kdir)*r3y
                rcz3 = dble(kdir)*r3z
!
                do while (.not.lallfound3)
!
!  Calculate square of distance
!
                  r2 = rx*rx + ry*ry + rz*rz
!
!  Check distance squared against cutoff squared
!
                  if (r2.le.r12.or.ldcellr) then
                    nvec = nvec + 1
!
!  Valid region 1 ion
!
                    nreg1old = nreg1old + 1
                    if (nreg1old.ge.maxr1at) then
                      maxr1at = nreg1old + 20
                      call changemaxr1at
                    endif
                    xperf(nreg1old) = rx + xdc
                    yperf(nreg1old) = ry + ydc
                    zperf(nreg1old) = rz + zdc
                    qp(nreg1old) = qf(j)
                    occp(nreg1old) = occuf(j)
                    natp(nreg1old) = nat(j)
                    ntypep(nreg1old) = nftype(j)
                    npsite(nreg1old) = j
                    qreg1old = qreg1old + qp(nreg1old)*occp(nreg1old)
!
!  Molecule handling
!
                    if (lmol) then
                      ixji = ii - ijx
                      iyji = jj - ijy
                      izji = kk - ijz
                      newind = (ixji + 5) + 10*(iyji + 5) + 100*(izji + 5)
                      ndefmolp(nreg1old) = nmj
                      if (nmj.eq.0) then
                        ndefindp(nreg1old) = 0
                      else
                        ndefindp(nreg1old) = newind
                      endif
                    endif
                  endif
!
!  Increment by third vector
!
                  kk = kk + kdir
                  rx = rx + rcx3
                  ry = ry + rcy3
                  rz = rz + rcz3
!
!  Check to see if this direction is complete
!
                  if (ldcellr) then
                    lallfound3 = (abs(kk).ge.1)
                  else
                    lallfound3 = (r2.gt.r2k.and.r2.gt.r12)
                  endif
                  r2k = r2
                enddo
              enddo
!
!  Increment by second vector
!
              jj = jj + jdir
              rxj = rxj + rcx2
              ryj = ryj + rcy2
              rzj = rzj + rcz2
!
!  Check to see if this direction is complete
!
              if (ldcellr) then
                lallfound2 = (abs(jj).ge.1)
              else
                lallfound2 = (r2.gt.r2j.and.r2.gt.r12.and.nvec.eq.nveck0)
              endif
              r2j = r2
            enddo
          enddo
!
!  Increment by first vector
!
          ii = ii + idir
          rxi = rxi + rcx1
          ryi = ryi + rcy1
          rzi = rzi + rcz1
!
!  Check to see if this direction is complete
!
          if (ldcellr) then
            lallfound1 = (abs(ii).ge.1)
          else
            lallfound1 = (r2.gt.r2i.and.r2.gt.r12.and.nvec.eq.nvecj0)
          endif
          r2i = r2
        enddo
      enddo
!
!  End of cell search
!
    endif
  enddo
  ncores = nreg1old
!*********************************
!  Locate shells to match cores  *
!*********************************
  do j = 1,numat
    if (nat(j).gt.maxele) then
      xcrd = xclat(j) - xdc
      ycrd = yclat(j) - ydc
      zcrd = zclat(j) - zdc
      if (lmol) then
        nmj = natmol(j)
        indmj = nmolind(j)
        izj = (indmj/100) - 5
        ind = indmj - 100*(izj + 5)
        iyj = (ind/10) - 5
        ind = ind - 10*(iyj + 5)
        ixj = ind - 5
!
!  Adjust molecule indices for coordinate shift
!
        ijx = ixj + imid
        ijy = iyj + jmid
        ijz = izj + kmid
      endif
!
!  Search unit cells iteratively for all images
!
      xij = xcrd
      yij = ycrd
      zij = zcrd
!
!  Find projection of cell vector 3 on to i - j vector
!
      rnorm = xij*xij + yij*yij + zij*zij
      if (rnorm.gt.1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
      proj3 = rnorm*recipc*(xij*r3x + yij*r3y + zij*r3z)
      kmid = nint(proj3)
      xij = xij - kmid*r3x
      yij = yij - kmid*r3y
      zij = zij - kmid*r3z
!
!  Find projection of cell vector 2 on to i - j vector
!
      rnorm = xij*xij + yij*yij + zij*zij
      if (rnorm.gt.1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
      proj2 = rnorm*recipb*(xij*r2x + yij*r2y + zij*r2z)
      jmid = nint(proj2)
      xij = xij - jmid*r2x
      yij = yij - jmid*r2y
      zij = zij - jmid*r2z
!
!  Find projection of cell vector 1 on to i - j vector
!
      rnorm = xij*xij + yij*yij + zij*zij
      if (rnorm.gt.1.0d-12) rnorm = 1.0_dp/sqrt(rnorm)
      proj1 = rnorm*recipa*(xij*r1x + yij*r1y + zij*r1z)
      imid = nint(proj1)
      xij = xij - imid*r1x
      yij = yij - imid*r1y
      zij = zij - imid*r1z
!
!  Adjust molecule indices for coordinate shift
!
      if (lmol) then
        ijx = ixj + imid
        ijy = iyj + jmid
        ijz = izj + kmid
      endif
!
!  Initialise number of distances to zero
!
      nvec = 0
!
!  Outer loop over first cell vector direction
!
      do idir = 1,-1,-2
!
!  Reinitialise distance squared
!
        r2i = 10000.0_dp*r12s
!
!  Loop over first cell vector
!
        lallfound1 = .false.
        if (idir.eq.1) then
          ii = 0
        else
          ii = - 1
        endif
!
!  Set initial coordinate vector
!
        rxi = xij + dble(ii)*r1x
        ryi = yij + dble(ii)*r1y
        rzi = zij + dble(ii)*r1z
!
!  Set increment vector
!
        rcx1 = dble(idir)*r1x
        rcy1 = dble(idir)*r1y
        rcz1 = dble(idir)*r1z
!
        do while (.not.lallfound1)
!
!  Save number of vectors before search over second direction
!
          nvecj0 = nvec
!
!  Outer loop over second cell vector direction
!
          do jdir = 1,-1,-2
!
!  Reinitialise saved distance squared
!
            r2j = 10000.0_dp*r12s
!
!  Loop over second cell vector
!
            lallfound2 = .false.
            if (jdir.eq.1) then
              jj = 0
            else
              jj = - 1
            endif
!
!  Set initial coordinate vector
!
            rxj = rxi + dble(jj)*r2x
            ryj = ryi + dble(jj)*r2y
            rzj = rzi + dble(jj)*r2z
!
!  Set increment vector
!
            rcx2 = dble(jdir)*r2x
            rcy2 = dble(jdir)*r2y
            rcz2 = dble(jdir)*r2z
!
            do while (.not.lallfound2)
!
!  Save number of vectors before search over third direction
!
              nveck0 = nvec
!
!  Outer loop over third cell vector direction
!
              do kdir = 1,-1,-2
!
!  Reinitialise saved distance squared
!
                r2k = 10000.0_dp*r12s
!
!  Loop over third cell vector
!
                lallfound3 = .false.
                if (kdir.eq.1) then
                  kk = 0
                else
                  kk = - 1
                endif
!
!  Set initial coordinate vector
!
                rx = rxj + dble(kk)*r3x
                ry = ryj + dble(kk)*r3y
                rz = rzj + dble(kk)*r3z
!
!  Set increment vector
!
                rcx3 = dble(kdir)*r3x
                rcy3 = dble(kdir)*r3y
                rcz3 = dble(kdir)*r3z
!
                do while (.not.lallfound3)
!
!  Calculate square of distance
!
                  r2 = rx*rx + ry*ry + rz*rz
!
!  Check distance squared against cutoff squared
!
                  if (r2.le.r12s.or.ldcellr) then
                    nvec = nvec + 1
!
!  Check that matching core is in region1 unless ldcellr
!
                    lfound = .false.
                    i = 0
                    do while (i.lt.ncores.and..not.lfound)
                      i = i + 1
                      xcs = rx - xperf(i) + xdc
                      ycs = ry - yperf(i) + ydc
                      zcs = rz - zperf(i) + zdc
                      rcs = xcs*xcs + ycs*ycs + zcs*zcs
                      lfound = (rcs.le.cut2s)
                    enddo
                    if (lfound) then
!
!  Valid region 1 shell
!
                      nreg1old = nreg1old + 1
                      if (nreg1old.ge.maxr1at) then
                        maxr1at = nreg1old + 20
                        call changemaxr1at
                      endif
                      xperf(nreg1old) = rx + xdc
                      yperf(nreg1old) = ry + ydc
                      zperf(nreg1old) = rz + zdc
                      qp(nreg1old) = qf(j)
                      occp(nreg1old) = occuf(j)
                      natp(nreg1old) = nat(j)
                      ntypep(nreg1old) = nftype(j)
                      npsite(nreg1old) = j
                      qreg1old = qreg1old + qp(nreg1old)*occp(nreg1old)
!
!  Molecule handling
!
                      if (lmol) then
                        ixji = ii - ijx
                        iyji = jj - ijy
                        izji = kk - ijz
                        newind = (ixji + 5) + 10*(iyji + 5) + 100*(izji + 5)
                        ndefmolp(nreg1old) = nmj
                        if (nmj.eq.0) then
                          ndefindp(nreg1old) = 0
                        else
                          ndefindp(nreg1old) = newind
                        endif
                      endif
                    endif
                  endif
!
!  Increment by third vector
!
                  kk = kk + kdir
                  rx = rx + rcx3
                  ry = ry + rcy3
                  rz = rz + rcz3
!
!  Check to see if this direction is complete
!
                  if (ldcellr) then
                    lallfound3 = (abs(kk).ge.1)
                  else
                    lallfound3 = (r2.gt.r2k.and.r2.gt.r12s)
                  endif
                  r2k = r2
                enddo
              enddo
!
!  Increment by second vector
!
              jj = jj + jdir
              rxj = rxj + rcx2
              ryj = ryj + rcy2
              rzj = rzj + rcz2
!
!  Check to see if this direction is complete
!
              if (ldcellr) then
                lallfound2 = (abs(jj).ge.1)
              else
                lallfound2 = (r2.gt.r2j.and.r2.gt.r12s.and.nvec.eq.nveck0)
              endif
              r2j = r2
            enddo
          enddo
!
!  Increment by first vector
!
          ii = ii + idir
          rxi = rxi + rcx1
          ryi = ryi + rcy1
          rzi = rzi + rcz1
!
!  Check to see if this direction is complete
!
          if (ldcellr) then
            lallfound1 = (abs(ii).ge.1)
          else
            lallfound1 = (r2.gt.r2i.and.r2.gt.r12s.and.nvec.eq.nvecj0)
          endif
          r2i = r2
        enddo
      enddo
!
!  End of cell search
!
    endif
  enddo
!********************************
!  Generate defective region 1  *
!********************************
  if (lexpand) then
!
!  Remove ions from excluded region for region 1 list
!
!  Locate cores first
!
    do i = 1,nreg1old
      if (natp(i).le.maxele) then
        xd = xperf(i) - xdc
        yd = yperf(i) - ydc
        zd = zperf(i) - zdc
        r = xd*xd + yd*yd + zd*zd
        if (r.gt.r1ex) then
          nreg1 = nreg1 + 1
          lr1created(nreg1) = .false.
          ldfix(nreg1) = .false.
          xdefe(nreg1) = xperf(i)
          ydefe(nreg1) = yperf(i)
          zdefe(nreg1) = zperf(i)
          qdefe(nreg1) = qp(i)
          radefe(nreg1) = radcfg(nrelat(npsite(i)) + nsft)
          occdefe(nreg1) = occp(i)
          qreg1 = qreg1 + qa(nreg1)*occua(nreg1)
          natdefe(nreg1) = natp(i)
          ntypdefe(nreg1) = ntypep(i)
          ldefbsmat(nreg1) = lbsmat(nrelat(npsite(i)) + nsft)
          ldqmatom(nreg1) = .false.
          if (lmol) then
            ndefmol(nreg1) = ndefmolp(i)
            ndefind(nreg1) = ndefindp(i)
          endif
          nreldef(nreg1) = npsite(i)
        endif
      endif
    enddo
    ncores = nreg1
!
!  Locate matching shells for cores
!
    do i = 1,nreg1old
      if (natp(i).gt.maxele) then
        xd = xperf(i)
        yd = yperf(i)
        zd = zperf(i)
        j = 0
        lfound = .false.
        do while (j.lt.ncores.and..not.lfound)
          j = j + 1
          xcs = xd - xdefe(j)
          ycs = yd - ydefe(j)
          zcs = zd - zdefe(j)
          rcs = xcs*xcs + ycs*ycs + zcs*zcs
          lfound = (rcs.le.cut2s)
        enddo
        if (lfound) then
          nreg1 = nreg1 + 1
          lr1created(nreg1) = .false.
          xdefe(nreg1) = xperf(i)
          ydefe(nreg1) = yperf(i)
          zdefe(nreg1) = zperf(i)
          qdefe(nreg1) = qp(i)
          radefe(nreg1) = radcfg(nrelat(npsite(i)) + nsft)
          occdefe(nreg1) = occp(i)
          qreg1 = qreg1 + qa(nreg1)*occua(nreg1)
          natdefe(nreg1) = natp(i)
          ntypdefe(nreg1) = ntypep(i)
          ldefbsmat(nreg1) = lbsmat(nrelat(npsite(i)) + nsft)
          ldqmatom(nreg1) = .false.
          if (lmol) then
            ndefmol(nreg1) = ndefmolp(i)
            ndefind(nreg1) = ndefindp(i)
          endif
          nreldef(nreg1) = npsite(i)
        endif
      endif
    enddo
  elseif (.not.lrestore) then
    nreg1 = nreg1old
    qreg1 = qreg1old
!
!  Copy perfect region 1
!
    do i = 1,nreg1
      lr1created(i) = .false.
      ldfix(i) = .false.
      xdefe(i) = xperf(i)
      ydefe(i) = yperf(i)
      zdefe(i) = zperf(i)
      qdefe(i) = qp(i)
      radefe(i) = radcfg(nrelat(npsite(i)) + nsft)
      occdefe(i) = occp(i)
      natdefe(i) = natp(i)
      ntypdefe(i) = ntypep(i)
      ldefbsmat(i) = lbsmat(nrelat(npsite(i)) + nsft)
      ldqmatom(nreg1) = .false.
      if (lmol) then
        ndefmol(i) = ndefmolp(i)
        ndefind(i) = ndefindp(i)
      endif
      nreldef(i) = npsite(i)
    enddo
!
!  Set optimisation flags prior to restore
!
    if (index(keyword,' brea').ne.0.or.index(keyword,'brea').eq.1) then
      do i = 1,3*nreg1
        idopt(i) = 0
      enddo
    elseif (index(keyword,'shel').ne.0) then
      ind = 0
      do i = 1,nreg1
        if (natdefe(i).gt.maxele) then
          idopt(ind + 1) = 1
          idopt(ind + 2) = 1
          idopt(ind + 3) = 1
        else
          idopt(ind + 1) = 0
          idopt(ind + 2) = 0
          idopt(ind + 3) = 0
        endif
        ind = ind + 3
      enddo
    else
      do i = 1,3*nreg1
        idopt(i) = 1
      enddo
    endif
  endif
!
  if (lrestore) then
!**************************************
!  Set based on configuration arrays  *
!**************************************
!
!  If defect list is present change mode to 1 or 2
!
    if (mode2a.ge.3.and..not.ldeflin(ncf).and..not.lindvacptr(ncf).and..not.lindintptr(ncf)) mode2a = mode2a - 2
    if (mode2a.ge.3.and..not.ldeflin(ncf).and..not.lindvacptr(ncf).and..not.lindintptr(ncf)) mode2a = 2
!
!  Find start of configuration in array
!
    ndind = 0
    niind = 0
    nvind = 0
    do i = 1,ncf-1
      ndind = ndind + nreg1cfg(i)
      niind = niind + nintecfg(i)
      nvind = nvind + nvacacfg(i)
    enddo
    nreg1 = nreg1cfg(ncf)
    ind = 0
    do i = 1,nreg1cfg(ncf)
      natdefe(i) = natdefecfg(ndind+i)
      ntypdefe(i) = ntypdefecfg(ndind+i)
      nreldef(i) = nreldefcfg(ndind+i)
      xdefe(i) = xdefecfg(ndind+i)
      ydefe(i) = ydefecfg(ndind+i)
      zdefe(i) = zdefecfg(ndind+i)
      qdefe(i) = qdefecfg(ndind+i)
      occdefe(i) = occdefecfg(ndind+i)
      radefe(i) = radefecfg(ndind+i)
      ndefmol(i) = ndefmolcfg(ndind+i)
      ndefind(i) = ndefindcfg(ndind+i)
      ldefbsmat(i) = ldefbsmatcfg(ndind+i)
      ldqmatom(i) = ldqmatomcfg(ndind+i)
      idopt(ind+1) = idoptcfg(1,ndind+i)
      idopt(ind+2) = idoptcfg(2,ndind+i)
      idopt(ind+3) = idoptcfg(3,ndind+i)
      ind = ind + 3
      qreg1 = qreg1 + qdefe(i)*occdefe(i)
      lr1created(i) = .false.
    enddo
!
    ninte = nintecfg(ncf)
    nvaca = nvacacfg(ncf)
    do i = 1,nvaca
      ndvacptr(i) = ndvacptrcfg(nvind+i)
    enddo
    do i = 1,ninte
      ndintptr(i) = ndintptrcfg(niind+i)
    enddo
!
!  Check all species are present
!
    do i = 1,nreg1
      nati = natdefe(i)
      ntypi = ntypdefe(i)
      lfind = .false.
      do j = 1,nspec
        if (nati.eq.natspec(j).and.(ntypi.eq.ntypspec(j).or.ntypspec(j).eq.0)) lfind = .true.
      enddo
      if (.not.lfind) then
        nspec = nspec + 1
        natspec(nspec) = nati
        ntypspec(nspec) = ntypi
        qlspec(nspec) = qa(i)
        radspec(nspec) = rada(i)
        lbrspec(nspec) = ldefbsmat(i)
      endif
    enddo
  else
!
!  Process defect list adjusting region 1
!
    do i = 1,nldef
      ndt = ndeftyp(ndefst + i)
      lcore = .true.
      lshel = .true.
      if (ndt.lt.10) then
!************
!  Vacancy  *
!************
        if (ndt.eq.1) then
          ni = nint(xdef(ndefst + i))
          xv = xclat(nrel2(ni))
          yv = yclat(nrel2(ni))
          zv = zclat(nrel2(ni))
          call nearesti(xdc,ydc,zdc,xv,yv,zv,r1xf,r1yf,r1zf,r2xf,r2yf,r2zf,r3xf,r3yf,r3zf,ix,iy,iz)
        elseif (ndt.eq.2) then
          inat = ndefnat(ndefst + i)
          itype = ndeftp(ndefst + i)
          if (inat.gt.2*maxele) then
            lcore = .true.
            lshel = .true.
            inat = inat - 2*maxele
          elseif (inat.gt.maxele) then
            lcore = .false.
            lshel = .true.
          else
            lcore = .true.
            lshel = .false.
          endif
          ni = 0
          ni2 = 0
          do k = 1,numat
            if (ni.eq.0) then
              if (nat(k).eq.inat.and.(itype.eq.0.or.itype.eq.nftype(k))) ni = k
            elseif (ni2.eq.0) then
              if (nat(k).eq.inat.and.(itype.eq.0.or.itype.eq.nftype(k))) then
                ni2 = k
                nwarn = nwarn + 1
                call outwarning('Ambiguous vacancy specifier used',0_i4)
              endif
            endif
          enddo
          xv = xclat(ni)
          yv = yclat(ni)
          zv = zclat(ni)
          call nearesti(xdc,ydc,zdc,xv,yv,zv,r1xf,r1yf,r1zf,r2xf,r2yf,r2zf,r3xf,r3yf,r3zf,ix,iy,iz)
        elseif (ndt.eq.3) then
          xi = xdef(ndefst + i)
          yi = ydef(ndefst + i)
          zi = zdef(ndefst + i)
          xv = xi*r1xf + yi*r2xf + zi*r3xf
          yv = xi*r1yf + yi*r2yf + zi*r3yf
          zv = xi*r1zf + yi*r2zf + zi*r3zf
          call nearesti(xdc,ydc,zdc,xv,yv,zv,r1xf,r1yf,r1zf,r2xf,r2yf,r2zf,r3xf,r3yf,r3zf,ix,iy,iz)
        elseif (ndt.eq.4) then
          xv = xdef(ndefst + i)
          yv = ydef(ndefst + i)
          zv = zdef(ndefst + i)
        else
!
!  Find atoms with correct molecule no. for central cell
!
          nmi = nint(xdef(ndefst + i))
          nv = 0
          do j = 1,nreg1
            if (ndefmol(j).eq.nmi.and.ndefind(j).eq.555.and..not.lr1created(j)) then
              nv = nv + 1
              nptrr1(nv) = j
            endif
          enddo
!
!  Find centroid of cell
!
          if (nv.gt.0) then
            xv = 0.0_dp
            yv = 0.0_dp
            zv = 0.0_dp
            do j = 1,nv
              xv = xv + xdefe(nptrr1(j))
              yv = yv + ydefe(nptrr1(j))
              zv = zv + zdefe(nptrr1(j))
            enddo
            xv = xv/nv
            yv = yv/nv
            zv = zv/nv
!
!  Find image of molecule nearest to defect centre
!
            call nearesti(xdc,ydc,zdc,xv,yv,zv,r1xf,r1yf,r1zf,r2xf,r2yf,r2zf,r3xf,r3yf,r3zf,ix,iy,iz)
            ind = (ix + 5) + 10*(iy + 5) + 100*(iz + 5)
!
!  Find atoms with correct molecule no. and index
!
            nv = 0
            do j = 1,nreg1
              if (ndefmol(j).eq.nmi.and.ndefind(j).eq.ind.and..not.lr1created(j)) then
                nv = nv + 1
                nptrr1(nv) = j
              endif
            enddo
          endif
        endif
        if (ndt.ne.5) then
!
!  Check that site is in region 1
!
          rtest = (xv-xdc)**2 + (yv-ydc)**2 + (zv-zdc)**2
          if (rtest.gt.r12) then
            call outerror('vacancy is outside region 1',0_i4)
            call stopnow('setdef')
          endif
          nv = 0
          do j = 1,nreg1
            nj = natdefe(j)
            xdiff = xdefe(j) - xv
            ydiff = ydefe(j) - yv
            zdiff = zdefe(j) - zv
            r2 = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff
            if (lcore.and.nj.le.maxele) then
              if (r2.le.rfind.and..not.lr1created(j)) then
                nv = nv + 1
                nptrr1(nv) = j
              endif
            elseif (lshel.and.nj.gt.maxele) then
              if (r2.le.rfind.and..not.lr1created(j)) then
                nv = nv + 1
                nptrr1(nv) = j
              endif
            endif
          enddo
        endif
        if (nv.gt.0) then
          do j = 1,nv
            call remove(nptrr1(j),qreg1,lr1created)
            do k = j + 1,nv
              nptrr1(k) = nptrr1(k) - 1
            enddo
          enddo
        endif
      elseif (ndt.lt.20) then
!*************
!  Impurity  *
!*************
        if (ndt.eq.11.or.ndt.eq.12) then
          ni = nint(xdef(ndefst + i))
          xv = xclat(nrel2(ni))
          yv = yclat(nrel2(ni))
          zv = zclat(nrel2(ni))
          call nearesti(xdc,ydc,zdc,xv,yv,zv,r1xf,r1yf,r1zf,r2xf,r2yf,r2zf,r3xf,r3yf,r3zf,ix,iy,iz)
        elseif (ndt.eq.13) then
          xi = xdef(ndefst + i)
          yi = ydef(ndefst + i)
          zi = zdef(ndefst + i)
          xv = xi*r1xf + yi*r2xf + zi*r3xf
          yv = xi*r1yf + yi*r2yf + zi*r3yf
          zv = xi*r1zf + yi*r2zf + zi*r3zf
          call nearesti(xdc,ydc,zdc,xv,yv,zv,r1xf,r1yf,r1zf,r2xf,r2yf,r2zf,r3xf,r3yf,r3zf,ix,iy,iz)
        elseif (ndt.eq.14) then
          xv = xdef(ndefst + i)
          yv = ydef(ndefst + i)
          zv = zdef(ndefst + i)
        endif
        nv = 0
!
!  Check that site is in region 1
!
        rtest = (xv-xdc)**2 + (yv-ydc)**2 + (zv-zdc)**2
        if (rtest.gt.r12) then
          call outerror('impurity is outside region 1',0_i4)
          call stopnow('setdef')
        endif
        do j = 1,nreg1
          xdiff = xdefe(j) - xv
          ydiff = ydefe(j) - yv
          zdiff = zdefe(j) - zv
          r2 = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff
          if (r2.le.rfind.and..not.lr1created(j)) then
            nv = nv + 1
            nptrr1(nv) = j
          endif
        enddo
        if (nv.gt.0) then
          do j = 1,nv
            call remove(nptrr1(j),qreg1,lr1created)
            do k = j + 1,nv
              nptrr1(k) = nptrr1(k) - 1
            enddo
          enddo
        endif
        nreg1 = nreg1 + 1
        if (nreg1.ge.maxr1at) then
          maxr1at = nreg1 + 20
          call changemaxr1at
        endif
        inat = ndefnat(ndefst + i)
        lbre = .false.
        if (inat.gt.3*maxele) then
          lbre = .true.
          inat = inat - 3*maxele
        endif
        if (inat.gt.2*maxele) then
          natdefe(nreg1) = inat - 2*maxele
        else
          natdefe(nreg1) = inat
        endif
        ntypdefe(nreg1) = ndeftp(ndefst + i)
        xdefe(nreg1) = xv
        ydefe(nreg1) = yv
        zdefe(nreg1) = zv
        lr1created(nreg1) = .true.
        ldfix(nreg1) = ldeffix(ndefst + i)
        inddfix(nreg1) = inddeffix(ndefst + i)
        q = 0.0_dp
        ra = 0.0_dp
        lfind = .false.
        do k = 1,nspec
          if (natdefe(nreg1).eq.natspec(k).and.(ntypspec(k).eq.0.or.ntypspec(k).eq.ntypdefe(nreg1))) then
            q = qlspec(k)
            ra = radspec(k)
            lbre = lbrspec(k)
            lfind = .true.
          endif
        enddo
        if (.not.lfind) then
          nspec = nspec + 1
          natspec(nspec) = natdefe(nreg1)
          ntypspec(nspec) = ntypdefe(nreg1)
          qlspec(nspec) = q
          radspec(nspec) = ra
          lbrspec(nspec) = lbre
        endif
        qdefe(nreg1) = q
        radefe(nreg1) = ra
        ldefbsmat(nreg1) = lbre
        ldqmatom(nreg1) = .false.
        occdefe(nreg1) = 1.0_dp
        nreldef(nreg1) = 0
        qreg1 = qreg1 + q
        if (inat.gt.2*maxele) then
          inat = inat - maxele
          lshel = .false.
          do k = 1,nspec
            if (inat.eq.natspec(k).and.(ntypspec(k).eq.0.or.ntypspec(k).eq.ntypdefe(nreg1))) then
              q = qlspec(k)
              ra = radspec(k)
              lbre = lbrspec(k)
              lshel = .true.
            endif
          enddo
          if (lshel) then
            nreg1 = nreg1 + 1
            if (nreg1.ge.maxr1at) then
              maxr1at = nreg1 + 20
              call changemaxr1at
            endif
            natdefe(nreg1) = inat
            ntypdefe(nreg1) = ndeftp(ndefst + i)
            xdefe(nreg1) = xv
            ydefe(nreg1) = yv
            zdefe(nreg1) = zv
            lr1created(nreg1) = .true.
            ldfix(nreg1) = ldeffix(ndefst + i)
            inddfix(nreg1) = inddeffix(ndefst + i)
            qdefe(nreg1) = q
            radefe(nreg1) = ra
            ldefbsmat(nreg1) = lbre
            ldqmatom(nreg1) = .false.
            occdefe(nreg1) = 1.0_dp
            nreldef(nreg1) = 0
            qreg1 = qreg1 + q
          endif
        endif
      else
!*****************
!  Interstitial  *
!*****************
        nreg1 = nreg1 + 1
        if (nreg1.ge.maxr1at) then
          maxr1at = nreg1 + 20
          call changemaxr1at
        endif
        inat = ndefnat(ndefst + i)
        lbre = .false.
        if (inat.gt.3*maxele) then
          lbre = .true.
          inat = inat-3*maxele
        endif
        if (inat.gt.2*maxele) then
          natdefe(nreg1) = inat - 2*maxele
        else
          natdefe(nreg1) = inat
        endif
        ntypdefe(nreg1) = ndeftp(ndefst + i)
        lr1created(nreg1) = .true.
        ldfix(nreg1) = ldeffix(ndefst + i)
        inddfix(nreg1) = inddeffix(ndefst + i)
        q = 0.0_dp
        ra = 0.0_dp
        lfind = .false.
        do k = 1,nspec
          if (natdefe(nreg1).eq.natspec(k).and.(ntypspec(k).eq.0.or.ntypspec(k).eq.ntypdefe(nreg1))) then
            q = qlspec(k)
            ra = radspec(k)
            lbre = lbrspec(k)
            lfind = .true.
          endif
        enddo
        if (.not.lfind) then
          nspec = nspec + 1
          natspec(nspec) = natdefe(nreg1)
          ntypspec(nspec) = ntypdefe(nreg1)
          qlspec(nspec) = q
          radspec(nspec) = ra
          lbrspec(nspec) = lbre
        endif
        qdefe(nreg1) = q
        radefe(nreg1) = ra
        ldefbsmat(nreg1) = lbre
        ldqmatom(nreg1) = .false.
        occdefe(nreg1) = 1.0_dp
        nreldef(nreg1) = 0
        qreg1 = qreg1 + q
        if (ndt.eq.21) then
!
!  Fractional coordinates
!
          xi = xdef(ndefst + i)
          yi = ydef(ndefst + i)
          zi = zdef(ndefst + i)
          xv = xi*r1xf + yi*r2xf + zi*r3xf
          yv = xi*r1yf + yi*r2yf + zi*r3yf
          zv = xi*r1zf + yi*r2zf + zi*r3zf
          call nearesti(xdc,ydc,zdc,xv,yv,zv,r1xf,r1yf,r1zf,r2xf,r2yf,r2zf,r3xf,r3yf,r3zf,ix,iy,iz)
          xdefe(nreg1) = xv
          ydefe(nreg1) = yv
          zdefe(nreg1) = zv
        elseif (ndt.eq.22) then
!
!  Cartesian coordinates
!
          xdefe(nreg1) = xdef(ndefst + i)
          ydefe(nreg1) = ydef(ndefst + i)
          zdefe(nreg1) = zdef(ndefst + i)
        elseif (ndt.eq.23) then
!
!  Bond directive - symbol
!  locate species nearest defect centre with valid symbol
!
          inat = nint(xdef(ndefst + i))
          ityp = nint(ydef(ndefst + i))
          rmax = 100.0_dp
          iptr = 0
          do j = 1,nreg1-1
            if (natdefe(j).eq.inat) then
              jtyp = ntypdefe(j)
              if (ityp.eq.0.or.ityp.eq.jtyp) then
                xci = xdefe(j) - xdc
                yci = ydefe(j) - ydc
                zci = zdefe(j) - zdc
                r = xci*xci + yci*yci + zci*zci
                if (r.lt.rmax) then
                  iptr = j
                  rmax = r
                endif
              endif
            endif
          enddo
          if (iptr.eq.0) then
            call outerror('no valid species in region 1 for bonded interstitial',0_i4)
            call stopnow('setdef')
          endif
          call bondint(iptr)
        elseif (ndt.eq.24.or.ndt.eq.25) then
!
!  Bond directive - fractional or cartesian coords
!
          if (ndt.eq.24) then
            xi = xdef(ndefst + i)
            yi = ydef(ndefst + i)
            zi = zdef(ndefst + i)
            xv = xi*r1xf + yi*r2xf + zi*r3xf
            yv = xi*r1yf + yi*r2yf + zi*r3yf
            zv = xi*r1zf + yi*r2zf + zi*r3zf
            call nearesti(xdc,ydc,zdc,xv,yv,zv,r1xf,r1yf,r1zf,r2xf,r2yf,r2zf,r3xf,r3yf,r3zf,ix,iy,iz)
          else
            xv = xdef(ndefst + i)
            yv = ydef(ndefst + i)
            zv = zdef(ndefst + i)
          endif
          rmax = 100.0_dp
          iptr = 1
          do j = 1,nreg1-1
            xci = xdefe(j) - xv
            yci = ydefe(j) - yv
            zci = zdefe(j) - zv
            r = xci*xci + yci*yci + zci*zci
            if (r.lt.rmax.and.natdefe(j).le.maxele) then
              iptr = j
              rmax = r
            endif
          enddo
          call bondint(iptr)
        endif
!
!  Check that site is in region 1
!
        rtest = (xdefe(nreg1)-xdc)**2 + (ydefe(nreg1)-ydc)**2 + (zdefe(nreg1)-zdc)**2
        if (rtest.gt.r12) then
          call outerror('interstitial is outside region 1',0_i4)
          call stopnow('setdef')
        endif
        if (inat.gt.2*maxele) then
          inat = inat-maxele
          lshel = .false.
          do k = 1,nspec
            if (inat.eq.natspec(k).and.(ntypspec(k).eq.0.or.ntypspec(k).eq.ntypdefe(nreg1))) then
              q = qlspec(k)
              ra = radspec(k)
              lbre = lbrspec(k)
              lshel = .true.
            endif
          enddo
          if (lshel) then
            nreg1 = nreg1 + 1
            if (nreg1.ge.maxr1at) then
              maxr1at = nreg1 + 20
              call changemaxr1at
            endif
            natdefe(nreg1) = inat
            ntypdefe(nreg1) = ndeftp(ndefst + i)
            lr1created(nreg1) = .true.
            ldfix(nreg1) = ldeffix(ndefst + i)
            inddfix(nreg1) = inddeffix(ndefst + i)
            qdefe(nreg1) = q
            ldefbsmat(nreg1) = lbre
            ldqmatom(nreg1) = .false.
            radefe(nreg1) = ra
            occdefe(nreg1) = 1.0_dp
            nreldef(nreg1) = 0
            qreg1 = qreg1 + q
            xdefe(nreg1) = xdefe(nreg1-1)
            ydefe(nreg1) = ydefe(nreg1-1)
            zdefe(nreg1) = zdefe(nreg1-1)
          endif
        endif
      endif
    enddo
  endif
!
!  Assign net charge of defect
!
  qdef = qreg1 - qreg1old
!
!  Assign molecule numbers to interstitial atoms
!
  if (lmol) call setmoldef(lr1created)
!
!  Set flag for BSM in region 1
!
  ldbsm = .false.
  i = 0
  do while (.not.ldbsm.and.i.lt.nreg1)
    i = i + 1
    ldbsm = ldefbsmat(i)
  enddo
!
!  Sort atoms in region 1
!
  call sort1(lnosort,lexpand)
!
!  For parallel runs set up distribution of region 1 atoms
!
  if (nprocs.gt.1) then
    nreg1onnode = 0
    nreg1onnodec = 0
    nreg1onnodes = 0
    node = 0
    icount = 0
!
!  If block size hasn't been input then choose a value based on the number of atoms versus processors
!
    if (nblocksize.eq.0) then
      nblocksize = 1
    endif
!
    do i = 1,nreg1
      icount = icount + 1
      reg12node(i) = node
      if (node.eq.procid) then
        nreg1onnode = nreg1onnode + 1
        node2reg1(nreg1onnode) = i
        reg12local(i) = nreg1onnode
        if (natdefe(i).le.maxele) then
          nreg1onnodec = nreg1onnodec + 1
        else
          nreg1onnodes = nreg1onnodes + 1
        endif
      else
        reg12local(i) = 0
      endif
      if (icount.eq.nblocksize) then
        icount = 0
        node = node + 1
        if (node.eq.nprocs) node = 0
      endif
    enddo
!
!  Add data for storage of region 2 block in nreg1 + 1
!
    i = nreg1 + 1
    reg12node(i) = node
    if (node.eq.procid) then
      node2reg1(nreg1onnode+1) = i
      reg12local(i) = nreg1onnode + 1
    else
      reg12local(i) = 0
    endif
  else
    nreg1onnode = nreg1
    do i = 1,nreg1
      node2reg1(i) = i
      reg12local(i) = i
      reg12node(i) = 0
    enddo
!
!  Add data for storage of region 2 block in nreg1 + 1
!
    i = nreg1 + 1
    reg12node(i) = 0
    node2reg1(nreg1+1) = i
    reg12local(i) = nreg1 + 1
  endif
!
!  Restore tasks
!
  if (lrestore) then
    allocate(itmp(3*nreg1),stat=status)
    if (status/=0) call outofmemory('setdef','itmp')
    do i = 1,3*nreg1
      itmp(i) = idopt(i)
    enddo
    nvar = 0
    do i = 1,3*nreg1
      if (itmp(i).eq.1) then
        nvar = nvar + 1
        idopt(nvar) = i
      endif
    enddo
    mvar = 3*nreg1
    do i = 1,nreg1
      if (ldefbsmat(i)) then
        nvar = nvar + 1
        idopt(nvar) = i + mvar
      endif
    enddo
!
!  Free local memory
!
    deallocate(itmp,stat=status)
    if (status/=0) call deallocate_error('setdef','itmp')
  endif
!**************************************
!  Set memory for defect calculation  *
!**************************************
  call realloc(x0,4_i4*nreg1,ierror)
  if (ierror.ne.0) call outofmemory('setdef','x0')
!***************************
!  Set optimisation flags  *
!***************************
!
!  Including correct for fixed impurities/interstitials
!
  if (.not.lflags.or..not.lrestore) then
    nfix = 0
    do i = 1,ndef
      if (ldeffix(i)) nfix = nfix + 1
    enddo
    if (index(keyword,' brea').ne.0.or.index(keyword,'brea').eq.1) then
      nvar = 0
    elseif (index(keyword,'shel').ne.0) then
      if (nfix.gt.0) then
        nvar = 0
        ix = 3*ncoreg1-2
        iy = ix + 1
        iz = ix + 2
        do i = 1,nreg1
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          if (.not.ldfix(i)) then
            nvar = nvar + 1
            idopt(nvar) = ix
            nvar = nvar + 1
            idopt(nvar) = iy
            nvar = nvar + 1
            idopt(nvar) = iz
          else
            if (inddfix(i).eq.1) then
              nvar = nvar + 1
              idopt(nvar) = iy
              nvar = nvar + 1
              idopt(nvar) = iz
            elseif (inddfix(i).eq.2) then
              nvar = nvar + 1
              idopt(nvar) = ix
              nvar = nvar + 1
              idopt(nvar) = iz
            elseif (inddfix(i).eq.3) then
              nvar = nvar + 1
              idopt(nvar) = ix
              nvar = nvar + 1
              idopt(nvar) = iy
            elseif (inddfix(i).eq.4) then
              nvar = nvar + 1
              idopt(nvar) = iz
            elseif (inddfix(i).eq.5) then
              nvar = nvar + 1
              idopt(nvar) = iy
            elseif (inddfix(i).eq.6) then
              nvar = nvar + 1
              idopt(nvar) = ix
            endif
          endif
        enddo
      else
        nvar = 3*nshreg1
        do i = 1,3*nshreg1
          idopt(i) = 3*ncoreg1 + i
        enddo
      endif
    else
      if (nfix.gt.0) then
        nvar = 0
        ix = -2
        iy = -1
        iz = 0
        do i = 1,nreg1
          ix = ix + 3
          iy = iy + 3
          iz = iz + 3
          if (.not.ldfix(i)) then
            nvar = nvar + 1
            idopt(nvar) = ix
            nvar = nvar + 1
            idopt(nvar) = iy
            nvar = nvar + 1
            idopt(nvar) = iz
          else
            if (inddfix(i).eq.1) then
              nvar = nvar + 1
              idopt(nvar) = iy
              nvar = nvar + 1
              idopt(nvar) = iz
            elseif (inddfix(i).eq.2) then
              nvar = nvar + 1
              idopt(nvar) = ix
              nvar = nvar + 1
              idopt(nvar) = iz
            elseif (inddfix(i).eq.3) then
              nvar = nvar + 1
              idopt(nvar) = ix
              nvar = nvar + 1
              idopt(nvar) = iy
            elseif (inddfix(i).eq.4) then
              nvar = nvar + 1
              idopt(nvar) = iz
            elseif (inddfix(i).eq.5) then
              nvar = nvar + 1
              idopt(nvar) = iy
            elseif (inddfix(i).eq.6) then
              nvar = nvar + 1
              idopt(nvar) = ix
            endif
          endif
        enddo
      else
        nvar = 3*nreg1
        do i = 1,3*nreg1
          idopt(i) = i
        enddo
      endif
    endif
  endif
  if (index(keyword,'nobr').eq.0) then
!
!  Need to add breathing flags
!
    do i = 1,nreg1
      if (ldefbsmat(i)) then
        nvar = nvar + 1
        idopt(nvar) = 3*nreg1 + i
      endif
    enddo
  endif
!
!  Place relevant quantities in configuration arrays
!
  if (nreg1tot+nreg1.gt.maxtotr1at) then
    maxtotr1at = nreg1tot + nreg1
    call changemaxtotr1at
  endif
  nreg1cfg(ncf) = nreg1
  nreg1tot = nreg1tot + nreg1
!
!  Find start of atoms in defect configuration arrays
!
  ndind = 0
  do i = 1,ncf-1
    ndind = ndind + nreg1cfg(i)
  enddo
  do i = 1,nreg1
    natdefecfg(ndind+i) = natdefe(i)
    ntypdefecfg(ndind+i) = ntypdefe(i)
    xdefecfg(ndind+i) = xdefe(i)
    ydefecfg(ndind+i) = ydefe(i)
    zdefecfg(ndind+i) = zdefe(i)
    qdefecfg(ndind+i) = qdefe(i)
    occdefecfg(ndind+i) = occdefe(i)
    radefecfg(ndind+i) = radefe(i)
    ndefmolcfg(ndind+i) = ndefmol(i)
    ndefindcfg(ndind+i) = ndefind(i)
    ldefbsmatcfg(ndind+i) = ldefbsmat(i)
    ldqmatomcfg(ndind+i) = ldqmatom(i)
  enddo
!
  return
  end
