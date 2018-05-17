  subroutine real1Dd(erealin,esregion12,esregion2,lgrad1,lgrad2)
!
!  This subroutine calculates the electrostatic energy of a 1-D system 
!  in real space based on the algorithm implemented in CRYSTAL. A 
!  neutralising uniform background charge density is applied and then 
!  subtracted again. Because a sum over neutral unit cells is required, 
!  this code is kept separate from the other real space routines. The 
!  conventional real space routines are used to handle the Coulomb
!  subtraction issues in order to keep this routine simple.
!  Distributed second derivatives version. Freezing not allowed.
!
!   5/13 Created from real1D
!   9/16 cputime renamed to g_cpu_time
!   9/16 constants renamed to g_constants
!   7/17 Calls to d1charge changed to d1charged
!   7/17 Self terms removed since they duplicate
!   7/17 Site energies added
!   7/17 Error in energy and convergence checking fixed
!   8/17 Incorrect change to handling of EM terms reversed
!   8/17 Site energies and eregion corrected for hfunc and emfunc terms
!   8/17 Factor of half introduced to correct sderv2 
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
!  No claim is made that this program is free from errors and
!  no liability will be accepted for any loss or damage that
!  may result. The user is responsible for checking the validity
!  of their results.
!
!  Copyright Curtin University 2017
!
!  Julian Gale, CIC, Curtin University, August 2017
!
  use g_constants,      only : angstoev
  use configurations,   only : nregionno, nregions, nregiontype, QMMMmode
  use control,          only : lnoreal, lseok, keyword, lDoQDeriv1, lDoQDeriv2, latomicstress
  use current
  use derivatives
  use element,          only : maxele
  use energies,         only : eregion2region, siteenergy
  use general,          only : nmaxcells, nemorder, smallself
  use iochannels,       only : ioout
  use kspace,           only : accf1D
  use parallel,         only : natomsonnode, node2atom, ioproc, nprocs
  use qmedata,          only : maxloop
  use shells,           only : cuts
  use symmetry,         only : lstr
  use times
  implicit none
!
!  Passed arguments
!
  real(dp), intent(inout) :: erealin
  real(dp), intent(inout) :: esregion12
  real(dp), intent(inout) :: esregion2
  logical,  intent(in)    :: lgrad1
  logical,  intent(in)    :: lgrad2
!
!  Local variables
!
  integer                 :: i
  integer                 :: ii
  integer                 :: ix
  integer                 :: iy
  integer                 :: iz 
  integer                 :: j 
  integer                 :: jx
  integer                 :: jy
  integer                 :: jz 
  integer                 :: m
  integer                 :: nati
  integer                 :: natj
  integer                 :: nregioni
  integer                 :: nregionj
  integer                 :: nregiontypi
  integer                 :: nregiontypj
  logical                 :: lconverged
  logical                 :: lcspair
  logical                 :: lopi
  logical                 :: lopj
  logical                 :: lreg2one
  logical                 :: lreg2pair
  real(dp)                :: accf
  real(dp)                :: acell
  real(dp)                :: g_cpu_time
  real(dp)                :: cut2s
  real(dp)                :: d0
  real(dp)                :: d0i
  real(dp)                :: d0j
  real(dp)                :: d0term
  real(dp)                :: d1
  real(dp)                :: d1i
  real(dp)                :: d1is
  real(dp)                :: d1ix
  real(dp)                :: d1iy
  real(dp)                :: d1iz
  real(dp)                :: d1j
  real(dp)                :: d1js
  real(dp)                :: d1jx
  real(dp)                :: d1jy
  real(dp)                :: d1jz
  real(dp)                :: d1s
  real(dp)                :: d1snoq
  real(dp)                :: d2
  real(dp)                :: d2i2
  real(dp)                :: d2ij
  real(dp)                :: d2j2
  real(dp)                :: dh1(3)
  real(dp)                :: dh2(3)
  real(dp)                :: dh1s
  real(dp)                :: dh2s
  real(dp)                :: d2h1(6)
  real(dp)                :: d2h2(6)
  real(dp)                :: d2h1m(3)
  real(dp)                :: d2h2m(3)
  real(dp)                :: d2h1s
  real(dp)                :: d2h2s
  real(dp)                :: d3h1(10)
  real(dp)                :: d3h1m(6)
  real(dp)                :: ediff
  real(dp)                :: elast
  real(dp)                :: ereal
  real(dp)                :: esum
  real(dp)                :: esumem
  real(dp)                :: esumh
  real(dp)                :: esum12
  real(dp)                :: esumem12
  real(dp)                :: esumh12
  real(dp)                :: esum2
  real(dp)                :: esumem2
  real(dp)                :: esumh2
  real(dp)                :: etrm
  real(dp)                :: e1
  real(dp)                :: e2
  real(dp)                :: h1
  real(dp)                :: h2
  real(dp)                :: lna
  real(dp)                :: oci
  real(dp)                :: ocj
  real(dp)                :: qi
  real(dp)                :: qj
  real(dp)                :: qij
  real(dp)                :: r
  real(dp)                :: rcut
  real(dp)                :: rr
  real(dp)                :: t1, t2
  real(dp)                :: u
  real(dp)                :: x
  real(dp)                :: y
  real(dp)                :: z
!
  t1 = g_cpu_time()
!
!  If noreal specified, return
!
  if (lnoreal) then
    erealin = 0.0_dp
    return
  endif
!********************************************************
!  Calculate Coulomb sum converged to desired accuracy  *
!********************************************************
!
!  Loop over number of cells in sum
!
  accf = 10.0**(-accf1D)
  lna = log(a)
  cut2s = cuts*cuts
  m = - 1
  lconverged = .false.
  elast = 0.0_dp
  esum = 0.0_dp
  esum12 = 0.0_dp
  esum2 = 0.0_dp
!
  if (index(keyword,'verb').ne.0) then
    if (ioproc) then
      write(ioout,'(/,'' Convergence of 1-D electrostatic summation:'',/)')
    endif
  endif
  do while (m.lt.nmaxcells.and..not.lconverged) 
    m = m + 1
!********************************************
!  Direct sum component over neutral cells  *
!********************************************
    acell = dble(m)*a
    ix = -2
    iy = -1
    iz =  0
    do ii = 1,natomsonnode
      i = node2atom(ii)
      nati = nat(i)
      oci = occuf(i)
      qi = qf(i)*oci
      nregioni = nregionno(nsft+i)
      nregiontypi = nregiontype(nregioni,ncf)
      lopi = .true.
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      jx = -2
      jy = -1
      jz =  0
      jloop: do j = 1,numat
        natj = nat(j)
        ocj = occuf(j)
        qj = qf(j)*ocj
        nregionj = nregionno(nsft+j)
        nregiontypj = nregiontype(nregionj,ncf)
        lopj = .true.
        lcspair = (abs(nati-natj).eq.maxele.or.(oci+ocj).lt.1.0001d0)
        if (lcspair) then
          rcut = cut2s
        else
          rcut = smallself
        endif
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3  
!
!  Set region 2 pair flags
!        
        lreg2one  = .false.
        lreg2pair = .false.
        if (lseok.and.nregions(ncf).ge.2) then
          lreg2pair = (nregioni.eq.2.and.nregionj.eq.2)
          if (.not.lreg2pair) lreg2one = (nregioni.eq.2.or.nregionj.eq.2)
        endif
!
!  QM/MM handling : i & j are both QM atoms => exclude
!
        if (QMMMmode(ncf).gt.0) then
          if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop
!           
!  QM/MM : Electrostatic embedding : If either i or j are QM atoms => exclude 
!       
          if (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1)) cycle jloop
        endif
!
        x = acell + xclat(j) - xclat(i)
        y = yclat(j) - yclat(i)
        z = zclat(j) - zclat(i)
        r = x*x + y*y + z*z
        if (r.gt.rcut) then
          r = sqrt(r)
          rr = 1.0_dp/r
          d0 = qi*qj*rr
          if (lreg2one) then
            esum12 = esum12 + d0
          elseif (lreg2pair) then
            esum2 = esum2 + d0
          else
            esum = esum + d0
          endif
!
          eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + 0.5_dp*d0*angstoev
!
          siteenergy(i) = siteenergy(i) + 0.5_dp*d0*angstoev
!
          if (lgrad1) then
            d1 = d0*rr*rr*angstoev
            xdrv(i) = xdrv(i) + d1*x
            ydrv(i) = ydrv(i) + d1*y
            zdrv(i) = zdrv(i) + d1*z
            if (nregioni.ne.nregionj) then
              xregdrv(nregioni) = xregdrv(nregioni) + d1*x
              yregdrv(nregioni) = yregdrv(nregioni) + d1*y
              zregdrv(nregioni) = zregdrv(nregioni) + d1*z
            endif
            if (lstr) then
              rstrd(1) = rstrd(1) - 0.5_dp*d1*x*x
              if (latomicstress) then
                atomicstress(1,i) = atomicstress(1,i) - 0.5_dp*d1*x*x
              endif
            endif
            if (lDoQDeriv1) then
              d0i = qj*rr*angstoev
              call d1charged(ii,i,lopi,1_i4,d0i)
            endif
!
            if (lgrad2) then
              d2 = - 3.0_dp*d1*rr*rr
              derv2(jx,ix) = derv2(jx,ix) + d2*x*x
              derv2(jy,ix) = derv2(jy,ix) + d2*y*x
              derv2(jz,ix) = derv2(jz,ix) + d2*z*x
              derv2(jx,iy) = derv2(jx,iy) + d2*x*y
              derv2(jy,iy) = derv2(jy,iy) + d2*y*y
              derv2(jz,iy) = derv2(jz,iy) + d2*z*y
              derv2(jx,iz) = derv2(jx,iz) + d2*x*z
              derv2(jy,iz) = derv2(jy,iz) + d2*y*z
              derv2(jz,iz) = derv2(jz,iz) + d2*z*z
              derv2(jx,ix) = derv2(jx,ix) + d1
              derv2(jy,iy) = derv2(jy,iy) + d1
              derv2(jz,iz) = derv2(jz,iz) + d1
              if (lstr) then
                derv3(ix,1) = derv3(ix,1) + d2*x*x*x
                derv3(iy,1) = derv3(iy,1) + d2*x*x*y
                derv3(iz,1) = derv3(iz,1) + d2*x*x*z
                sderv2(1,1) = sderv2(1,1) - 0.5_dp*d2*x*x*x*x
              endif
!*******************************************************
!  Variable charge contribution to second derivatives  *
!*******************************************************
              if (lDoQDeriv2) then
                d0i  = qj*rr*angstoev
                d0j  = qi*rr*angstoev
                d1i  = - d0i*rr*rr
                d1j  = - d0j*rr*rr
                d2i2 = 0.0_dp
                d2ij = rr*angstoev
                d2j2 = 0.0_dp
                d1ix = d1i*x
                d1iy = d1i*y
                d1iz = d1i*z
                d1jx = d1j*x
                d1jy = d1j*y
                d1jz = d1j*z
                d1is = d1i*x*x
                d1js = d1j*x*x
                call d2charge(i,j,1_i4,ix,iy,iz,jx,jy,jz,lopi,lopj,d0i,d0j,d1ix,d1iy,d1iz, &
                              d1jx,d1jy,d1jz,d1is,d1js,d2i2,d2ij,d2j2,0.0_dp,0.0_dp,0.0_dp, &
                              .true.,.false.)
              endif
            endif
          endif
        endif
        if (m.gt.0) then
          x = - acell + xclat(j) - xclat(i)
          r = x*x + y*y + z*z
          if (r.gt.rcut) then
            r = sqrt(r)
            rr = 1.0_dp/r
            d0 = qi*qj*rr
            if (lreg2one) then
              esum12 = esum12 + d0
            elseif (lreg2pair) then
              esum2 = esum2 + d0
            else
              esum = esum + d0
            endif
!
            eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + 0.5_dp*d0*angstoev
!
            siteenergy(i) = siteenergy(i) + 0.5_dp*d0*angstoev
!
            if (lgrad1) then
              d1 = d0*rr*rr*angstoev
              xdrv(i) = xdrv(i) + d1*x
              ydrv(i) = ydrv(i) + d1*y
              zdrv(i) = zdrv(i) + d1*z
              if (nregioni.ne.nregionj) then
                xregdrv(nregioni) = xregdrv(nregioni) + d1*x
                yregdrv(nregioni) = yregdrv(nregioni) + d1*y
                zregdrv(nregioni) = zregdrv(nregioni) + d1*z
              endif
              if (lstr) then
                rstrd(1) = rstrd(1) - 0.5_dp*d1*x*x
                if (latomicstress) then
                  atomicstress(1,i) = atomicstress(1,i) - 0.5_dp*d1*x*x
                endif
              endif
              if (lDoQDeriv1) then
                d0i = qj*rr*angstoev
                call d1charged(ii,i,lopi,1_i4,d0i)
              endif
              if (lgrad2) then
                d2 = - 3.0_dp*d1*rr*rr
                derv2(jx,ix) = derv2(jx,ix) + d2*x*x
                derv2(jy,ix) = derv2(jy,ix) + d2*y*x
                derv2(jz,ix) = derv2(jz,ix) + d2*z*x
                derv2(jx,iy) = derv2(jx,iy) + d2*x*y
                derv2(jy,iy) = derv2(jy,iy) + d2*y*y
                derv2(jz,iy) = derv2(jz,iy) + d2*z*y
                derv2(jx,iz) = derv2(jx,iz) + d2*x*z
                derv2(jy,iz) = derv2(jy,iz) + d2*y*z
                derv2(jz,iz) = derv2(jz,iz) + d2*z*z
                derv2(jx,ix) = derv2(jx,ix) + d1
                derv2(jy,iy) = derv2(jy,iy) + d1
                derv2(jz,iz) = derv2(jz,iz) + d1
                if (lstr) then
                  derv3(ix,1) = derv3(ix,1) + d2*x*x*x
                  derv3(iy,1) = derv3(iy,1) + d2*x*x*y
                  derv3(iz,1) = derv3(iz,1) + d2*x*x*z
                  sderv2(1,1) = sderv2(1,1) - 0.5_dp*d2*x*x*x*x
                endif
!*******************************************************
!  Variable charge contribution to second derivatives  *
!*******************************************************
                if (lDoQDeriv2) then
                  d0i  = qj*rr*angstoev
                  d0j  = qi*rr*angstoev
                  d1i  = - d0i*rr*rr
                  d1j  = - d0j*rr*rr
                  d2i2 = 0.0_dp
                  d2ij = rr*angstoev
                  d2j2 = 0.0_dp
                  d1ix = d1i*x
                  d1iy = d1i*y
                  d1iz = d1i*z
                  d1jx = d1j*x
                  d1jy = d1j*y
                  d1jz = d1j*z
                  d1is = d1i*x*x
                  d1js = d1j*x*x
                  call d2charge(i,j,1_i4,ix,iy,iz,jx,jy,jz,lopi,lopj,d0i,d0j,d1ix,d1iy,d1iz, &
                                d1jx,d1jy,d1jz,d1is,d1js,d2i2,d2ij,d2j2,0.0_dp,0.0_dp,0.0_dp, &
                                .true.,.false.)
                endif
              endif
            endif
          endif
        endif
      enddo jloop
    enddo
!
!  Neutralising terms
!
!  Background
!
!  and
!
!  Euler-MacLaurin component
!
    esumh = 0.0_dp
    esumh12 = 0.0_dp
    esumh2 = 0.0_dp
    esumem = 0.0_dp
    esumem12 = 0.0_dp
    esumem2 = 0.0_dp
!
    if (m.gt.0) then
      u = (dble(m)+0.5_dp)*a
      do ii = 1,natomsonnode
        i = node2atom(ii)
        oci = occuf(i)
        qi = qf(i)*oci
        nregioni = nregionno(nsft+i)
        nregiontypi = nregiontype(nregioni,ncf)
        lopi = .true.
        jloop2: do j = 1,numat
          ocj = occuf(j)
          qj = qf(j)*ocj
          nregionj = nregionno(nsft+j)
          nregiontypj = nregiontype(nregionj,ncf)
          lopj = .true.
!     
!  Set region 2 pair flags
!     
          lreg2one  = .false.
          lreg2pair = .false.
          if (lseok.and.nregions(ncf).ge.2) then
            lreg2pair = (nregioni.eq.2.and.nregionj.eq.2)
            if (.not.lreg2pair) lreg2one = (nregioni.eq.2.or.nregionj.eq.2)
          endif
!
!  QM/MM handling : i & j are QM atoms => exclude
!
          if (QMMMmode(ncf).gt.0) then
            if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop2
!
!  QM/MM : Electrostatic embedding : If either i or j are QM atoms => exclude 
!
            if (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1)) cycle jloop2
          endif
!
          qij = qi*qj
!
          x = xclat(j) - xclat(i)
          y = yclat(j) - yclat(i)
          z = zclat(j) - zclat(i)
          call hfunc(u,+x,1.0_dp,y,z,h1,dh1,d2h1,d3h1,.false.,.false.,.false.)
          call hfunc(u,-x,-1.0_dp,y,z,h2,dh2,d2h2,d3h1,.false.,.false.,.false.)
          if (lreg2one) then
            esumh12 = esumh12 - qij*(h1 + h2 - 2.0_dp*lna)/a
          elseif (lreg2pair) then
            esumh2 = esumh2 - qij*(h1 + h2 - 2.0_dp*lna)/a
          else
            esumh = esumh - qij*(h1 + h2 - 2.0_dp*lna)/a
          endif
!
          call emfunc(nemorder,u,+x,1.0_dp,y,z,a,e1,dh1,d2h1,dh1s,d2h1s,d2h1m,d3h1,d3h1m,.false.,.false.,.false.)
          call emfunc(nemorder,u,-x,-1.0_dp,y,z,a,e2,dh2,d2h2,dh2s,d2h2s,d2h2m,d3h1,d3h1m,.false.,.false.,.false.)
          if (lreg2one) then
            esumem12 = esumem12 + qij*(e1 + e2)
          elseif (lreg2pair) then
            esumem2 = esumem2 + qij*(e1 + e2)
          else
            esumem = esumem + qij*(e1 + e2)
          endif
        enddo jloop2
      enddo
    endif
!
!  Global sum of ereal to ensure consistent convergence
!
    call sumall(esum+esumh+esumem,ereal,1_i4,"real1Dd","ereal")
!
!  Correct for double counting
!
    ereal = 0.5_dp*ereal
!
!  Compare to energy with previous number of cells
!  and check for convergence to the required accuracy.
!
    if (abs(ereal).lt.accf) lconverged = .true.
    if (.not.lconverged) then
      ediff = abs((ereal - elast)/ereal)
      lconverged = (ediff.lt.accf)
    endif
    if (index(keyword,'verb').ne.0) then
      if (ioproc) then
        write(ioout,'('' Mcell = '',i4,''  Ereal(1-D) = '',f15.6,'' eV'')') m,ereal*angstoev
      endif
    endif
    elast = ereal
  enddo
  if (index(keyword,'verb').ne.0.and.ioproc) write(ioout,'(/)')
!
!  Divide ereal by number of processors to compensate for the fact that it has already been summed
!
  ereal = ereal*angstoev/dble(nprocs)
  erealin = erealin + ereal
  esregion12 = esregion12 + (esum12 + esumh12 + esumem12)*angstoev
  esregion2 = esregion2 + (esum2 + esumh2 + esumem2)*angstoev
!
!  Save number of cells needed
!
  maxloop(1) = m
!
!  Derivatives of Euler-MacLaurin terms since these are not cumulative
!
  if (lgrad1.and.m.gt.0) then
    u = (dble(m)+0.5_dp)*a
    ix = -2
    iy = -1
    iz =  0
    do ii = 1,natomsonnode
      i = node2atom(ii)
      oci = occuf(i)
      qi = qf(i)*oci
      nregioni = nregionno(nsft+i)
      nregiontypi = nregiontype(nregioni,ncf)
      lopi = .true.
      ix = ix + 3
      iy = iy + 3
      iz = iz + 3
      jx = - 2
      jy = - 1
      jz =   0
      jloop3: do j = 1,numat
        ocj = occuf(j)
        qj = qf(j)*ocj
        nregionj = nregionno(nsft+j)
        nregiontypj = nregiontype(nregionj,ncf)
!
!  QM/MM handling : i & j are QM atoms => exclude
!
        if (QMMMmode(ncf).gt.0) then
          if (nregiontypi.eq.1.and.nregiontypj.eq.1) cycle jloop3
!
!  QM/MM : Electrostatic embedding : If either i or j are QM atoms => exclude 
!
          if (QMMMmode(ncf).eq.2.and.(nregiontypi.eq.1.or.nregiontypj.eq.1)) cycle jloop3
        endif
        lopj = .true.
        jx = jx + 3
        jy = jy + 3
        jz = jz + 3  
!
        qij = qi*qj
        x = xclat(j) - xclat(i)
        y = yclat(j) - yclat(i)
        z = zclat(j) - zclat(i)
!
        call hfunc(u,+x,1.0_dp,y,z,h1,dh1,d2h1,d3h1,lgrad1,lgrad2,.false.)
        call hfunc(u,-x,-1.0_dp,y,z,h2,dh2,d2h2,d3h1,lgrad1,lgrad2,.false.)
!
        etrm = qij*angstoev*(h1 + h2 - 2.0_dp*lna)/a
        eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) - 0.5_dp*etrm
!
        siteenergy(i) = siteenergy(i) - 0.5_dp*etrm
!
        d1 = qij*angstoev/a
        xdrv(i) = xdrv(i) + d1*(dh1(1) + dh2(1))
        ydrv(i) = ydrv(i) + d1*(dh1(2) + dh2(2))
        zdrv(i) = zdrv(i) + d1*(dh1(3) + dh2(3))
        if (nregioni.ne.nregionj) then
          xregdrv(nregioni) = xregdrv(nregioni) + d1*(dh1(1) + dh2(1))
          yregdrv(nregioni) = yregdrv(nregioni) + d1*(dh1(2) + dh2(2))
          zregdrv(nregioni) = zregdrv(nregioni) + d1*(dh1(3) + dh2(3))
        endif
        if (lstr) then
          d1s = - (dh1(1)*(u+x)-dh2(1)*(u-x))
          d1s = d1s + (h1+h2)
          d1s = d1s + 2.0_dp*(1.0d0-lna)
          d1s = d1s*angstoev/a
          d1snoq = d1s
          d1s = 0.5_dp*d1s*qij
          rstrd(1) = rstrd(1) + d1s
          if (latomicstress) then
            atomicstress(1,i) = atomicstress(1,i) + d1s
          endif
        endif
!
!  Compute variable charge term = E without charges
!
        d0term = - angstoev*(h1 + h2 - 2.0_dp*lna)/a
        if (lgrad2) then
!
!  Second derivatives
!
          d2 = qij*angstoev/a
          derv2(jx,ix) = derv2(jx,ix) + d2*(d2h1(1)+d2h2(1))
          derv2(jy,ix) = derv2(jy,ix) + d2*(d2h1(6)+d2h2(6))
          derv2(jz,ix) = derv2(jz,ix) + d2*(d2h1(5)+d2h2(5))
          derv2(jx,iy) = derv2(jx,iy) + d2*(d2h1(6)+d2h2(6))
          derv2(jy,iy) = derv2(jy,iy) + d2*(d2h1(2)+d2h2(2))
          derv2(jz,iy) = derv2(jz,iy) + d2*(d2h1(4)+d2h2(4))
          derv2(jx,iz) = derv2(jx,iz) + d2*(d2h1(5)+d2h2(5))
          derv2(jy,iz) = derv2(jy,iz) + d2*(d2h1(4)+d2h2(4))
          derv2(jz,iz) = derv2(jz,iz) + d2*(d2h1(3)+d2h2(3))
          if (lstr) then
            derv3(ix,1) = derv3(ix,1) + d2*(d2h1(1)*(u+x) - d2h2(1)*(u-x))  &
                                      - d1*(dh1(1) + dh2(1)) - d1*(dh1(1) + dh2(1))
            derv3(iy,1) = derv3(iy,1) + d2*(d2h1(6)*(u+x) - d2h2(6)*(u-x))  &
                                      - d1*(dh1(2) + dh2(2))
            derv3(iz,1) = derv3(iz,1) + d2*(d2h1(5)*(u+x) - d2h2(5)*(u-x))  &
                                      - d1*(dh1(3) + dh2(3))
            sderv2(1,1) = sderv2(1,1) - 0.5_dp*d2*(d2h1(1)*(u+x)*(u+x) + d2h2(1)*(u-x)*(u-x))  &
                                      - d2 - 1.5_dp*d1s
          endif
          if (lDoQDeriv2) then
!
!  Accumulate terms that contribute to the variable charge second derivatives
!
            d1ix = - qj*angstoev*(dh1(1) + dh2(1))/a
            d1iy = - qj*angstoev*(dh1(2) + dh2(2))/a
            d1iz = - qj*angstoev*(dh1(3) + dh2(3))/a
            d1jx = - qi*angstoev*(dh1(1) + dh2(1))/a
            d1jy = - qi*angstoev*(dh1(2) + dh2(2))/a
            d1jz = - qi*angstoev*(dh1(3) + dh2(3))/a
            if (lstr) then
              d1is = d1snoq*qj
              d1js = d1snoq*qi
            else
              d1is = 0.0_dp
              d1js = 0.0_dp
            endif
          endif
        endif
!
        call emfunc(nemorder,u,+x,1.0_dp,y,z,a,e1,dh1,d2h1,dh1s, &
                    d2h1s,d2h1m,d3h1,d3h1m,lgrad1,lgrad2,.false.)
        call emfunc(nemorder,u,-x,-1.0_dp,y,z,a,e2,dh2,d2h2,dh2s, &
                    d2h2s,d2h2m,d3h1,d3h1m,lgrad1,lgrad2,.false.)
!
        etrm = qij*angstoev*(e1 + e2)
        eregion2region(nregionj,nregioni) = eregion2region(nregionj,nregioni) + 0.5_dp*etrm
!
        siteenergy(i) = siteenergy(i) + 0.5_dp*etrm
!
        d1 = qij*angstoev
        xdrv(i) = xdrv(i) - d1*(dh1(1) + dh2(1))
        ydrv(i) = ydrv(i) - d1*(dh1(2) + dh2(2))
        zdrv(i) = zdrv(i) - d1*(dh1(3) + dh2(3))
        if (nregioni.ne.nregionj) then
          xregdrv(nregioni) = xregdrv(nregioni) - d1*(dh1(1) + dh2(1))
          yregdrv(nregioni) = yregdrv(nregioni) - d1*(dh1(2) + dh2(2))
          zregdrv(nregioni) = zregdrv(nregioni) - d1*(dh1(3) + dh2(3))
        endif
        if (lstr) then
          rstrd(1) = rstrd(1) + 0.5_dp*d1*(dh1s + dh2s)
          if (latomicstress) then
            atomicstress(1,i) = atomicstress(1,i) + 0.5_dp*d1*(dh1s + dh2s)
          endif
        endif
        d0term = d0term + (e1 + e2)*angstoev 
        if (lDoQDeriv1) then
          d0i = qj*d0term
          call d1charged(ii,i,lopi,1_i4,d0i)
        endif
        if (lgrad2) then
          d2 = qij*angstoev
          derv2(jx,ix) = derv2(jx,ix) - d2*(d2h1(1)+d2h2(1))
          derv2(jy,ix) = derv2(jy,ix) - d2*(d2h1(6)+d2h2(6))
          derv2(jz,ix) = derv2(jz,ix) - d2*(d2h1(5)+d2h2(5))
          derv2(jx,iy) = derv2(jx,iy) - d2*(d2h1(6)+d2h2(6))
          derv2(jy,iy) = derv2(jy,iy) - d2*(d2h1(2)+d2h2(2))
          derv2(jz,iy) = derv2(jz,iy) - d2*(d2h1(4)+d2h2(4))
          derv2(jx,iz) = derv2(jx,iz) - d2*(d2h1(5)+d2h2(5))
          derv2(jy,iz) = derv2(jy,iz) - d2*(d2h1(4)+d2h2(4))
          derv2(jz,iz) = derv2(jz,iz) - d2*(d2h1(3)+d2h2(3))
          if (lstr) then
            derv3(ix,1) = derv3(ix,1) + d2*(d2h1m(1)+d2h2m(1))
            derv3(iy,1) = derv3(iy,1) + d2*(d2h1m(2)+d2h2m(2))
            derv3(iz,1) = derv3(iz,1) + d2*(d2h1m(3)+d2h2m(3))
            sderv2(1,1) = sderv2(1,1) + 0.5_dp*d2*(d2h1s + d2h2s)
          endif
!*******************************************************
!  Variable charge contribution to second derivatives  *
!*******************************************************
          if (lDoQDeriv2) then
            d0i  = qj*d0term
            d0j  = qi*d0term
            d2i2 = 0.0_dp
            d2ij = d0term
            d2j2 = 0.0_dp
            d1ix = d1ix + qj*angstoev*(dh1(1) + dh2(1))
            d1iy = d1iy + qj*angstoev*(dh1(2) + dh2(2))
            d1iz = d1iz + qj*angstoev*(dh1(3) + dh2(3))
            d1jx = d1jx + qi*angstoev*(dh1(1) + dh2(1))
            d1jy = d1jy + qi*angstoev*(dh1(2) + dh2(2))
            d1jz = d1jz + qi*angstoev*(dh1(3) + dh2(3))
            if (lstr) then
              d1is = d1is + (dh1s + dh2s)*qj*angstoev
              d1js = d1js + (dh1s + dh2s)*qi*angstoev
            endif
            call d2charge(i,j,1_i4,ix,iy,iz,jx,jy,jz,lopi,lopj,d0i,d0j,d1ix,d1iy,d1iz, &
                          d1jx,d1jy,d1jz,d1is,d1js,d2i2,d2ij,d2j2,0.0_dp,0.0_dp,0.0_dp, &
                          .true.,.false.)
          endif
        endif
      enddo jloop3
    enddo
  endif
!
  t2 = g_cpu_time()
  tatom = tatom + t2 - t1
!
  return
  end
