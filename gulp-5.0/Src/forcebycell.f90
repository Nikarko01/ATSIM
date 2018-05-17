  subroutine forcebycell(maxlim,d1cell,matom,vmatom)
!
!  Routine for generating forces stored as separate images.
!  These are used for finite differences to obtain force constants
!  for later phasing in phonon calculations.
!
!  11/14 Created from dynamic_fc
!  12/14 Displaced atom coordinates/radius passed in separately
!   1/15 Bondorder and brenner potentials added
!   2/15 Modified so that density1fc is called for EAM case since
!        real1fc can't provide density with current algorithm
!   2/15 Criterion for calling getBOcharge changed by adding nboQ0
!  10/15 ReaxFF added
!  11/17 EDIP added
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
!  Julian Gale, CIC, Curtin University, November 2017
!
  use bondorderdata, only : nbopot, nboQ, nboQ0
  use control
  use current
  use derivatives
  use eam,           only : lMEAM, maxmeamcomponent
  use four
  use six
  use sutton
  use m_bondorder,   only : bondorder1fc2
  use m_edip,        only : EDIP1fc2
  use m_three
  implicit none
!
!  Passed variables
!
  integer(i4),                  intent(in)  :: maxlim
  integer(i4),                  intent(in)  :: matom
  real(dp),                     intent(out) :: d1cell(4,maxlim,*)
  real(dp),                     intent(in)  :: vmatom(4)
!
!  Local variables
!
  integer(i4)                               :: i
!
!  Zero derivatives
!
  d1cell(1:4,1:maxlim,1:nd2cells) = 0.0_dp
!
!  Zero many-body terms since finite difference will change the densities
!
  if (lsuttonc) then
    if (lMEAM) then
      do i = 1,numat
        scrho(1:maxmeamcomponent,i) = 0.0_dp
      enddo
    else
      do i = 1,numat
        scrho(1,i) = 0.0_dp
      enddo
    endif
  endif
!**********************************************
!  EEM/QEq calculation of charge derivatives  *
!**********************************************
  if (leem) then
    call dcharge(.false.,.true.,.false.)
  endif
!*************************************************
!  Bond Order calculation of charge derivatives  *
!*************************************************
  if ((nboQ+nboQ0).gt.0) then
    call getBOcharge(.true.,.false.)
  endif
!*************************
!  Real space component  *
!*************************
  call real1fc(maxlim,d1cell,matom,vmatom)
!*************************
!  Three-body component  *
!*************************
  if (nthb.gt.0) call three1fc(maxlim,d1cell,matom,vmatom)
!************************
!  Four-body component  *
!************************
  if (nfor.gt.0) call four1fc(maxlim,d1cell,matom,vmatom)
!***********************
!  Six-body component  *
!***********************
  if (nsix.gt.0) call six1fc(maxlim,d1cell,matom,vmatom)
!************************
!  Many-body component  *
!************************
  if (lsuttonc) then
    call density1fc(matom,vmatom)
    call many1fc(maxlim,d1cell,matom,vmatom)
  endif
!**********************
!  Brenner potential  *
!**********************
  if (lbrenner) then
    call brenner1fc(maxlim,d1cell,matom,vmatom)
  endif
!************************
!  Bondorder potential  *
!************************
  if (nbopot.gt.0) then
    if (lfastfd) then
      call bondorder1fc2(maxlim,d1cell,matom,vmatom)
    else
      call bondorder1fc(maxlim,d1cell,matom,vmatom)
    endif
  endif
!********************
!  EDIP forcefield  *
!********************
  if (lEDIP) then
    if (lfastfd) then
      call EDIP1fc2(maxlim,d1cell,matom,vmatom)
    else
      call EDIP1fc(maxlim,d1cell,matom,vmatom)
    endif
  endif
!**********************
!  ReaxFF forcefield  *
!**********************
  if (lreaxFF) then
    call reaxff1fc(maxlim,d1cell,matom,vmatom)
  endif
  return
  end
