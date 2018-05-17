  subroutine kimmd(nkim,ekim,lgrad1)
!
!  Calculates the energy and up to first derivatives for OpenKIM potentials.
!
!  On entry : 
!
!  nkim            = integer reference to number of KIM model for this call
!  lgrad1          = if .true. calculate the first derivatives
!
!  On exit :
!
!  ekim            = the value of the energy contribution from the nkim'th model
!
!  10/12 Created from bondordermd.f90
!   4/14 Modifications for OpenKIM version 1.4 made
!   7/14 Modifications for OpenKIM version 1.7.3 made
!   8/16 Forces now subtracted to change sign to be derivative
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
!  Copyright Curtin University 2016
!
!  Julian Gale, CIC, Curtin University, August 2016
!
#ifdef KIM_F03
#include "KIM_API_status.h"
  use, intrinsic :: iso_c_binding
#elif KIM
#include "KIM_API_status.h"
#endif
  use datatypes
  use current
  use iochannels
#if KIM_F03
  use derivatives,    only : xdrv, ydrv, zdrv, rstrd
  use kim_models,     only : kim_energy, kim_coord, kim_forces, kim_virial
  use kim_models,     only : pkim_model
  use kim_models,     only : kim_nbc
  use symmetry,       only : lstr
#elif KIM
  use derivatives,    only : xdrv, ydrv, zdrv, rstrd
  use kim_models,     only : kim_energy, kim_coord, kim_forces, kim_virial
  use kim_models,     only : pkim_model
  use kim_models,     only : kim_nbc
  use symmetry,       only : lstr
#endif
  use kim_functions,  only : set_kim_neighbours
  use parallel
  use times
#ifdef KIM_F03
  use KIM_API_F03
#elif KIM
  use KIM_API
#endif
  implicit none
!
!  Passed variables
!
  integer(i4),   intent(in)                        :: nkim               ! KIM model number for this call
  real(dp),      intent(out)                       :: ekim               ! KIM model energy on return
  logical,       intent(in)                        :: lgrad1             ! If true, then compute the forces
!
!  Local variables
!
#ifdef KIM_F03
  integer(i4)                                      :: i
  integer(c_int)                                   :: KIMerror           ! KIM error flag
  integer(c_int)                                   :: KIMerror_d         ! KIM error flag
#elif KIM
  integer(i4)                                      :: i
  integer(i4)                                      :: KIMerror           ! KIM error flag
  integer(i4)                                      :: KIMerror_d         ! KIM error flag
#endif
  real(dp)                                         :: g_cpu_time
  real(dp)                                         :: t1
  real(dp)                                         :: t2
!
  t1 = g_cpu_time()
!
#if KIM_F03
  if (kim_nbc(nkim,ncf).eq.1) then
!
!  Call neighbour list set up for KIM potentials
!
    call set_kim_neighbours(nkim)
  endif
#elif KIM
  if (kim_nbc(nkim,ncf).eq.1) then
!
!  Call neighbour list set up for KIM potentials
!
    call set_kim_neighbours(nkim)
  endif
#endif
!
!  Initialise energy
!
  ekim = 0.0_dp

#ifdef KIM_F03
!*********************************************
!  Compute KIM energy and optionally forces  *
!*********************************************
!
!  Copy coordinates of configuration to KIM arrays
!
  do i = 1,numat
    kim_coord(1,i) = xclat(i)
    kim_coord(2,i) = yclat(i)
    kim_coord(3,i) = zclat(i)
  enddo
!
!  Set flags to KIM that tell it whether the energy/forces are needed
!
  call kim_api_set_compute(pkim_model(nkim,ncf),"energy",KIM_COMPUTE_TRUE,KIMerror)
  if (KIMerror.lt.KIM_STATUS_OK) then
    if (ioproc) then
      KIMerror_d = kim_api_report_error(__LINE__,"kimmd","kim_api_set_compute",KIMerror)
    endif
    call stopnow('kimmd')
  endif
  if (lgrad1) then
    call kim_api_set_compute(pkim_model(nkim,ncf),"forces",KIM_COMPUTE_TRUE,KIMerror)
  else
    call kim_api_set_compute(pkim_model(nkim,ncf),"forces",KIM_COMPUTE_FALSE,KIMerror)
  endif
  if (KIMerror.lt.KIM_STATUS_OK) then
    if (ioproc) then
      KIMerror_d = kim_api_report_error(__LINE__,"kimmd","kim_api_set_compute",KIMerror)
    endif
    call stopnow('kimmd')
  endif
  if (ndim.gt.0) then
    if (lstr) then
      call kim_api_set_compute(pkim_model(nkim,ncf),"virial",KIM_COMPUTE_TRUE,KIMerror)
    else
      call kim_api_set_compute(pkim_model(nkim,ncf),"virial",KIM_COMPUTE_FALSE,KIMerror)
    endif
    if (KIMerror.lt.KIM_STATUS_OK) then
      if (ioproc) then
        KIMerror_d = kim_api_report_error(__LINE__,"kimmd","kim_api_set_compute",KIMerror)
      endif
      call stopnow('kimmd')
    endif
  endif
!
!  Call KIM model compute
!
  KIMerror = kim_api_model_compute(pkim_model(nkim,ncf))
  if (KIMerror.lt.KIM_STATUS_OK) then
    if (ioproc) then
      KIMerror_d = kim_api_report_error(__LINE__,"kimmd","kim_api_model_compute",KIMerror)
    endif
    call stopnow('kimmd')
  endif
!
!  Copy energy, forces and virial back
!
  ekim = ekim + kim_energy
  if (lgrad1) then
    do i = 1,numat
      xdrv(i) = xdrv(i) - kim_forces(1,i)
      ydrv(i) = ydrv(i) - kim_forces(2,i)
      zdrv(i) = zdrv(i) - kim_forces(3,i)
    enddo
    if (lstr) then
      select case(ndim)
        case(1)
          rstrd(1) = rstrd(1) + kim_virial(1)
        case(2)
          rstrd(1) = rstrd(1) + kim_virial(1)
          rstrd(2) = rstrd(2) + kim_virial(2)
          rstrd(3) = rstrd(3) + kim_virial(3)
        case(3)
          rstrd(1) = rstrd(1) + kim_virial(1)
          rstrd(2) = rstrd(2) + kim_virial(2)
          rstrd(3) = rstrd(3) + kim_virial(3)
          rstrd(4) = rstrd(4) + kim_virial(4)
          rstrd(5) = rstrd(5) + kim_virial(5)
          rstrd(6) = rstrd(6) + kim_virial(6)
      end select
    endif
  endif
!
!  NB: Need to handle regions & eattach!
!
#elif KIM
!*********************************************
!  Compute KIM energy and optionally forces  *
!*********************************************
!
!  Copy coordinates of configuration to KIM arrays
!
  do i = 1,numat
    kim_coord(1,i) = xclat(i)
    kim_coord(2,i) = yclat(i)
    kim_coord(3,i) = zclat(i)
  enddo
!
!  Set flags to KIM that tell it whether the energy/forces are needed
!
  call kim_api_set_compute_f(pkim_model(nkim,ncf),"energy",KIM_COMPUTE_TRUE,KIMerror)
  if (KIMerror.lt.KIM_STATUS_OK) then
    if (ioproc) then
      KIMerror_d = kim_api_report_error_f(__LINE__,"kimmd","kim_api_set_compute_f",KIMerror)
    endif
    call stopnow('kimmd')
  endif
  if (lgrad1) then
    call kim_api_set_compute_f(pkim_model(nkim,ncf),"forces",KIM_COMPUTE_TRUE,KIMerror)
  else
    call kim_api_set_compute_f(pkim_model(nkim,ncf),"forces",KIM_COMPUTE_FALSE,KIMerror)
  endif
  if (KIMerror.lt.KIM_STATUS_OK) then
    if (ioproc) then
      KIMerror_d = kim_api_report_error_f(__LINE__,"kimmd","kim_api_set_compute_f",KIMerror)
    endif
    call stopnow('kimmd')
  endif
  if (ndim.gt.0) then
    if (lstr) then
      call kim_api_set_compute_f(pkim_model(nkim,ncf),"virial",KIM_COMPUTE_TRUE,KIMerror)
    else
      call kim_api_set_compute_f(pkim_model(nkim,ncf),"virial",KIM_COMPUTE_FALSE,KIMerror)
    endif
    if (KIMerror.lt.KIM_STATUS_OK) then
      if (ioproc) then
        KIMerror_d = kim_api_report_error_f(__LINE__,"kimmd","kim_api_set_compute_f",KIMerror)
      endif
      call stopnow('kimmd')
    endif
  endif
!
!  Call KIM model compute
!
  KIMerror = kim_api_model_compute_f(pkim_model(nkim,ncf))
  if (KIMerror.lt.KIM_STATUS_OK) then
    if (ioproc) then
      KIMerror_d = kim_api_report_error_f(__LINE__,"kimmd","kim_api_model_compute_f",KIMerror)
    endif
    call stopnow('kimmd')
  endif
!
!  Copy energy, forces and virial back
!
  ekim = ekim + kim_energy
  if (lgrad1) then
    do i = 1,numat
      xdrv(i) = xdrv(i) - kim_forces(1,i)
      ydrv(i) = ydrv(i) - kim_forces(2,i)
      zdrv(i) = zdrv(i) - kim_forces(3,i)
    enddo
    if (lstr) then
      select case(ndim)
        case(1)
          rstrd(1) = rstrd(1) + kim_virial(1)
        case(2)
          rstrd(1) = rstrd(1) + kim_virial(1)
          rstrd(2) = rstrd(2) + kim_virial(2)
          rstrd(3) = rstrd(3) + kim_virial(3)
        case(3)
          rstrd(1) = rstrd(1) + kim_virial(1)
          rstrd(2) = rstrd(2) + kim_virial(2)
          rstrd(3) = rstrd(3) + kim_virial(3)
          rstrd(4) = rstrd(4) + kim_virial(4)
          rstrd(5) = rstrd(5) + kim_virial(5)
          rstrd(6) = rstrd(6) + kim_virial(6)
      end select
    endif
  endif
!
!  NB: Need to handle regions & eattach!
!
#endif
  t2 = g_cpu_time()
  tkim = tkim + t2 - t1
!
  return
  end
