  subroutine changemaxkimmodel
!
!  Alters the size of the arrays associated with maxkimmodel
!
!  10/12 Created
!   7/14 Modified for changes to OpenKIM v1.4
!   7/16 Modified for changes to OpenKIM v1.7.3
!   1/18 Reallocation of pkim_model removed for F03 version
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
!  Copyright Curtin University 2018
!
!  Julian Gale, CIC, Curtin University, January 2018
!
  use configurations, only : maxcfg
  use kim_models
  use reallocate
#if KIM_F03
  use, intrinsic :: iso_c_binding
  use KIM_API_F03
  use reallocate_kim
#include "KIM_API_status.h"
#elif KIM
  use KIM_API
  use reallocate_kim
#endif
  implicit none
!
#if KIM_F03
  integer(c_int)    :: ierror
#elif KIM
  integer(i4)       :: ierror
#else
  integer(i4)       :: ierror
#endif
  integer(i4)       :: i, j
  integer(i4), save :: oldmaxkimmodel = 0
!
#if KIM_F03
  call realloc_ch80(kim_model,maxkimmodel,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_model')
  call realloc(kim_cutoff,maxkimmodel,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_cutoff')
  call realloc(kim_nbc,maxkimmodel,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_nbc')
#elif KIM
  call realloc_ch80(kim_model,maxkimmodel,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_model')
  call realloc(kim_cutoff,maxkimmodel,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_cutoff')
  call realloc(kim_nbc,maxkimmodel,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_nbc')
#endif
  call realloc(lkim_model_cfg_OK,maxkimmodel,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','lkim_model_cfg_OK')
#if KIM_F03
  call realloc_kim_len(kim_test_descriptor,maxkimmodel,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_test_descriptor')
#elif KIM
  call realloc_kim_len(kim_test_descriptor,maxkimmodel,maxcfg,ierror)
  if (ierror.ne.0) call outofmemory('changemaxkimmodel','kim_test_descriptor')
#endif
!
!  Initialise new part of data array
!
  if (maxkimmodel.gt.oldmaxkimmodel) then
    do i = 1,maxcfg
      do j = oldmaxkimmodel+1,maxkimmodel
        lkim_model_cfg_OK(j,i) = .true.
#if KIM_F03
        kim_nbc(j,i) = 0
#elif KIM
        kim_nbc(j,i) = 0
        pkim_model(j,i) = 0
#endif
      enddo
    enddo
#if KIM_F03
    do i = oldmaxkimmodel+1,maxkimmodel
      kim_model(i) = ' '
      kim_cutoff(i) = 0.0_dp
    enddo
#elif KIM
    do i = oldmaxkimmodel+1,maxkimmodel
      kim_model(i) = ' '
      kim_cutoff(i) = 0.0_dp
    enddo
#endif
  endif
!
!  Save current value of maxkimmodel for next call
!
  oldmaxkimmodel = maxkimmodel
!
  return
  end
