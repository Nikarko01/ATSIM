  subroutine secondorderfcd(lshelleliminate)
!
!  Calculates the second order force constant matrix for the cores
!  Distributed memory parallel version.
!
!  The block diagonal elements of derv2 must be stored, as these are
!  not recalculated in dynamic as the summing of off diagonals is no
!  longer applicable.
!
!   6/17 Created from secondorderfc
!   8/17 Routine wrapped in ifdef since it won't work or be called without MPI
!
!  NB: Partial occupancy and breathing shells not currently supported
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
!  Julian Gale, CIC, Curtin University, August 2017
!
  use configurations
  use g_constants
  use control
  use current
  use derivatives
  use element
  use gulp_files
  use frequencies
  use general
  use iochannels
  use parallel
  use partial
  use phononatoms
  use properties
  use shells
  use times
  implicit none
!
!  Passed variables
!
  logical,                           intent(in)     :: lshelleliminate
#ifdef MPI
!
!  Local variables
!
  integer(i4)                                       :: i
  integer(i4)                                       :: ifail
  integer(i4)                                       :: iloc
  integer(i4)                                       :: indi
  integer(i4)                                       :: indiloc
  integer(i4)                                       :: indii
  integer(i4)                                       :: indj
  integer(i4)                                       :: indjj
  integer(i4)                                       :: ix
  integer(i4)                                       :: iy
  integer(i4)                                       :: iz
  integer(i4)                                       :: j
  integer(i4)                                       :: m
  integer(i4)                                       :: mm
  integer(i4)                                       :: maxlim
  integer(i4)                                       :: mcv
  integer(i4)                                       :: mcvloc
  integer(i4),  dimension(:),     allocatable       :: mcvrptr
  integer(i4)                                       :: mint
  integer(i4)                                       :: mintloc
  integer(i4)                                       :: msv
  integer(i4)                                       :: msvloc
  integer(i4)                                       :: status
!
  integer                                           :: idesc(9)
  integer                                           :: ifails
  integer                                           :: ld
  integer                                           :: nb
  integer                                           :: ncs
!
  logical                                           :: lnoanald2loc
  real(dp)                                          :: fc
!
!  Set logicals
!
  lnoanald2loc = (lnoanald2)
!************************************************
!  Find number of atoms for phonon calculation  *
!************************************************
!
!  Set up points for atoms involve in phonon calculations
!
  call setphonptr
!
!  Allocate local pointer arrays
!
  lpocc = (nsfoc+ncfoc.ne.nphonat)
  if (lpocc) then
    call outerror('partial occupancy not compatible with distributed 2nd derivs',0_i4)
    call stopnow('secondorderfcd')
  endif
!
!  Calculate a few constants to do with the size of the problem
!
  mint = 3*nphonat
  mintloc = 3*nphonatonnode
!
  maxlim = mint
  if (nbsmat.gt.0) maxlim = maxlim + nphonat
  msv = 3*nsfoc + nbfoc
  mcv = 3*ncfoc
  msvloc = 3*nphonatonnodes + nphonatonnodeb
  mcvloc = 3*nphonatonnodec
!
  allocate(mcvrptr(mcv),stat=status)
  if (status/=0) call outofmemory('secondorderfcd','mcvrptr')
!
  mcvrptr(1:mcv) = 0
  do iloc = 1,nphonatonnodec
    i = nphonatonnodecptr(iloc)
    m = 3*(iloc-1)
    mm = 3*(i-1)
    mcvrptr(mm+1) = m + 1
    mcvrptr(mm+2) = m + 2
    mcvrptr(mm+3) = m + 3
  enddo
!
!  Check that maxd2 is greater than or equal to mcv
!
  if (maxd2.lt.mcv) then
    maxd2 = mcv
    call changemaxd2
  endif
  if (maxd2u.lt.mcvloc) then
    maxd2u = mcvloc
    call changemaxd2
  endif
!
!  Evaluate second derivatives
!
  fc = 0.0_dp
  call energy(fc,.true.,.true.)
!
  if (.not.lnoanald2loc.and..not.lfinitediff) then
!
!  Store diagonal blocks in derv3 to avoid recalculation
!
    do i = 1,nphonat
      iloc = atom2local(nphonatptr(i))
      if (iloc.gt.0) then
!
!  Atom is local to this node
!
        indi = 3*(nphonatptr(i)-1)
        indiloc = 3*(iloc-1)
        derv3(indiloc+1,1) = derv2(indi+1,indiloc+1)
        derv3(indiloc+2,1) = derv2(indi+2,indiloc+1)
        derv3(indiloc+3,1) = derv2(indi+3,indiloc+1)
        derv3(indiloc+1,2) = derv2(indi+1,indiloc+2)
        derv3(indiloc+2,2) = derv2(indi+2,indiloc+2)
        derv3(indiloc+3,2) = derv2(indi+3,indiloc+2)
        derv3(indiloc+1,3) = derv2(indi+1,indiloc+3)
        derv3(indiloc+2,3) = derv2(indi+2,indiloc+3)
        derv3(indiloc+3,3) = derv2(indi+3,indiloc+3)
      endif
    enddo
  endif
!
!  Generate second derivatives
!
  if (lnoanald2loc.or.lfinitediff) then
    call dynamicn
  else
    call dynamic(0.0_dp,0.0_dp,0.0_dp)
  endif
!
  if (.not.lnoanald2loc.and..not.lfinitediff) then
!
!  Include diagonal blocks, stored in derv3
!
    do i = 1,nphonat
      iloc = atom2local(nphonatptr(i))
      if (iloc.gt.0) then
!
!  Atom is local to this node
!
        indi = 3*(nphonatptr(i)-1)
        indiloc = 3*(iloc-1)
        derv2(indi+1,indiloc+1) = derv2(indi+1,indiloc+1) + derv3(indiloc+1,1)
        derv2(indi+2,indiloc+1) = derv2(indi+2,indiloc+1) + derv3(indiloc+2,1)
        derv2(indi+3,indiloc+1) = derv2(indi+3,indiloc+1) + derv3(indiloc+3,1)
        derv2(indi+1,indiloc+2) = derv2(indi+1,indiloc+2) + derv3(indiloc+1,2)
        derv2(indi+2,indiloc+2) = derv2(indi+2,indiloc+2) + derv3(indiloc+2,2)
        derv2(indi+3,indiloc+2) = derv2(indi+3,indiloc+2) + derv3(indiloc+3,2)
        derv2(indi+1,indiloc+3) = derv2(indi+1,indiloc+3) + derv3(indiloc+1,3)
        derv2(indi+2,indiloc+3) = derv2(indi+2,indiloc+3) + derv3(indiloc+2,3)
        derv2(indi+3,indiloc+3) = derv2(indi+3,indiloc+3) + derv3(indiloc+3,3)
      endif
    enddo
  endif
!*****************************************************************
!  Compress second derivatives according to partial occupancies  *
!*****************************************************************
  if (lpocc) then
! DEBUG - not yet implemented
    ncsfoc = ncfoc + nsfoc
    call compressd2(derv2,maxd2,ncfoc,nsfoc,nbfoc,numat,iocptr,ibocptr)
  elseif (numat.ne.nphonat) then
!*********************************************************************
!  Compress full second derivatives down to region 1 only if needed  *
!*********************************************************************
! DEBUG - not yet implemented - unless phonon atoms are all first
    if (nphonatptr(nphonat).ne.nphonat) then
      call outerror('compression of region 1 not implemented with distributed 2nd derivs',0_i4)
      call stopnow('secondorderfcd')
!
      do i = 1,nphonat
        indi  = 3*(i-1)
        indii = 3*(nphonatptr(i)-1)
        do j = 1,nphonat
          indj  = 3*(j-1)
          indjj = 3*(nphonatptr(j)-1)
          derv2(indj+1,indi+1) = derv2(indjj+1,indii+1)
          derv2(indj+2,indi+1) = derv2(indjj+2,indii+1)
          derv2(indj+3,indi+1) = derv2(indjj+3,indii+1)
          derv2(indj+1,indi+2) = derv2(indjj+1,indii+2)
          derv2(indj+2,indi+2) = derv2(indjj+2,indii+2)
          derv2(indj+3,indi+2) = derv2(indjj+3,indii+2)
          derv2(indj+1,indi+3) = derv2(indjj+1,indii+3)
          derv2(indj+2,indi+3) = derv2(indjj+2,indii+3)
          derv2(indj+3,indi+3) = derv2(indjj+3,indii+3)
        enddo
      enddo
    endif
  endif
!
!  Output the uncompressed second derivatives for debugging
!
  if (index(keyword,'force').ne.0.and.index(keyword,'debu').ne.0) then
    if (ioproc) then
      write(ioout,'(/,''  Uncompressed Force Constant matrix :'',/)')
    endif
    call mpbarrier
    do i = 1,numat
      iloc = atom2local(i)
      if (iloc.gt.0) then
        ix = 3*(iloc-1) + 1
        iy = ix + 1
        iz = ix + 2
        write(ioout,'(12f11.6)')(derv2(j,ix),j=1,mcv+msv)
        write(ioout,'(12f11.6)')(derv2(j,iy),j=1,mcv+msv)
        write(ioout,'(12f11.6)')(derv2(j,iz),j=1,mcv+msv)
      endif
      call mpbarrier
    enddo
  endif
!**********************************
!  Eliminate shell contributions  *
!**********************************
  if (lshelleliminate.and.msv.gt.0) then
!**************************
!  Real Matrix Inversion  *
!**************************
    ifail = 0                   
!
!  Call library to invert matrix 
!
    call matrix_inversion_shells(msv,mcv+1_i4,maxd2,derv2,nshell,nshellonnode,ifail)
!
!  Check return flag
!
    if (ifail.ne.0) then
      call outerror('inversion of shell 2nd derivatives failed',0_i4)
      call stopnow('secondorderfcd')
    endif
!***********************************************
!  Corrected second derivatives = R - T*S-1*T  *
!***********************************************
!
!  Set up Blacs descriptors for matrices
!
    nb = nblocksize
    ifails = 0
    ncs = mcv + msv
    ld = maxd2
    call descinit( idesc, ncs, ncs, 3*nb, 3*nb, 0, 0, iBlacsContext, ld, ifails )
!
    if (ifails.ne.0) then
      call outerror('initialisation in descinit failed',0_i4)
      call stopnow('secondorderfcd')
    endif
!
!  First pass : S-1*T stored in diagonally opposite copy of T
!
    call pdgemm('N','T',msv,mcv,msv,1.0d0,derv2,mcv+1,mcv+1,idesc,derv2,1,mcv+1,idesc,0.0d0,derv2,mcv+1,1,idesc)
!
!  Second pass : T*(S-1*T)
!
    call pdgemm('N','N',mcv,mcv,msv,-1.0d0,derv2,1,mcv+1,idesc,derv2,mcv+1,1,idesc,1.0d0,derv2,1,1,idesc)
  endif
!****************************
!  End of shell correction  *
!****************************
!
!  If debugging print out dynamical matrix
!
  if (index(keyword,'force').ne.0) then
    if (ioproc) then
      write(ioout,'(/,''  Real Force Constant Matrix :'',/)')
    endif
    call mpbarrier
    do i = 1,mcv
      iloc = mcvrptr(i)
      if (iloc.gt.0) then
        write(ioout,'(12f11.6)')(derv2(j,iloc),j=1,mcv)
      endif
      call mpbarrier
    enddo
  endif
!
  deallocate(mcvrptr,stat=status)
  if (status/=0) call deallocate_error('secondorderfcd','mcvrptr')
#else
  call outerror('secondorderfcd called when not compiled with MPI',0_i4)
  call stopnow('secondorderfcd')
#endif
!
  return
  end
