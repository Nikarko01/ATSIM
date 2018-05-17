  module kim_functions
!
!  Module that contains the data structures and associated KIM functions within GULP
!  including the neighbour lists.
!
!   7/16 Modified for OpenKIM version 1.7.3
!
!  Julian Gale, Curtin University, July 2016
!
    use datatypes
#if KIM_F03
    use, intrinsic :: iso_c_binding
    integer(c_int), dimension(:),  pointer, public,  save :: kim_neighobject => null()
#elif KIM
    integer(i4),                            public,  save :: kim_neighobject = 0
#endif
    integer(i4),                            private, save :: maxneigh = 12
    integer(i4), dimension(:),     pointer, private, save :: nneigh => null()
    integer(i4), dimension(:,:),   pointer, private, save :: neighno => null()
    real(dp),    dimension(:,:),   pointer, private, save :: rneigh => null()
    real(dp),    dimension(:,:,:), pointer, private, save :: xyzneigh => null()

  contains

  subroutine set_kim_neighbours(nkim)
!
!  Sets up the neighbour list for a KIM model calculation
!
!  On entry : 
!
!  nkim            = integer reference to number of KIM model for this call
!
!  10/12 Created from bondordermd.f90
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
!  Copyright Curtin University 2012
!
!  Julian Gale, CIC, Curtin University, October 2012
!
#if KIM_F03
#include "KIM_API_status.h"
#elif KIM
#include "KIM_API_status.h"
#endif
  use datatypes
  use control,        only : keyword
  use current
  use iochannels
#if KIM_F03
  use kim_models,     only : kim_cutoff
#elif KIM
  use kim_models,     only : kim_cutoff
#endif
  use neighbours
  use parallel
  use spatial,        only : lspatialok
  use spatial,        only : natomcell
  use spatial,        only : natomnodeptr
  use spatial,        only : natompernode
  use spatial,        only : ncellsearch
  use spatial,        only : nspcell
  use spatial,        only : nspcellat
  use spatial,        only : nspcellatptr
  use spatial,        only : nspcellat1ptr
  use spatial,        only : nspcellatptrcell
  use spatial,        only : xinbox
  use spatial,        only : yinbox
  use spatial,        only : zinbox
  use times
!
  implicit none
!
!  Passed variables
!
  integer(i4),   intent(in)                        :: nkim
!
!  Local variables
!
  integer(i4)                                      :: i
  integer(i4)                                      :: ic
  integer(i4)                                      :: ii
  integer(i4)                                      :: imx
  integer(i4)                                      :: imy
  integer(i4)                                      :: imz
  integer(i4)                                      :: ind
  integer(i4)                                      :: ind2
  integer(i4)                                      :: indn
  integer(i4)                                      :: ix
  integer(i4)                                      :: iy
  integer(i4)                                      :: iz
  integer(i4)                                      :: j
  integer(i4)                                      :: jc
  integer(i4)                                      :: jj
  integer(i4)                                      :: k
  integer(i4)                                      :: kc
  integer(i4)                                      :: kmax
  integer(i4)                                      :: l
  integer(i4)                                      :: maxxy
  integer(i4)                                      :: maxx
  integer(i4)                                      :: n1j
  integer(i4)                                      :: nati
  integer(i4)                                      :: natj
  integer(i4)                                      :: ndone
  integer(i4), dimension(:),     allocatable, save :: ndoneptr
  integer(i4)                                      :: nj
  integer(i4)                                      :: nn
  integer(i4)                                      :: nnshell
  integer(i4)                                      :: nsplower(3)
  integer(i4)                                      :: nspupper(3)
  integer(i4)                                      :: ntypi
  integer(i4)                                      :: ntypj
  integer(i4)                                      :: status
  logical,     dimension(:),     allocatable, save :: lalreadydone
  logical,     dimension(:),     allocatable, save :: latomdone
  logical,     dimension(:),     allocatable, save :: latomdone2
  real(dp)                                         :: cut2
  real(dp)                                         :: rij
  real(dp)                                         :: r2
  real(dp)                                         :: xi
  real(dp)                                         :: yi     
  real(dp)                                         :: zi
  real(dp)                                         :: xji
  real(dp)                                         :: yji
  real(dp)                                         :: zji
  real(dp)                                         :: xji0
  real(dp)                                         :: yji0
  real(dp)                                         :: zji0
!
!  Allocate local memory that does not depend on maxneigh             
!
  allocate(lalreadydone(numat),stat=status)
  if (status/=0) call outofmemory('kimmd','lalreadydone')
  allocate(latomdone(numat),stat=status)
  if (status/=0) call outofmemory('kimmd','latomdone')
  allocate(latomdone2(numat),stat=status)
  if (status/=0) call outofmemory('kimmd','latomdone2')
  allocate(ndoneptr(numat),stat=status)
  if (status/=0) call outofmemory('kimmd','ndoneptr')
!
!  Call memory allocation routine to ensure that neighbour list is initialise with the right size
!
  call changemaxneigh
!********************************
!  Set square of cutoff radius  *
!********************************
#if KIM_F03
  cut2 = kim_cutoff(nkim)**2
#elif KIM
  cut2 = kim_cutoff(nkim)**2
#else
  cut2 = 0.0_dp
#endif
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
!
!  Set up logical array of atoms done, so that only those needed are done in parallel
!
  latomdone(1:numat) = .false.
!
!  Compute neighbour list
!
  call get_kim_neighbours(cut2,latomdone)
!**********************************************************************
!  If maxneigh has been exceeded, increase value and return to start  *
!**********************************************************************
  if (nprocs.gt.1) then
!******************************
!  Parallel additional setup  *
!******************************
!
!  Loop over atoms again to do those that are neighbours of the main atoms
!  since these will be needed in the energy evaluation. This process has to
!  be done once (in contrast to Brenner potential) since there are no
!  torsional terms. However, construct is left in here in case we want to 
!  add torsions later!
!
    if (lspatialok) then 
!******************** 
!  Spatial version  *
!********************
      maxxy = nspcell(1)*nspcell(2)
      maxx  = nspcell(1)
!                 
!  Initialise atom done once pointer
!                   
      ndone = 0
      do i = 1,numat
        lalreadydone(i) = .false.
      enddo
!
!  Loop over shells
!
      do nnshell = 1,1
        latomdone2(1:numat) = latomdone(1:numat)
        if (nnshell.eq.1) then
          kmax = natompernode
        else
          kmax = numat
        endif
        do kc = 1,kmax
          if (nnshell.eq.1) then
            k = natomnodeptr(kc)
          else
            k = kc
          endif
          if (latomdone2(k)) then
            do l = 1,nneigh(k)
              i = neighno(l,k)
              if (.not.latomdone(i)) then
                nati = nat(i)
                ntypi = nftype(i)
                nneigh(i) = 0
!
!  Find cell containing central image of i 
!
                ind = natomcell(i)
                ind2 = ind - 1
                iz = ind2/maxxy
                ind2 = ind2 - maxxy*iz
                iy = ind2/maxx
                ix = ind2 - maxx*iy + 1
                iy = iy + 1 
                iz = iz + 1
!
                xi = xinbox(i)
                yi = yinbox(i)
                zi = zinbox(i)
!
!  Set cell search bounds
!
                nspupper(1) = min(ix+ncellsearch(1),nspcell(1))
                nspupper(2) = min(iy+ncellsearch(2),nspcell(2))
                nspupper(3) = min(iz+ncellsearch(3),nspcell(3))
                nsplower(1) = max(ix-ncellsearch(1),1)
                nsplower(2) = max(iy-ncellsearch(2),1)
                nsplower(3) = max(iz-ncellsearch(3),1)
!
!  Loop over neighbouring cells
!
                do imz = nsplower(3),nspupper(3)
                  do imy = nsplower(2),nspupper(2)
                    do imx = nsplower(1),nspupper(1)
                      indn = (imz-1)*maxxy + (imy-1)*maxx + imx
!
!  Loop over atoms within neighbouring cells
!
                      nj = nspcellat(indn)
                      n1j = nspcellat1ptr(indn)
                      do jj = 1,nj
                        j = nspcellatptr(n1j+jj)
!
!  Exclude self term
!               
                        if (.not.lalreadydone(j)) then
                          if (i.ne.j.or.ind.ne.indn) then
                            if (latomdone(j)) then
!
!  Atom has already been done and therefore information can be copied
!
                              do ii = 1,nneigh(j)
                                if (neighno(ii,j).eq.i) then
!
!  Check size of neighbour list arrays
!
                                  if (nneigh(i).ge.maxneigh) then
                                    maxneigh = nneigh(i) + 10
                                    call changemaxneigh
                                  endif
!
!  Add neighbour to list
!
                                  nneigh(i) = nneigh(i) + 1
                                  neighno(nneigh(i),i) = j
                                  rneigh(nneigh(i),i) = rneigh(ii,j)
                                  xyzneigh(1,nneigh(i),i) = - xyzneigh(1,ii,j)
                                  xyzneigh(2,nneigh(i),i) = - xyzneigh(2,ii,j)
                                  xyzneigh(3,nneigh(i),i) = - xyzneigh(3,ii,j)
                                endif
                              enddo
!
!  Set pointer to avoid repetition of this atom in the copy phase
!
                              ndone = ndone + 1
                              ndoneptr(ndone) = j
                              lalreadydone(j) = .true.
                            else
                              natj = nat(j)
                              ntypj = nftype(j)
!
!  Set centre cell coordinate differences
!
                              jc = nspcellatptrcell(n1j+jj)
                              xji = xvec2cell(jc) + xinbox(j) - xi
                              yji = yvec2cell(jc) + yinbox(j) - yi
                              zji = zvec2cell(jc) + zinbox(j) - zi
                              r2 = xji*xji + yji*yji + zji*zji
!
!  Check distance squared against cutoff
!
                              if (r2.lt.cut2) then
!
!  Check size of neighbour list arrays
!
                                if (nneigh(i).ge.maxneigh) then
                                  maxneigh = nneigh(i) + 10
                                  call changemaxneigh
                                endif
!
!  Add neighbour to list
!
                                rij = sqrt(r2)
                                nneigh(i) = nneigh(i) + 1
                                neighno(nneigh(i),i) = j
                                rneigh(nneigh(i),i) = rij
                                xyzneigh(1,nneigh(i),i) = xji
                                xyzneigh(2,nneigh(i),i) = yji
                                xyzneigh(3,nneigh(i),i) = zji
!
!  Set pointer to avoid repetition of this atom in the search phase
!
                                ndone = ndone + 1
                                ndoneptr(ndone) = j
                                lalreadydone(j) = .true.
                              endif
                            endif
                          endif
                        endif
                      enddo
                    enddo
                  enddo
                enddo
!
!  Set flag for this atom to indicate that it has been done
!
                latomdone(i) = .true.
!
!  Clear already done pointers for this atom
!
                do j = 1,ndone
                  lalreadydone(ndoneptr(j)) = .false.
                enddo
                ndone = 0
              endif
            enddo
          endif
        enddo
      enddo
    else
!************************
!  Non-spatial version  *
!************************
      do nnshell = 1,1
        latomdone2(1:numat) = latomdone(1:numat)
        if (nnshell.eq.1) then
          kmax = natompernode   
        else
          kmax = numat
        endif
        do kc = 1,kmax
          if (nnshell.eq.1) then
            k = natomnodeptr(kc)
          else
            k = kc   
          endif
          if (latomdone2(k)) then
            do l = 1,nneigh(k)
              i = neighno(l,k)
!
              if (.not.latomdone(i)) then
                nneigh(i) = 0
                nati = nat(i)
                ntypi = nftype(i)
!
!  Loop over atoms
!  
                do j = 1,numat
                  if (latomdone(j)) then
!  
!  Atom has already been done and therefore information can be copied      
!  
                    do ii = 1,nneigh(j)
                      if (neighno(ii,j).eq.i) then
!
!  Check size of neighbour list arrays
!
                        if (nneigh(i).ge.maxneigh) then
                          maxneigh = nneigh(i) + 10
                          call changemaxneigh
                        endif
!
!  Add neighbour to list
!
                        nneigh(i) = nneigh(i) + 1
                        neighno(nneigh(i),i) = j
                        rneigh(nneigh(i),i) = rneigh(ii,j)
                        xyzneigh(1,nneigh(i),i) = - xyzneigh(1,ii,j)
                        xyzneigh(2,nneigh(i),i) = - xyzneigh(2,ii,j)
                        xyzneigh(3,nneigh(i),i) = - xyzneigh(3,ii,j)
                      endif     
                    enddo       
                  else
                    natj = nat(j)  
                    ntypj = nftype(j)
!
!  Set centre cell coordinate differences
!
                    xji0 = xclat(j) - xclat(i)
                    yji0 = yclat(j) - yclat(i)
                    zji0 = zclat(j) - zclat(i)
!
!  Loop over unit cells
!
                    do ii = 1,iimax2
!
!  Exclude self term
!
                      if (i.ne.j.or.ii.ne.iimid2) then
                        xji = xji0 + xvec2cell(ii)
                        yji = yji0 + yvec2cell(ii)
                        zji = zji0 + zvec2cell(ii)
                        r2 = xji*xji + yji*yji + zji*zji
!
!  Check distance squared against cutoff
!
                        if (r2.lt.cut2) then
!
!  Check size of neighbour list arrays
!
                          if (nneigh(i).ge.maxneigh) then
                            maxneigh = nneigh(i) + 10
                            call changemaxneigh
                          endif
!
!  Add neighbour to list
!
                          rij = sqrt(r2)
                          nneigh(i) = nneigh(i) + 1
                          neighno(nneigh(i),i) = j
                          rneigh(nneigh(i),i) = rij
                          xyzneigh(1,nneigh(i),i) = xji
                          xyzneigh(2,nneigh(i),i) = yji
                          xyzneigh(3,nneigh(i),i) = zji
                        endif
                      endif
                    enddo
                  endif
                enddo
!
!  Set flag for this atom to indicate that it has been done
!
                latomdone(i) = .true.
              endif
            enddo
          endif
        enddo
      enddo
    endif
  endif
!*********************************************
!  Print debugging information if requested  *
!*********************************************
  if (ioproc.and.index(keyword,'debu').ne.0) then
    write(ioout,'(/,''  Number of neighbours for atoms :'',/)')
    write(ioout,'(''  i'',2x,''Nat'',3x,''N'')')
    do ic = 1,natompernode
      i = natomnodeptr(ic)
      write(ioout,'(i4,1x,i3,1x,i6)') i,nat(i),nneigh(i)
    enddo
    write(ioout,'(/,''  Neighbours of atoms :'',/)')
    do ic = 1,natompernode
      i = natomnodeptr(ic)
      write(ioout,'(i4,8(1x,i4))') i,(neighno(nn,i),nn=1,nneigh(i))
    enddo
  endif
!
!  Free local memory
!
  deallocate(ndoneptr,stat=status)
  if (status/=0) call deallocate_error('kimmd','ndoneptr')
  deallocate(latomdone2,stat=status)
  if (status/=0) call deallocate_error('kimmd','latomdone2')
  deallocate(latomdone,stat=status)
  if (status/=0) call deallocate_error('kimmd','latomdone')
  deallocate(lalreadydone,stat=status)
  if (status/=0) call deallocate_error('kimmd','lalreadydone')
!
  end subroutine set_kim_neighbours
!
!  get_kim_neighbours
!
  subroutine get_kim_neighbours(cut2,latomdone)
!
!  Finds neighbour list for a KIM potential
!
!  On entry : 
!
!  cut2          = maximum potential cutoff from KIM squared
!
!  On exit :
!
!  nneigh        = number of neighbours for each atom
!  neighno       = pointer to atom numbers of neighbours for each atom
!  rneigh        = distances of neighbours for each atom
!  xyzneigh      = x/y/z components of distances of neighbours for each atom
!  latomdone     = logical indicating whether atom was done locally
!
!  10/12 Created from getBOneighbour
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
!  Copyright Curtin University 2012
!
!  Julian Gale, CIC, Curtin University, October 2012
!
  use datatypes
  use current,        only : numat, iimax, iimid
  use current,        only : xclat, yclat, zclat
  use current,        only : xvec1cell, yvec1cell, zvec1cell
  use current,        only : xvec2cell, yvec2cell, zvec2cell
  use spatial,        only : lbuffercell
  use spatial,        only : lspatialok
  use spatial,        only : ncellsearch
  use spatial,        only : ncellpernode
  use spatial,        only : ncellnodeptr
  use spatial,        only : nspcell
  use spatial,        only : nspcellat
  use spatial,        only : nspcellatptr
  use spatial,        only : nspcellat1ptr
  use spatial,        only : nspcellatptrcell
  use spatial,        only : xinbox
  use spatial,        only : yinbox
  use spatial,        only : zinbox
  implicit none
!
!  Passed variables
!
  real(dp),    intent(in)                        :: cut2
  logical,     intent(out)                       :: latomdone(numat)
!
!  Local variables
!
  integer(i4)                                    :: i
  integer(i4)                                    :: ic
  integer(i4)                                    :: ii
  integer(i4)                                    :: imx
  integer(i4)                                    :: imy
  integer(i4)                                    :: imz
  integer(i4)                                    :: ind1
  integer(i4)                                    :: ind2
  integer(i4)                                    :: ix
  integer(i4)                                    :: ixyz
  integer(i4)                                    :: iy
  integer(i4)                                    :: iz
  integer(i4)                                    :: j
  integer(i4)                                    :: jc
  integer(i4)                                    :: jj
  integer(i4)                                    :: maxxy
  integer(i4)                                    :: maxx
  integer(i4)                                    :: n1
  integer(i4)                                    :: n1j
  integer(i4)                                    :: ni
  integer(i4)                                    :: nj
  integer(i4)                                    :: nsplower(3)
  integer(i4)                                    :: nspupper(3)
  real(dp)                                       :: rij
  real(dp)                                       :: r2
  real(dp)                                       :: xi
  real(dp)                                       :: yi     
  real(dp)                                       :: zi
  real(dp)                                       :: xji
  real(dp)                                       :: yji
  real(dp)                                       :: zji
  real(dp)                                       :: xji0
  real(dp)                                       :: yji0
  real(dp)                                       :: zji0
!
!  Initialise nneigh
!
  nneigh(1:numat) = 0
!****************************************
!  Calculate neighbour lists for atoms  *
!****************************************
  if (lspatialok) then
    maxxy = nspcell(1)*nspcell(2)
    maxx  = nspcell(1)
!  
!  Loop over all spatial cells looking for non-buffer cells
!     
    do ixyz = 1,ncellpernode
      if (.not.lbuffercell(ixyz)) then
        ind1 = ncellnodeptr(ixyz)
        ind2 = ind1 - 1
        iz = ind2/maxxy
        ind2 = ind2 - maxxy*iz
        iy = ind2/maxx
        ix = ind2 - maxx*iy + 1
        iy = iy + 1
        iz = iz + 1
! 
!  Set cell search bounds
!  
        nspupper(1) = min(ix+ncellsearch(1),nspcell(1))
        nspupper(2) = min(iy+ncellsearch(2),nspcell(2))
        nspupper(3) = min(iz+ncellsearch(3),nspcell(3))
        nsplower(1) = max(ix-ncellsearch(1),1)
        nsplower(2) = max(iy-ncellsearch(2),1)
        nsplower(3) = max(iz-ncellsearch(3),1)
!     
!  Get number of atoms in this cell
!       
        ni = nspcellat(ind1)
        n1 = nspcellat1ptr(ind1)
!  
!  Loop over atoms in the cell finding neighbours
!     
        do ii = 1,ni
          i = nspcellatptr(n1+ii)
          ic = nspcellatptrcell(n1+ii)
          latomdone(i) = .true.
!       
!  Set coordinates
!
          xi = xinbox(i) + xvec2cell(ic)
          yi = yinbox(i) + yvec2cell(ic)
          zi = zinbox(i) + zvec2cell(ic)
!
!  Loop over neighbouring cells
!
          do imz = nsplower(3),nspupper(3)
            do imy = nsplower(2),nspupper(2)
              do imx = nsplower(1),nspupper(1)
                ind2 = (imz-1)*maxxy + (imy-1)*maxx + imx
!                       
!  Loop over atoms within neighbouring cells  
!                         
                nj = nspcellat(ind2)
                n1j = nspcellat1ptr(ind2)
                do jj = 1,nj
                  j = nspcellatptr(n1j+jj)
!                     
!  Exclude self term    
!                         
                  if (i.ne.j.or.ind1.ne.ind2) then
                    jc = nspcellatptrcell(n1j+jj)
!                             
!  Set centre cell coordinate differences
!  
                    xji = xvec2cell(jc) + xinbox(j) - xi
                    yji = yvec2cell(jc) + yinbox(j) - yi
                    zji = zvec2cell(jc) + zinbox(j) - zi
!  
                    r2 = xji*xji + yji*yji + zji*zji
!
!  Check whether j is within two body cut-off
!
                    if (r2.lt.cut2) then
!
!  Check whether neighbour list arrays are large enough
!
                      if (nneigh(i).ge.maxneigh) then
                        maxneigh = nneigh(i) + 10
                        call changemaxneigh
                      endif
!
!  Save information for new neighbour
!
                      rij = sqrt(r2)
                      nneigh(i) = nneigh(i) + 1
                      neighno(nneigh(i),i) = j
                      rneigh(nneigh(i),i) = rij
                      xyzneigh(1,nneigh(i),i) = xji
                      xyzneigh(2,nneigh(i),i) = yji
                      xyzneigh(3,nneigh(i),i) = zji
                    endif
                  endif
!
                enddo
              enddo
            enddo
          enddo             
!
        enddo                 
!
      endif
    enddo                       
  else
    do i = 1,numat
      latomdone(i) = .true.
!
!  Loop over atoms
!
      do j = 1,numat
!
!  Set centre cell coordinate differences
!
        xji0 = xclat(j) - xclat(i)
        yji0 = yclat(j) - yclat(i)
        zji0 = zclat(j) - zclat(i)
!
!  Loop over unit cells
!
        do ii = 1,iimax
!
!  Exclude self term
!
          if (i.ne.j.or.ii.ne.iimid) then
            xji = xji0 + xvec1cell(ii)
            yji = yji0 + yvec1cell(ii)
            zji = zji0 + zvec1cell(ii)
            r2 = xji*xji + yji*yji + zji*zji
!
!  Check distance squared against cutoff
!
            if (r2.lt.cut2) then
!
!  Check whether neighbour list arrays are large enough
!
              if (nneigh(i).ge.maxneigh) then
                maxneigh = nneigh(i) + 10
                call changemaxneigh
              endif
!
!  Add neighbour to lists
!
              rij = sqrt(r2)
              nneigh(i) = nneigh(i) + 1
              neighno(nneigh(i),i) = j
              rneigh(nneigh(i),i) = rij
              xyzneigh(1,nneigh(i),i) = xji
              xyzneigh(2,nneigh(i),i) = yji
              xyzneigh(3,nneigh(i),i) = zji
            endif
          endif
        enddo
      enddo
    enddo
  endif
!
  end subroutine get_kim_neighbours
!
!  changemaxneigh - reallocates memory for neighbour lists
!
  subroutine changemaxneigh
!
!  Alters the size of the arrays associated with maxneigh and maxat relevant to neighbour list
!
!  10/12 Created
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
!  Copyright Curtin University 2012
!
!  Julian Gale, CIC, Curtin University, October 2012
!
  use datatypes
  use current,        only : maxat
  use reallocate
  implicit none
!
!  Local variables
!
  integer(i4)       :: ierror, i
  integer(i4), save :: oldmaxat = 0
!
  call realloc(nneigh,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxneigh','nneigh')
  call realloc(neighno,maxneigh,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxneigh','neighno')
  call realloc(rneigh,maxneigh,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxneigh','rneigh')
  call realloc(xyzneigh,3_i4,maxneigh,maxat,ierror)
  if (ierror.ne.0) call outofmemory('changemaxneigh','xyzneigh')
!
!  Initialise new parts of data arrays
!
  if (maxat.gt.oldmaxat) then
    do i = oldmaxat+1,maxat
      nneigh(i) = 0
    enddo
  endif
!
!  Save current value of maxat for next call
!
  oldmaxat = maxat
!
  return
  end subroutine changemaxneigh
#if KIM_F03
!*****************************************
!  get_neigh function needed by OpenKIM  *
!*****************************************
  integer(c_int) function get_kim_neigh(pkim,mode,request,atom,numneigh_atom,neigh_atom,vector_atom) bind(c)

  use, intrinsic :: iso_c_binding
  use KIM_API_F03

  use current,     only : numat
  implicit none
!
!  Passed variables
!
  type(c_ptr),              intent(in)  :: pkim             ! C pointer used to reference pkim_model
  integer(c_int),           intent(in)  :: mode
  integer(c_int),           intent(in)  :: request
  integer(c_int),           intent(out) :: atom
  integer(c_int),           intent(out) :: numneigh_atom
  type(c_ptr),              intent(out) :: neigh_atom
  type(c_ptr),              intent(out) :: vector_atom
!  integer(c_int),           intent(out) :: neigh_atom
!  integer(c_int),           intent(out) :: vector_atom
!
!  Local variables
!
  integer(i4),                     save :: atom_next = 0
  integer(i4)                           :: atomToReturn
  integer(c_int)                        :: KIMerror
!
!  If mode = 0 then this is iterator mode
!
  if (mode.eq.0) then
    if (request.eq.0) then ! reset iterator
!
!  If request = 0 then reset atom number to zero
!
      atom_next = 0
      get_kim_neigh = KIM_STATUS_NEIGH_ITER_INIT_OK
      return
    elseif (request.eq.1) then ! increment iterator
!
!  If request = 1 then increment atom number
!
      atom_next = atom_next + 1
      if (atom_next.gt.numat) then
        get_kim_neigh = KIM_STATUS_NEIGH_ITER_PAST_END
        return
      else
        atomToReturn = atom_next
      endif
    else
!
!  Invalid request
!
      KIMerror = kim_api_report_error(__LINE__,"get_kim_neigh","Invalid request",KIM_STATUS_NEIGH_INVALID_REQUEST)
      get_kim_neigh = KIM_STATUS_NEIGH_INVALID_REQUEST
      return
    endif
  elseif (mode.eq.1) then
!
!  Locator mode
!
    if ((request.gt.numat).or.(request.lt.1)) then
      KIMerror = kim_api_report_error(__LINE__,"get_kim_neigh","Invalid atom number",KIM_STATUS_PARTICLE_INVALID_ID)
      get_kim_neigh = KIM_STATUS_PARTICLE_INVALID_ID
      return
    else
      atomToReturn = request
    endif
  else 
!
!  Unknown mode!
!
    KIMerror = kim_api_report_error(__LINE__,"get_kim_neigh","Invalid mode",KIM_STATUS_NEIGH_INVALID_MODE)
    get_kim_neigh = KIM_STATUS_NEIGH_INVALID_MODE
    return
  endif
!
!  Set the return atom
!
  atom = atomToReturn
!
!  Set the number of neighbours for the return atom
!
  numneigh_atom = nneigh(atom)
!
!  Set the neighbours of the atom
!
  neigh_atom = c_loc(neighno(1,atom))
!
!  Set the displacement vectors of the atom
!
  if (associated(xyzneigh)) then
    vector_atom = c_loc(xyzneigh(1,1,atom))
  else
    vector_atom = c_null_ptr
  endif

  get_kim_neigh = KIM_STATUS_OK

  return
  end function get_kim_neigh
#elif KIM
!*****************************************
!  get_neigh function needed by OpenKIM  *
!*****************************************
  integer function get_kim_neigh(pkim,mode,request,atom,numneigh_atom,neigh_atom,vector_atom)

  use KIM_API
  use current,     only : numat
  implicit none
!
!  Passed variables
!
  integer(kind=kim_intptr), intent(in)  :: pkim
  integer,                  intent(in)  :: mode
  integer,                  intent(in)  :: request
  integer,                  intent(out) :: atom
  integer,                  intent(out) :: numneigh_atom
  integer(kind=kim_intptr), intent(out) :: neigh_atom
  integer(kind=kim_intptr), intent(out) :: vector_atom
!
!  Local variables
!
  integer(i4),                     save :: atom_next = 0
  integer(i4)                           :: atomToReturn
  integer(i4)                           :: KIMerror
!
!  If mode = 0 then this is iterator mode
!
  if (mode.eq.0) then
    if (request.eq.0) then ! reset iterator
!
!  If request = 0 then reset atom number to zero
!
      atom_next = 0
      get_kim_neigh = KIM_STATUS_NEIGH_ITER_INIT_OK
      return
    elseif (request.eq.1) then ! increment iterator
!
!  If request = 1 then increment atom number
!
      atom_next = atom_next + 1
      if (atom_next.gt.numat) then
        get_kim_neigh = KIM_STATUS_NEIGH_ITER_PAST_END
        return
      else
        atomToReturn = atom_next
      endif
    else
!
!  Invalid request
!
      KIMerror = kim_api_report_error_f(__LINE__,"get_kim_neigh","Invalid request",KIM_STATUS_NEIGH_INVALID_REQUEST)
      get_kim_neigh = KIM_STATUS_NEIGH_INVALID_REQUEST
      return
    endif
  elseif (mode.eq.1) then
!
!  Locator mode
!
    if ((request.gt.numat).or.(request.lt.1)) then
      KIMerror = kim_api_report_error_f(__LINE__,"get_kim_neigh","Invalid atom number",KIM_STATUS_PARTICLE_INVALID_ID)
      get_kim_neigh = KIM_STATUS_PARTICLE_INVALID_ID
      return
    else
      atomToReturn = request
    endif
  else 
!
!  Unknown mode!
!
    KIMerror = kim_api_report_error_f(__LINE__,"get_kim_neigh","Invalid mode",KIM_STATUS_NEIGH_INVALID_MODE)
    get_kim_neigh = KIM_STATUS_NEIGH_INVALID_MODE
    return
  endif
!
!  Set the return atom
!
  atom = atomToReturn
!
!  Set the number of neighbours for the return atom
!
  numneigh_atom = nneigh(atom)
!
!  Set the neighbours of the atom
!
  neigh_atom = loc(neighno(1,atom))
!
!  Set the displacement vectors of the atom
!
  vector_atom = loc(xyzneigh(1,1,atom))

  get_kim_neigh = KIM_STATUS_OK

  return
  end function get_kim_neigh
#endif
!*********************
!  KIM destroy call  *
!*********************
  subroutine kim_destroy(nkim,icfg)

#if KIM_F03
  use KIM_API_F03
  use kim_models,     only : pkim_model
  use parallel,       only : ioproc
#elif KIM
  use KIM_API
  use kim_models,     only : pkim_model
  use parallel,       only : ioproc
#endif

  implicit none
!
!  Passed variables
!
  integer(i4),   intent(in)                        :: nkim
  integer(i4),   intent(in)                        :: icfg
!
!  Local variables
!
#if KIM_F03
  integer(c_int)                                   :: KIMerror
  integer(c_int)                                   :: KIMerror_d

  KIMerror = kim_api_model_destroy(pkim_model(nkim,icfg))
  if (KIMerror.lt.KIM_STATUS_OK) then
    if (ioproc) then
      KIMerror_d = kim_api_report_error(__LINE__,"kim_destroy","kim_api_model_destroy",KIMerror)
    endif
    call stopnow('kim_destroy')
  endif
#elif KIM
  integer                                          :: KIMerror
  integer                                          :: KIMerror_d

  KIMerror = kim_api_model_destroy_f(pkim_model(nkim,icfg))
  if (KIMerror.lt.KIM_STATUS_OK) then
    if (ioproc) then
      KIMerror_d = kim_api_report_error_f(__LINE__,"kim_destroy","kim_api_model_destroy",KIMerror)
    endif
    call stopnow('kim_destroy')
  endif
#endif

  end subroutine kim_destroy

  end module kim_functions
