subroutine sort_particles(ilevel, use_histograms)
  use pm_commons,     only: npart, part_level_offset, &
                            nbins, bin_keys, part_hkey
  use amr_commons,    only: myid, levelmin, ncpu
  use amr_parameters, only: nhilbert
  use sort,           only: lsd_radix_sort_particles, apply_particle_permutation
  use hilbert,        only: hilbert_for_particle 
  use coordinates,    only: check_refinements
#ifndef WITHOUTMPI
  use particle_communication, only: build_communicator
#endif
  implicit none

  integer, intent(in) :: ilevel
  logical, intent(in) :: use_histograms

  ! This routine moves all particles that sit in refined cells at level ilevel
  ! to level ilevel + 1. It then sorts the remaining ilevel particles 
  ! by hilbert key.
  
  integer :: ndata_remote, ndata_local, local_oft
  integer :: ilev, offset, np, ip
  integer, dimension(1:ncpu, 1:4) :: communicator
  integer, allocatable, dimension(:) :: refined
    
  offset = part_level_offset(ilevel)
  np = npart - part_level_offset(ilevel)

  ! Compute hilbert keys (probably move outside of this routine)
  call hilbert_for_particle(offset, np, 0, ilevel)
  
  ! Compute a permutation that sorts ALL particles starting from offset
  call lsd_radix_sort_particles(offset, np, ilevel, ilevel, .true.)
  
  if (use_histograms)then
     call compute_particle_histogram(offset, np)          
#ifndef WITHOUTMPI
     call build_communicator(communicator, ndata_remote, &
          nbins, ndata_local, local_oft, bin_keys, ilevel)
#endif
     
     allocate(refined(1:nbins))     
     call check_refinements(communicator, ndata_remote, nbins, ndata_local, &
          local_oft, refined, bin_keys(1:nbins, 1:nhilbert), ilevel)
     call levelsort_particles(ilevel, np, nbins, refined, use_histograms)
  else

     ! Need to apply particle permutation here already to have parts sorted
     ! in memory. Using the index insidet build_communicator and
     ! communicate_refinements is possible but will let the code deviate more
     ! from the histogrammed case.
     call apply_particle_permutation(offset, np, ilevel) 
#ifndef WITHOUTMPI
     call build_communicator(communicator, ndata_remote, np, ndata_local, local_oft, &
          part_hkey(offset + 1: offset + np, 1:nhilbert), ilevel)
#endif
     
     allocate(refined(1:np))
     call check_refinements(communicator, ndata_remote, &
          np, ndata_local, local_oft, refined, part_hkey(offset + 1: offset + np, 1:nhilbert), ilevel)

     call levelsort_particles(ilevel, np, np, refined, use_histograms)
  end if

  ! Compute NEW number of particles in ilevel
  np = part_level_offset(ilevel + 1) - part_level_offset(ilevel)  

  ! Re-sort remaining (ilevel particles) (maybe oversort to gain for the cic step!)
  !  call hilbert_for_particle(offset, np, 0, ilevel + 1)
  call lsd_radix_sort_particles(offset, np, ilevel, ilevel, .true.)
  call apply_particle_permutation(offset, np, ilevel)
  !  call hilbert_for_particle(offset, np, 0, ilevel)
  deallocate(refined)
  
end subroutine sort_particles
!################################################################


!################################################################
subroutine levelsort_particles(ilevel, np, ndata, refined, use_histograms)
  use pm_commons,     only: part_level_offset, part_ind_permutation, part_hkey, &
                            bin_keys, part_ind_permutation2
  use sort,           only: gt_keys, apply_particle_permutation
  use amr_parameters, only: nlevelmax, nhilbert
  implicit none

  
  integer, intent(in) :: ilevel, np, ndata
  integer, dimension(1:ndata), intent(in) :: refined
  logical, intent(in) :: use_histograms
  
  ! Sort particles in memory according to their level
  ! by applying a couting sort on the particles.

  integer  :: offset, ibin, ip, ipart
  integer  :: unrefined_pos, refined_pos
  logical  :: unrefined

  if (ilevel == nlevelmax) return
  if (ndata == 0) return
  
  offset = part_level_offset(ilevel)

  ! Find starting indices for refined particles
  unrefined_pos = offset; refined_pos = offset          

  if (use_histograms)then
     ibin = 1; unrefined = (refined(ibin) == 0)     
     do ip = offset + 1, offset + np  
        ipart = part_ind_permutation(ip)
        if (gt_keys(part_hkey(ipart, 1:nhilbert), bin_keys(ibin, 1:nhilbert)))then
           ibin=ibin + 1
        end if
        if (refined(ibin) == 0) then 
           refined_pos = refined_pos + 1
        end if
     end do
  else
     refined_pos = refined_pos + np - sum(refined)
  end if

  ! Set "level boundary" in particle array and rearrange particles
  part_level_offset(ilevel + 1) = refined_pos

  if (use_histograms)then
     ibin = 1; unrefined = (refined(ibin) == 0)
     do ip = offset + 1, offset + np
        ipart = part_ind_permutation(ip)
        if (gt_keys(part_hkey(ipart,1:nhilbert), bin_keys(ibin,1:nhilbert))) then
           ibin = ibin + 1
        end if
        if (refined(ibin) == 1)then
           refined_pos = refined_pos + 1
           part_ind_permutation2(refined_pos) = ipart
        else
           unrefined_pos = unrefined_pos + 1
           part_ind_permutation2(unrefined_pos) = ipart
        end if
     end do
  else
     do ip = offset + 1, offset + np
        ipart = part_ind_permutation(ip)
        if (refined(ipart - offset) == 1)then
           refined_pos = refined_pos + 1
           part_ind_permutation2(refined_pos) = ipart
        else
           unrefined_pos = unrefined_pos + 1
           part_ind_permutation2(unrefined_pos) = ipart
        end if
     end do
  end if
  
  part_ind_permutation(offset + 1:offset + np) = &
       part_ind_permutation2(offset + 1:offset + np)
  call apply_particle_permutation(offset, np, ilevel)

end subroutine levelsort_particles
!#########################################################################


!#########################################################################
subroutine compute_particle_histogram(offset, np)
  use pm_commons, only: part_hkey, bin_keys, bin_count, bin_start_offset, bin_mass, nbins, part_ind_permutation
  use amr_parameters, only: nhilbert
  use sort,        only: gt_keys
  implicit none
  integer, intent(in) :: offset, np

  ! This routine computes a particle histogram for np particles in memory, 
  ! starting from offset + 1 to offset + np.
  ! IMPORTANT: There must be a precomputed array part_ind_permutation which sorts
  ! the particles by hilbert key.

  integer :: ibin, ipart, ip
  integer(kind=8), dimension(1:nhilbert) :: current_bin_key

  ! if there is nothing to do...
  nbins = 0
  if (.not. np > 0)return
  
  ! Count the number of bins
  nbins = 1
  current_bin_key(1:nhilbert) = part_hkey(part_ind_permutation(offset + 1), 1:nhilbert)
  do ip = offset + 2, offset + np
     ipart = part_ind_permutation(ip)
     if (gt_keys(part_hkey(ipart,1:nhilbert), current_bin_key(1:nhilbert)))then
        nbins=nbins+1
        current_bin_key(1:nhilbert) = part_hkey(ipart, 1:nhilbert)
     end if
  end do
  
  ! Allocate histograms
  if(allocated(bin_count))then
     deallocate(bin_keys)
     deallocate(bin_count)
     deallocate(bin_start_offset)
     deallocate(bin_mass)
  end if
  if (.not. allocated(bin_keys))then
     allocate(bin_keys(nbins, 1:nhilbert))
     allocate(bin_count(nbins))
     allocate(bin_start_offset(nbins+1))
     allocate(bin_mass(nbins))
  end if
  bin_mass = 0; bin_count = 0.d0

  ! Label every bin by a key and sum up the particles per bin, store 
  ! the offset of the first particle in each bin in the particle array

  ! First particle in first bin
  ibin=1
  bin_keys(ibin, 1:nhilbert) = part_hkey(part_ind_permutation(offset + 1), 1:nhilbert)  
  bin_count(ibin) = 1.d0
  bin_start_offset(ibin) = offset

  ! All other particles/bins
  current_bin_key(1:nhilbert) = part_hkey(part_ind_permutation(offset + 1), 1:nhilbert)
  do ip = offset + 2, offset + np
     ipart = part_ind_permutation(ip)
     if (gt_keys(part_hkey(ipart,1:nhilbert), current_bin_key(1:nhilbert)))then
        ibin = ibin + 1
        bin_start_offset(ibin) = ip - 1 
        bin_keys(ibin,1:nhilbert) = part_hkey(ipart, 1:nhilbert)
        current_bin_key(1:nhilbert) = part_hkey(ipart, 1:nhilbert)
     end if
     bin_count(ibin) = bin_count(ibin) + 1.d0
  end do
  bin_start_offset(nbins+1) = offset + np

end subroutine compute_particle_histogram
!################################################################




!#########################################################################
subroutine remove_escaped_particles(ilevel)
  use pm_commons,     only: npart, kill_one_particle, &
                            xp, array_pop, part_level_offset
  use pm_parameters,  only: npart
  use amr_commons,    only: nlevelmax
  use amr_parameters, only: boxlen, ndim
  implicit none
  
  integer, intent(in) :: ilevel

  !----------------------------------------------------------------------------
  ! Remove (level >= ilevel) particles that sit in the physical boundary cells.
  ! Periodic boundaries have been taken care of right after the drift, so ALL
  ! particles outside the comp. domain are removed!!!
  !----------------------------------------------------------------------------

  integer :: ipart, ip, idim, np, offset
  logical, allocatable, dimension(:) :: kill

  offset = part_level_offset(ilevel)
  np = npart - part_level_offset(ilevel)

  allocate(kill(offset + 1:offset + np))
  kill = .false.

  ! Flag particles that have left the box
  do idim = 1, ndim
     do ip = offset + 1, offset + np
        if (xp(ip, idim) >= boxlen) kill(ip) = .true.
        if (xp(ip, idim) < 0.d0) kill(ip) = .true.
     end do
  end do

  ! Kill particles, don't forget to update the "kill" flag too...
  ipart = offset + 1 
  do ip = 1, np
     if (kill(ipart))then
        call kill_one_particle(ipart)
        call array_pop(kill(offset + 1:offset + np), ipart - offset)
     else
        ipart = ipart + 1
     end if
  end do

  deallocate(kill)

  ! Reset offsets for all levels above ilevel to
  ! the updated number of particles
  part_level_offset(ilevel + 1: nlevelmax + 1) = npart

end subroutine remove_escaped_particles
!#########################################################################









! ---------------------------------------- UNUSED ROUTINES ---------------------------------------
! These are functioning but currently unused routines which might be valuable for debugging etc...
! -------------------------------------------------------------------------------------------------   
subroutine get_cell_index(cell_index,cell_levl,xpart,ilevel,n)
  use amr_commons
  implicit none

  integer::n,ilevel
  integer,dimension(1:nvector)::cell_index,cell_levl
  real(dp),dimension(1:nvector,1:3)::xpart

  !----------------------------------------------------------------------------
  ! This routine returns the index and level of the cell, (at maximum level
  ! ilevel), in which the input the position specified by xpart lies
  !----------------------------------------------------------------------------

  real(dp)::xx,yy,zz
  integer::i,j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0

  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     call clean_stop
  end if

  ind_cell=0
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,n
     xx = xpart(i,1)/boxlen + (nx-1)/2.0
     yy = xpart(i,2)/boxlen + (ny-1)/2.0
     zz = xpart(i,3)/boxlen + (nz-1)/2.0

     if(xx<0.)xx=xx+dble(nx)
     if(xx>dble(nx))xx=xx-dble(nx)
     if(yy<0.)yy=yy+dble(ny)
     if(yy>dble(ny))yy=yy-dble(ny)
     if(zz<0.)zz=zz+dble(nz)
     if(zz>dble(nz))zz=zz-dble(nz)

     igrid=igrid0
     do j=1,ilevel 
        ii=1; jj=1; kk=1
        if(xx<xg(igrid,1))ii=0
        if(yy<xg(igrid,2))jj=0
        if(zz<xg(igrid,3))kk=0
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index(i)=ind_cell
     cell_levl(i)=j
  end do
end subroutine get_cell_index
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_cell_index_from_cartesian(cell_index,cell_levl,xx,yy,zz,ilevel,n,bit_length)
  use amr_commons
  use amr_parameters, only: int_pre
  implicit none

  integer, intent(in)::n,ilevel,bit_length
  integer,dimension(1:nvector)::cell_index,cell_levl
  integer(int_pre),dimension(1:nvector)::xx,yy,zz
  !----------------------------------------------------------------------------
  ! Routine to obtain the cell index from the cartesina key.
  ! Not used at the moment, keep for debug purpose
  !----------------------------------------------------------------------------
  integer::i,j,ind,iskip,igrid,ind_cell,igrid0
  integer(int_pre)::ii,jj,kk

  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     stop
  end if

  if (bit_length>21)then
     print*, 'bit length too big for now'
  end if
  
  ind_cell=0
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,n
     igrid=igrid0
     do j=1,ilevel 
        ii=ISHFT(xx(i),-bit_length+j)
        jj=ISHFT(yy(i),-bit_length+j)
        kk=ISHFT(zz(i),-bit_length+j)
        ii=mod(ii,2)
        jj=mod(jj,2)
        kk=mod(kk,2)
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)exit
     end do
     cell_index(i)=ind_cell
     cell_levl(i)=j
  end do
end subroutine get_cell_index_from_cartesian
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine check_sorted(offset, np)
  use amr_parameters, only: nhilbert
  use pm_commons,   only : part_hkey, part_ind_permutation
  use sort,         only : ge_keys
  use amr_commons,  only : myid
  implicit none
  integer, intent(in) :: offset, np

  !----------------------------------------------------------------------------
  ! Simple helper routine to check whether the np particles starting form offset
  ! are sorted by hilbert key
  !----------------------------------------------------------------------------
  logical,                         save :: ok
  integer,                         save :: ipart, ip
  integer(kind=8), dimension(1:nhilbert), save :: current_key

  ok =.true.
  
  !  current_key(1:nhilbert) = part_hkey(part_ind_permutation(offset + 1),1:nhilbert)
    current_key(1:nhilbert) = part_hkey(offset + 1,1:nhilbert)
  do ip = offset + 2, offset + np
     !ipart = part_ind_permutation(ip)
     ipart = ip
     if (.not. ge_keys(part_hkey(ipart, 1:nhilbert), current_key(1:nhilbert)))then
        ok=.false.
        print*, "Detected unsorted particles on process", myid
        print*, part_hkey(ipart,1:nhilbert),current_key(1:nhilbert)
     end if
     current_key(1:nhilbert) = part_hkey(ipart,1:nhilbert)
  end do
  if (.not. ok)stop
end subroutine check_sorted

