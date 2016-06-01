


!################################################################
subroutine levelsort_particles(ilevel)
  use pm_commons,     only: part_level_offset, part_ind_permutation, part_ind_permutation2
  use sort,           only: apply_particle_permutation, lsd_radix_sort_particles
  use amr_parameters, only: nlevelmax, int_pre, ndim, dp
  use pm_commons,     only: xp, levelp, boxlen
  use pm_utils,     only: patched_particle_loop
  use hilbert,      only: hilbert_for_particle
  implicit none

  integer, intent(in) :: ilevel
  
  ! This routine sorts resorts the ilevel and ilevel + 1 particles

  integer  :: offset, ip, patch_size, nparts, nbits_patch, ilev2
  integer  :: unrefined_pos, refined_pos
  logical  :: unrefined
  integer, allocatable, dimension(:,:,:) :: refmap_tmp
  real(dp) :: dx
  
  
  if (ilevel == nlevelmax) return
  offset = part_level_offset(ilevel)
  nparts = part_level_offset(ilevel + 2) - offset
  if (nparts == 0) return


  nbits_patch = 3
  patch_size = 2 ** nbits_patch
  dx = boxlen * 0.5d0 ** ilevel
  
  ! Allocate two cell-thick boundaries to make the depostion onto the AMR grid
  ! simpler.
  allocate(refmap_tmp(-2: patch_size + 1, -2: patch_size + 1, -2: patch_size + 1))

  call patched_particle_loop(xp(offset + 1: offset + nparts, 1: ndim), nparts, ilevel, nbits_patch, levelsort_particles_callback)

  deallocate(refmap_tmp) 

  refined_pos = offset + nparts
  do ip = offset + 1, offset + nparts
     if (btest(levelp(ip), 31)) refined_pos = refined_pos - 1
  end do
  
  ! Find starting indices for refined particles
  unrefined_pos = offset

  ! Set "level boundary" in particle array and rearrange particles
  part_level_offset(ilevel + 1) = refined_pos

  do ip = offset + 1, offset + nparts
     if (btest(levelp(ip), 31))then
        refined_pos = refined_pos + 1
        part_ind_permutation(refined_pos) = ip
     else
        unrefined_pos = unrefined_pos + 1
        part_ind_permutation(unrefined_pos) = ip
     end if
  end do

  do ip = offset + 1, offset + nparts
     levelp(ip) = ibclr(levelp(ip), 31)
  end do
 
  call apply_particle_permutation(offset, nparts, ilevel)

  nparts = part_level_offset(ilevel + 1) - part_level_offset(ilevel)
  ilev2 = ilevel - 2
  ! Compute hilbert keys (probably move outside of this routine)
  call hilbert_for_particle(offset, nparts, 0, ilevel - 3)
  call lsd_radix_sort_particles(offset, nparts, ilevel - 3, ilevel - 3, .true.)
  call apply_particle_permutation(offset, nparts, ilevel)

  
contains
  subroutine levelsort_particles_callback(oft, np, grid_offset)
    use amr_parameters,         only: nvector
    use particle_interpolation, only: ngp_nvector
    use pm_utils,               only: patch_to_AMR
    implicit none
    integer, intent(in), value :: oft, np
    integer(int_pre), dimension(1: ndim) :: grid_offset

    integer(int_pre), dimension(1: nvector, 1:ndim) :: ix
    integer :: sweep_offset, sweep_nparts, ip, idim
    
    call patch_to_AMR(grid_offset, patch_size, ilevel, load_refmap_tmp_callback)

    ! Loop particles in nvector sweeps
    do sweep_offset = 0, np - 1, nvector
       sweep_nparts = min(np - sweep_offset, nvector)
       
       ! Get cloud corner integer coordinates and cloud fractions
       call ngp_nvector(xp(offset + oft + sweep_offset + 1: offset + oft + sweep_offset + sweep_nparts, 1:ndim), ix, sweep_nparts, dx)

       do idim = 1, ndim
          do ip = 1, sweep_nparts
             ix(ip, idim) = ix(ip, idim) - grid_offset(idim)
          end do
       end do
       
       do ip = 1, sweep_nparts
          if (.not. refmap_tmp(ix(ip, 1), ix(ip, 2), ix(ip, 3)) == 0)then
             levelp(offset + oft + sweep_offset + ip) = ibset(levelp(offset + oft + sweep_offset + ip), 31)
          end if
       end do
    end do
  end subroutine levelsort_particles_callback

  subroutine load_refmap_tmp_callback(ix, grid_index)
    use amr_parameters,  only: ndim, int_pre, ngridmax
    use amr_commons,     only: ind_table2, ncoarse, son
    implicit none
    integer(int_pre), dimension(1:ndim) :: ix, ixg
    integer                             :: grid_index, cell_index
    
    integer :: icell
    if (grid_index > 0) then
       !          do icell = 0, 7
       !             ixg(1:ndim) = ix(1:ndim) + ind_table2(1:ndim, icell)
       !refmap_tmp(ixg(1), ixg(2), ixg(3)) = son(ncoarse + icell * ngridmax + grid_index)
       !          end do
       refmap_tmp(ix(1)    , ix(2)    , ix(3)    ) = son(ncoarse + grid_index               )
       refmap_tmp(ix(1) + 1, ix(2)    , ix(3)    ) = son(ncoarse + grid_index +     ngridmax)
       refmap_tmp(ix(1)    , ix(2) + 1, ix(3)    ) = son(ncoarse + grid_index + 2 * ngridmax)
       refmap_tmp(ix(1) + 1, ix(2) + 1, ix(3)    ) = son(ncoarse + grid_index + 3 * ngridmax)
       refmap_tmp(ix(1)    , ix(2)    , ix(3) + 1) = son(ncoarse + grid_index + 4 * ngridmax)
       refmap_tmp(ix(1) + 1, ix(2)    , ix(3) + 1) = son(ncoarse + grid_index + 5 * ngridmax)
       refmap_tmp(ix(1)    , ix(2) + 1, ix(3) + 1) = son(ncoarse + grid_index + 6 * ngridmax)
       refmap_tmp(ix(1) + 1, ix(2) + 1, ix(3) + 1) = son(ncoarse + grid_index + 7 * ngridmax)
    else
       do icell = 0, 7
          ixg(1:ndim) = ix(1:ndim) + ind_table2(1:ndim, icell)
          refmap_tmp(ixg(1), ixg(2), ixg(3)) = 0
       end do
    end if
    
  end subroutine load_refmap_tmp_callback  
end subroutine levelsort_particles
!#########################################################################


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


!#########################################################################
subroutine balance_particles(ilevel)
  use pm_commons,     only: npart, kill_one_particle, &
       xp, vp, mp, idp, levelp, part_hkey, part_level_offset, npartmax, &
       current_state, part_ind_permutation, part_ind_permutation2
  use pm_parameters,  only: npart
  use amr_parameters, only: ndim, dp, nhilbert, id_pre
  use amr_commons,    only: myid, nlevelmax
  use particle_communication, only: hilbert_comm, build_dummy_communicator, part_data_to_domain
  use sort,           only: lsd_radix_sort_particles, apply_particle_permutation
  implicit none
  
  integer, intent(in) :: ilevel

  !----------------------------------------------------------------------------

  !----------------------------------------------------------------------------

  type(hilbert_comm) :: comm
  integer :: np, offset, ioft, i1oft, i1oft_new, npart_new, delta_npart, i, idim, ihilbert, irecv

  real(dp),allocatable,dimension(:,:)       ::xp_recv
  real(dp),allocatable,dimension(:,:)       ::vp_recv
  real(dp),allocatable,dimension(:)         ::mp_recv

  integer(kind=8),allocatable,dimension(:,:) :: part_hkey_recv
  integer ,       allocatable,dimension(:)   :: levelp_recv  
  integer(id_pre),allocatable,dimension(:)   :: idp_recv

  integer :: i1, i2, j1, j2

#ifndef WITHOUTMPI
  
  ioft = part_level_offset(ilevel)
  i1oft = part_level_offset(ilevel + 1)
  np = i1oft - ioft

  call build_dummy_communicator(comm, np)

  ! Compute the delta in number of particles 
  delta_npart = comm%nlocal - comm%ndata + comm%nrecv
  npart_new = npart + delta_npart
  i1oft_new = i1oft + delta_npart
  
  ! Check if there is enough space
  if (npart_new > npartmax)then
     write(*,*) 'too many particles'
     call clean_stop
  end if


  allocate(xp_recv(1:comm%nrecv, 1:ndim))
  allocate(vp_recv(1:comm%nrecv, 1:ndim))
  allocate(mp_recv(1:comm%nrecv))
  allocate(part_hkey_recv(1:comm%nrecv, 1:nhilbert))
  allocate(levelp_recv(1:comm%nrecv))
  allocate(idp_recv(1:comm%nrecv))

  do idim = 1, ndim
     call part_data_to_domain(comm, xp(ioft + 1: ioft + np, idim), xp_recv(:, idim))
     call part_data_to_domain(comm, vp(ioft + 1: ioft + np, idim), vp_recv(:, idim))
  end do
  call part_data_to_domain(comm, mp(ioft + 1: ioft + np), mp_recv)
  call part_data_to_domain(comm, idp(ioft + 1: ioft + np), idp_recv)
  call part_data_to_domain(comm, levelp(ioft + 1: ioft + np), levelp_recv)
  do ihilbert = 1, nhilbert
     call part_data_to_domain(comm, part_hkey(ioft + 1: ioft + np, ihilbert), part_hkey_recv(:, ihilbert))
  end do

  
  ! Make space by moving all higher level particles to the right/left
  ! Important: Removing excess space can only be done after communication!
  if (delta_npart > 0)then
     do i = npart, i1oft + 1, -1 
        xp(i + delta_npart, 1:ndim) = xp(i, 1:ndim)
        vp(i + delta_npart, 1:ndim) = vp(i, 1:ndim)
        mp(i + delta_npart) = mp(i)
        idp(i + delta_npart) = idp(i)
        levelp(i + delta_npart) = levelp(i)
        part_hkey(i + delta_npart, 1:nhilbert) = part_hkey(i, 1:nhilbert)
     end do
  end if

  ! Move particles that stay to the left by their local offset
  do i = ioft + 1, ioft + comm%nlocal
     xp(i, 1:ndim) = xp(i + comm%local_oft, 1:ndim)
     vp(i, 1:ndim) = vp(i + comm%local_oft, 1:ndim)
     mp(i) = mp(i + comm%local_oft)
     idp(i) = idp(i + comm%local_oft)
     levelp(i) = levelp(i + comm%local_oft)
     part_hkey(i, 1:nhilbert) = part_hkey(i + comm%local_oft, 1:nhilbert)
  end do
  
  ! Fill in received particles
  do i = 1, comm%nrecv
     xp(ioft + comm%nlocal + i, 1:ndim) = xp_recv(i, 1:ndim)
     vp(ioft + comm%nlocal + i, 1:ndim) = vp_recv(i, 1:ndim)
     mp(ioft + comm%nlocal + i) = mp_recv(i)
     idp(ioft + comm%nlocal + i) = idp_recv(i)
     levelp(ioft + comm%nlocal + i) = levelp_recv(i)
     part_hkey(ioft + comm%nlocal + i, 1:nhilbert) = part_hkey_recv(i, 1:nhilbert)
  end do
  deallocate(xp_recv, vp_recv, mp_recv, part_hkey_recv, levelp_recv, idp_recv)

  ! Move higher level particles to the left
  if (delta_npart < 0)then
     do i = i1oft + 1, npart
        xp(i + delta_npart, 1:ndim) = xp(i, 1:ndim)
        vp(i + delta_npart, 1:ndim) = vp(i, 1:ndim)
        mp(i + delta_npart) = mp(i)
        idp(i + delta_npart) = idp(i)
        levelp(i + delta_npart) = levelp(i)
        part_hkey(i + delta_npart, 1:nhilbert) = part_hkey(i, 1:nhilbert)
     end do
  end if

  
  ! Update total number of particles and level offset for ilevel + 1
  part_level_offset(ilevel + 1) = i1oft_new
  npart = npart_new
  if (ilevel < nlevelmax)part_level_offset(ilevel + 2 : nlevelmax + 1) = npart
  
  ! Resort particles
  call lsd_radix_sort_particles(ioft, np + delta_npart, ilevel, ilevel, .true.)
  call apply_particle_permutation(ioft, np + delta_npart, ilevel)

#else
  return
#endif

end subroutine balance_particles
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
  use hilbert,      only : ge_keys
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

