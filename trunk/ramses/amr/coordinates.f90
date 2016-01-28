module coordinates
  use amr_parameters, only: ndim, nvector, dp, int_pre, ngridmax
  use amr_commons, only: ncoarse, twotondim, ind_table, xg
  implicit none

  real(dp), dimension(1:3) :: skip_loc
  real(dp) :: scale
  
contains

  subroutine init_coords
    use amr_parameters, only: ndim, boxlen
    use amr_commons,    only: icoarse_max, icoarse_min, jcoarse_min, kcoarse_min 
    implicit none
    skip_loc=(/0.0d0,0.0d0,0.0d0/)
    if(ndim>0)skip_loc(1)=dble(icoarse_min)
    if(ndim>1)skip_loc(2)=dble(jcoarse_min)
    if(ndim>2)skip_loc(3)=dble(kcoarse_min)
    scale  = boxlen / dble(icoarse_max - icoarse_min + 1)
    
  end subroutine init_coords
  
  function grid_to_integer_nvector(xgrid, ilevel, n)
    real(dp), intent(in), dimension(:,:)  :: xgrid
    integer, intent(in) :: ilevel, n
    integer(int_pre), dimension(1:nvector, 1:ndim) :: grid_to_integer_nvector

    integer  :: idim, i
    real(dp) :: fact

    fact = 2.0_dp ** ilevel
    
    do idim = 1, ndim
       do i = 1, n
          grid_to_integer_nvector(i, idim) = &
               floor(fact * (xgrid(i, idim) - skip_loc(idim)), kind=dp)
       end do
    end do

  end function grid_to_integer_nvector

  function grid_to_integer(xgrid, ilevel)
    real(dp),        intent(in), dimension(:)  :: xgrid
    integer, intent(in) :: ilevel
    integer(int_pre), dimension(1:ndim) :: grid_to_integer
    
    integer  :: idim       
    do idim = 1, ndim
       grid_to_integer(idim) = &
            floor(2.0_dp ** ilevel * (xgrid(idim) - skip_loc(idim)), kind=dp)
    end do
    
  end function grid_to_integer

  subroutine get_cell_cartesian_key(cell_index, ix, ilevel, n)
    implicit none
    integer, dimension(1:nvector), intent(in) :: cell_index
    integer(int_pre), dimension(1:nvector, 1:ndim), intent(inout) :: ix
    integer, intent(in) :: ilevel, n

    integer, dimension(1:nvector) :: grid_index, ind
    integer :: i, idim
    real(kind=8) :: fact
    
    fact = real(2_int_pre**ilevel, kind = 8)
    
    ! Get cell position in the grid
    do i = 1, n
       ind(i) = (cell_index(i) - ncoarse - 1) / ngridmax
    end do
    ! Get grid
    do i = 1, n
       grid_index(i) = cell_index(i) - ncoarse - ind(i) * ngridmax
    end do
    ! Combine to optain cartesian key
    do idim = 1, ndim
       do i = 1, n
          ix(i, idim) = floor((xg(grid_index(i), idim) - skip_loc(idim)) &
               * fact - 0.5_dp, int_pre) + ind_table(ind(i), idim)
       end do
    end do
    
  end subroutine get_cell_cartesian_key
  
!################################################################
!################################################################
!################################################################
!################################################################
subroutine get_cell_index_from_hilbertkey(cell_index, cell_levl, hkey, np, ilevel)
  use amr_parameters, only: int_pre, nvector, nhilbert, ndim
  use hilbert,     only: hilbert_nd_reverse
  implicit none
  integer, intent(in) :: np, ilevel
  integer(kind=8),dimension(:,:), intent(in) :: hkey
  integer, dimension(:), intent(inout) :: cell_levl, cell_index

  integer(int_pre),dimension(1:nvector, 1:ndim) :: ix
  
  call hilbert_nd_reverse(ix, hkey, ilevel, np)

  call get_cell_index_from_cartesian_hash(cell_index, cell_levl, ix, ilevel, np)     

end subroutine get_cell_index_from_hilbertkey
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_cell_index_from_cartesian_hash(cell_index, cell_levl, ix, ilevel, n)
  use amr_commons
  use hash, only: hash_get
  use amr_parameters, only: int_pre
  implicit none

  integer, intent(in) :: n, ilevel
  integer(int_pre), intent(in), dimension(:,:) :: ix
  integer, intent(inout), dimension(:) :: cell_index, cell_levl

  !----------------------------------------------------------------------------
  !----------------------------------------------------------------------------
  integer :: i, idim
  integer(int_pre), dimension(0:ndim) :: hash_key
  integer, dimension(1:nvector) :: ind, igrid
!  logical, dimension(1:nvector) :: same, same2
!  integer, dimension(1:nvector) :: sort_ind
!  integer, save :: tot = 0, skipped = 0
!  integer:: dummy
  
  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     stop
  end if

!  tot = tot + n
  
  ! Construct ind from last digits
  do i = 1, n
     ind(i) = IAND(ix(i, 1), 1_int_pre)
  end do
#if NDIM>1
  do i = 1, n
     ind(i) = ind(i) + IAND(ix(i, 2), 1_int_pre) * 2
  end do
#endif
#if NDIM>2
  do i = 1, n
     ind(i) = ind(i) + IAND(ix(i, 3), 1_int_pre) * 4
  end do
#endif
 

  ! do i = 1, n
  !    sort_ind(i) = i
  ! end do


  ! do i = 2, n
  !    j = i
  !    do while (j > 1 .and. ix(sort_ind(j),1))
  !       if (j
  !    end do
  ! end do
  
  
  ! Check if two cell belong to the same grid -> hash table can be avoided
  ! TODO: what if one of the values is negative? bitwise exclusive or can be negative and thus smaller than 2... 
  ! same(1) = .false.
  ! do i = 2, n   
  !    same(i) = IEOR(ix(i, 1), ix(i - 1, 1)) < 2
  ! end do
  ! do i = 2, n
  !    same(i) = same(i) .and. IEOR(ix(i, 2), ix(i - 1, 2)) < 2
  ! end do
  ! do i = 2, n
  !    same(i) = same(i) .and. IEOR(ix(i, 3), ix(i - 1, 3)) < 2
  ! end do

  ! Probe for grid starting from ilevel, if not present, try coarser
  cell_levl(1:n) = ilevel
  do i = 1, n
     
     ! Check if I can skip accessing the hash table
     ! Access the hash table only if necessary
     ! if (same(i)) then
     !    igrid(i) = igrid(i - 1)
     !    cell_levl(i) = cell_levl(i - 1)
     !    skipped = skipped + 1
     !    cycle
     ! end if
     
     ! Initial hash key
     hash_key(0) = ilevel
     do idim = 1, ndim
        hash_key(idim) = ISHFT(ix(i, idim), -1)
        hash_key(idim) = hash_key(idim) - IBITS(hash_key(idim), 63, 1)
     end do
     igrid(i) = hash_get(grid_dict, hash_key(0:ndim))
     
     ! If nothing found, try coarser
     do while (igrid(i) == 0 .and. cell_levl(i) > 2)
        cell_levl(i) = cell_levl(i) - 1
        hash_key(0) = cell_levl(i)
        do idim = 1, ndim
           hash_key(idim) = ISHFT(hash_key(idim), -1)
           hash_key(idim) = hash_key(idim) - IBITS(hash_key(idim), 63, 1)
        end do
        igrid(i) = hash_get(grid_dict, hash_key(0:ndim))
     end do
  end do

  ! Check if all went well
  ! do i = 1, n
  !    if (igrid(i) == 0) then
  !       write(*,*)"Problem in get_cell_index_from_cartesian_hash"
  !       stop
  !    end if
  ! end do
  
  do i = 1, n
     cell_index(i) = ncoarse + igrid(i) + ind(i) * ngridmax 
  end do
  
!  if (mod(skipped,32768)==0) print*, tot, 1.0 * skipped / tot

end subroutine get_cell_index_from_cartesian_hash

!################################################################
subroutine check_refinements(refined, hkeys, ilevel)
  use amr_parameters,         only: nhilbert, nvector
  use amr_commons,            only: son, bound_key_level
  use particle_communication, only: hilbert_comm, build_communicator, &
                                    part_data_to_domain, domain_data_to_part
  implicit none

  integer,                         intent(in) :: ilevel
  integer(kind=8), dimension(:,:), intent(in) :: hkeys
  integer,        dimension(:), intent(inout) :: refined

  ! This routine checks for the provided hilbert keys if they correspond to 
  ! refined cells at level ilevel.

  integer  :: idata, ioft, nd, ihilbert
  integer, dimension(1:nvector) :: dummy_int
  type(hilbert_comm) :: comm
  
  integer,         allocatable, dimension(:  ) :: remote_refined
  integer(kind=8), allocatable, dimension(:,:) :: remote_keys

  call build_communicator(comm, hkeys, bound_key_level(:, ilevel))
  
#ifndef WITHOUTMPI
  allocate(remote_refined(1:comm%nrecv))
  allocate(remote_keys(1:comm%nrecv, 1:nhilbert))
  do ihilbert = 1, nhilbert
     call part_data_to_domain(comm, hkeys(:, ihilbert), remote_keys(:, ihilbert))
  end do
#endif
  
  ! Probe local cells for refinement (abuse refined to store cell index)
  do ioft = comm%local_oft, comm%local_oft + comm%nlocal - 1, nvector
     nd = min(comm%local_oft + comm%nlocal - ioft, nvector) 
     call get_cell_index_from_hilbertkey(refined(ioft + 1:ioft + nd), &
          dummy_int(1:nd), hkeys(ioft + 1:ioft + nd, 1:nhilbert), nd, ilevel)
  end do

  ! Mark data corresponding to refined cells
  do idata = comm%local_oft + 1, comm%local_oft + comm%nlocal
     if (son(refined(idata)) > 0)then
        refined(idata) = 1
     else
        refined(idata) = 0
     end if
  end do

#ifndef WITHOUTMPI
  ! Probe remote cells for refinement (abuse remote_refined to store cell index) 
  do ioft = 0, comm%nrecv - 1 , nvector
     nd = min(comm%nrecv - ioft, nvector) 
     call get_cell_index_from_hilbertkey(remote_refined(ioft + 1:ioft + nd), &
          dummy_int(1:nd), remote_keys(ioft + 1:ioft + nd, 1:nhilbert), nd, ilevel)
  end do
#endif

  ! Mark bins corresponding to refined cells
  do idata = 1, comm%nrecv
     if (son(remote_refined(idata))>0)then
        remote_refined(idata) = 1
     else
        remote_refined(idata) = 0
     end if
  end do

#ifndef WITHOUTMPI
  ! Send refinement information back
  call domain_data_to_part(comm, remote_refined, refined)
  deallocate(remote_refined)
  deallocate(remote_keys)
#endif
end subroutine check_refinements
!#########################################################################

  
end module coordinates
