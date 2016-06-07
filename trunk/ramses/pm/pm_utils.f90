module pm_utils

contains
  subroutine patched_particle_loop(xpart, nparts, grid_level, nbits_patch, particle_callback)
    use amr_parameters, only: dp, int_pre, ndim
    use amr_commons, only: boxlen
    implicit none
    real(dp), dimension(:,:), intent(in) :: xpart
    integer, intent(in) :: nparts, grid_level
    integer, value, intent(in) :: nbits_patch
    
    ! Minimal interface block for the actual callback
    interface
       subroutine callback(oft, np, grid_offset)
         use amr_parameters, only: ndim, int_pre
         implicit none
         integer, intent(in), value :: oft, np
         integer(int_pre), dimension(1:ndim) :: grid_offset
       end subroutine callback
    end interface

    procedure(callback) :: particle_callback

    logical :: evaluate_patch
    integer :: patch_size, idim, ip, ip_offset
    integer(int_pre), dimension(1:ndim) :: grid_offset
    integer(int_pre), dimension(1: ndim) :: ix_next, ix_current
    real(dp) :: part_to_grid


    part_to_grid = 2.d0 ** grid_level / boxlen
    patch_size = 2 ** nbits_patch

    ip_offset = 0 
    do idim = 1, ndim
       ix_next(idim) = xpart(1, idim) * part_to_grid
    end do
    
    do ip = 1, nparts
       ix_current(1: ndim) = ix_next(1: ndim)
       
       ! If current particle is last particle, evaluate!
       if (ip == nparts) then
          evaluate_patch = .true.        
       else
          do idim = 1, ndim
             ix_next(idim) = xpart(ip + 1, idim) * part_to_grid
          end do
          
          ! Check if current and next key are in same patch 
          ! (i.e. each integer key ix differs only by last n bits from the current)         
          evaluate_patch = .false.
          do idim = 1, ndim
             evaluate_patch = evaluate_patch .or. (IEOR(ix_current(idim), ix_next(idim)) > patch_size - 1)
          end do
       end if
       
       if (evaluate_patch)then
          do idim = 1, ndim
             grid_offset(idim) = ISHFT(ISHFT(ix_current(idim), -nbits_patch), nbits_patch)
          end do
          call particle_callback(ip_offset, ip - ip_offset, grid_offset)
          ip_offset = ip
       end if
    end do
  end subroutine patched_particle_loop
  
  
  
  subroutine patch_to_AMR(grid_offset, patch_size, grid_level, interaction_callback)
    use amr_parameters,  only: ndim, dp, int_pre
    use amr_commons,     only: ind_table2, grid_dict
    implicit none
    
    integer(int_pre), dimension(1:ndim), intent(in) :: grid_offset
    integer,value,                       intent(in) :: grid_level, patch_size

    interface
       subroutine callback(cartesian_index, amr_index)
         use amr_parameters, only: ndim, int_pre
         integer(int_pre), dimension(1:ndim) :: cartesian_index
         integer :: amr_index
       end subroutine callback
    end interface
    procedure(callback) :: interaction_callback

    ! Take a temporary regular grid patch of a given size at a given
    ! level with a given grid offset and add the content of it to the
    ! permanent rho/phi in memory.

    integer(int_pre), dimension(0:ndim) :: hash_key
    integer(int_pre), dimension(1:ndim) :: ix
    integer(int_pre) :: key_space_size, bitmask, i, j, k
    integer :: grid_index, idim, get_grid


    key_space_size = 2 ** (grid_level - 1)
    bitmask = key_space_size - 1
    hash_key(0) = int(grid_level, kind=int_pre)

    ! Loop over all octs in the grid patch
    do i = grid_offset(1) / 2 - 1, grid_offset(1) / 2 +  patch_size / 2
       do j = grid_offset(2) / 2  - 1, grid_offset(2) / 2 +  patch_size / 2
          do k = grid_offset(3) / 2 - 1, grid_offset(3) / 2 +  patch_size / 2
             
             ! Construct the hash key
             hash_key(1: ndim) = (/ i, j, k /)
             
             ! Take care of periodic boundaries!
             do idim = 1, ndim
                hash_key(idim) = IAND(key_space_size + hash_key(idim), bitmask)
             end do

             ! Dump the actual mass onto the grid
             grid_index = get_grid(hash_key, grid_dict, .true., .false.)
             ix(1:3) = 2 * (/ i, j, k /) - grid_offset(1:3)
             call interaction_callback(ix, grid_index)
          end do
       end do
    end do
  end subroutine patch_to_AMR
  

end module pm_utils
