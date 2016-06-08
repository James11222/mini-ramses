!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine kick_part(xpart, vpart, levelp, nparts, grid_level, nbits_patch, previous_timestep)
  use amr_parameters, only: ndim, dp, int_pre, MASK_VALUE
  use amr_commons,    only: ncpu, ind_table2, boxlen, dtnew, dtold, operation_kick, domain_decompos_amr, grid_dict
  use pm_utils,       only: patched_particle_loop
  implicit none
  integer, intent(in) :: grid_level, nparts
  integer, value, intent(in) :: nbits_patch
  integer, dimension(:), intent(inout) :: levelp
  real(dp), dimension(:, :), intent(inout) :: xpart, vpart
  logical, intent(in), value :: previous_timestep
  

  integer :: idim, patch_size, patch_size_coarse
  real(dp), allocatable, dimension(:,:,:,:), target :: f_tmp_fine, f_tmp_coarse
  real(dp), dimension(:,:,:,:), pointer :: f_tmp
  real(dp) :: dx

  if (nparts==0)return  

  call open_cache(operation_kick,domain_decompos_amr)

  patch_size = 2 ** nbits_patch
  patch_size_coarse  = max(patch_size / 2, 2)
  dx = boxlen * 0.5d0 ** grid_level
  
  ! Allocate two cell-thick boundaries to make the depostion onto the AMR grid simpler.
  allocate(f_tmp_fine(-2: patch_size + 1, -2: patch_size + 1, -2: patch_size + 1, 1:ndim))
  allocate(f_tmp_coarse(-2: patch_size_coarse + 1, -2: patch_size_coarse + 1, -2: patch_size_coarse + 1, 1:ndim))

  call patched_particle_loop(xpart, nparts, grid_level, 3, kick_part_callback)

  call close_cache(grid_dict)

  deallocate(f_tmp_fine, f_tmp_coarse)

contains

  subroutine kick_part_callback(oft, np, grid_offset)
    use amr_parameters,         only: nvector
    use particle_interpolation, only: cic_nvector
    use pm_utils,               only: patch_to_AMR
    implicit none
    integer, intent(in), value :: oft, np
    integer(int_pre), dimension(1: ndim) :: grid_offset

    integer(int_pre), dimension(1:nvector, 1:ndim, 0:7) :: ix
    real(dp),         dimension(1:nvector, 0:7)         :: vol
    real(dp),         dimension(1:nvector, 1:ndim)      :: ap
    real(dp),         dimension(1:nvector)              :: dteff
    logical,          dimension(1:nvector)              :: repeat_coarser
    integer(int_pre), dimension(1: ndim)                :: grid_offset_coarse
      
    integer :: idim, ip, ipart, icell, sweep_offset, sweep_nparts
    logical :: all_ok
  
    grid_offset_coarse = grid_offset / 2
    f_tmp => f_tmp_fine
    call patch_to_AMR(grid_offset, patch_size, grid_level, load_f_tmp_callback)
    f_tmp => f_tmp_coarse
    call patch_to_AMR(grid_offset_coarse, patch_size_coarse, grid_level - 1, load_f_tmp_callback)

    ! Loop particles in nvector sweeps
    do sweep_offset = 0, np - 1, nvector
       sweep_nparts = min(np - sweep_offset, nvector)

       ! Get cloud corner integer coordinates and cloud fractions
       call cic_nvector(xpart(oft + sweep_offset + 1: oft + sweep_offset + sweep_nparts, 1:ndim), ix, vol, sweep_nparts, dx)

       do icell = 0, 7
          do idim = 1, ndim
             do ip = 1, sweep_nparts
                ix(ip, idim, icell) = ix(ip, idim, icell) - grid_offset(idim)
             end do
          end do
       end do

       repeat_coarser = .false.
       all_ok = .true.
       ap = 0.d0
       do icell = 0, 7
          do ip = 1, sweep_nparts
             if (f_tmp_fine(ix(ip, 1, icell), ix(ip, 2, icell), ix(ip, 3, icell), 1) == MASK_VALUE) then
                repeat_coarser(ip) = .true.
                all_ok = .false.
             else
                ap(ip, 1:ndim) =  ap(ip, 1:ndim) + vol(ip, icell) * f_tmp_fine(ix(ip, 1, icell), ix(ip, 2, icell), ix(ip, 3, icell), 1:ndim)
             end if
          end do
       end do

       ! For particles which are partially in a coarser level, repeat at coarse level
       do ip = 1, sweep_nparts
          if(repeat_coarser(ip)) then
             ! Maybe write cic_one subroutine...
             call cic_nvector(xpart(oft + sweep_offset + ip: oft + sweep_offset + ip, 1:ndim), ix(1:1, 1:ndim, 0:7), vol(1:1, 0:7), 1, 2 * dx)
             ap(1:3, ip) = 0.d0
             do icell = 0, 7
                ix(1, 1:ndim, icell) = ix(1, 1:ndim, icell) - grid_offset_coarse(1:ndim)
                ap(ip, 1:ndim) =  ap(ip, 1:ndim) + vol(1, icell) * f_tmp_coarse(ix(1, 1, icell), ix(1, 2, icell), ix(1, 3, icell), 1:ndim)
             end do
          end if
       end do

       ! Compute individual time step
       if (previous_timestep)then
          do ip = 1, sweep_nparts
             ipart = oft + sweep_offset + ip 
             if(levelp(ipart) >= grid_level)then
                dteff(ip) = 0.5d0 * dtnew(levelp(ipart))
             else
                dteff(ip) = 0.5d0 * dtold(levelp(ipart))
             endif
          end do
       else
          dteff = 0.5d0 * dtnew(grid_level)
       end if

       ! Finally, apply the kick
       do idim = 1, ndim
          do ip = 1, sweep_nparts
             ipart = oft + sweep_offset + ip
             vpart(ipart, idim) = vpart(ipart, idim) + ap(ip, idim) * dteff(ip)
          end do
       end do
    end do
  end subroutine kick_part_callback

  subroutine load_f_tmp_callback(ix, grid_index)
    use amr_parameters,  only: ndim, int_pre, ngridmax
    use amr_commons,     only: ind_table2, grid
    implicit none
    integer(int_pre), dimension(1:ndim)      :: ix
    integer                                  :: grid_index
    
    if (grid_index > 0) then
       f_tmp(ix(1): ix(1) + 1, ix(2): ix(2) + 1, ix(3): ix(3) + 1, 1:ndim) = RESHAPE(grid(grid_index)%f(1:8, 1:ndim), (/2,2,2,ndim/))
    else
       f_tmp(ix(1): ix(1) + 1, ix(2): ix(2) + 1, ix(3): ix(3) + 1, 1:ndim) = MASK_VALUE
    end if
  end subroutine load_f_tmp_callback  
end subroutine kick_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine drift(ilevel)
  use pm_commons,     only: xp, vp, headp, tailp
  use amr_parameters, only: ndim, boxlen
  use amr_commons,    only: dtnew
  implicit none

  integer, intent(in) :: ilevel
  integer :: idim, ipart
  
  ! Update position
  do idim = 1, ndim
     do ipart = headp(ilevel), tailp(ilevel)
        xp(ipart, idim) = xp(ipart, idim) &
             + vp(ipart, idim) * dtnew(ilevel)
     end do
  end do

  ! Fix periodic boundaries
  do idim = 1, ndim
     do ipart = headp(ilevel), tailp(ilevel)
        if(xp(ipart,idim) > boxlen)then
           xp(ipart,idim) = xp(ipart,idim) - boxlen
        end if
        if(xp(ipart, idim) < 0.d0)then
           xp(ipart, idim) = xp(ipart,idim) + boxlen
        end if   
     end do
  end do
  
end subroutine drift
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine update_levelp(ilevel)
  use pm_commons, only: levelp, headp, tailp
  implicit none

  integer, intent(in) :: ilevel
  levelp(headp(ilevel): tailp(ilevel)) = ilevel  

end subroutine update_levelp
