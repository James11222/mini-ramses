module cic_stuff
contains

   subroutine cic_nvector(xp, ix, vol, np, dx)
     use amr_parameters, only: dp, ndim, int_pre, nvector
     use amr_commons,    only: ind_table2
     implicit none
     
     real(dp),         dimension(:, :),       intent(in)    :: xp
     integer(int_pre), dimension(1:, 1:, 0:), intent(inout) :: ix
     real(dp),         dimension(1:, 0:),     intent(inout) :: vol
     integer,                                 intent(in)    :: np
     real(dp),                                intent(in)    :: dx

     ! Worker subroutine to compute the fractional CIC cloud volumes and the
     ! integer coordinates of the respective cells for a nvector bunch of particles.
     
     
     real(dp) :: one_over_dx
     integer :: icloud, idim, ip
     integer, dimension(1:ndim) :: ind
     ! Auxiliary workspace arrays
     real(dp), dimension(1:nvector, 0:1, 1:ndim), save :: cloud_boundary
     real(dp), dimension(1:nvector, 1:ndim),      save :: xpart
     
     one_over_dx = 1.0D0 / dx

     ! Scale to grid-spacing coordinates
     do idim = 1, ndim
        do ip = 1, np
           xpart(ip, idim) = xp(ip, idim) * one_over_dx
        end do
     end do

     do idim = 1, ndim
        do ip = 1, np
           cloud_boundary(ip, 1, idim) = xpart(ip, idim) + 0.5D0
        end do
        do ip = 1, np
           ! upper/rigt/front boundary rel to nearest integer
           cloud_boundary(ip, 1, idim) = cloud_boundary(ip, 1, idim) - floor(cloud_boundary(ip, 1, idim), kind=dp)       
        end do
        do ip = 1, np
           ! lower/left/back boundary rel to nearest integer 
           cloud_boundary(ip, 0, idim) = 1.0D0 - cloud_boundary(ip, 1, idim)
        end do
     end do
     
     do icloud = 0, 7
        ind(1:3) = ind_table2(1:3, icloud)

        ! Compute cloud volumes
        do ip = 1, np
           vol(ip, icloud) = cloud_boundary(ip, ind(1), 1) * &
                             cloud_boundary(ip, ind(2), 2) * &
                             cloud_boundary(ip, ind(3), 3)
        end do
        
        ! Compute integer coordinates of each cic cloud corner
        do idim = 1, ndim
           do ip = 1, np
              ix(ip, idim, icloud) = floor(xpart(ip, idim) + ind(idim) - 0.5D0, kind=int_pre)
           end do
        end do
     end do
   end subroutine cic_nvector
   

   subroutine dump_rho_tmp(rho_tmp, grid_offset, patch_size, grid_level)
     use amr_parameters,  only: ndim, dp, int_pre
     use pm_utils,        only: patch_to_AMR
     implicit none
     
     ! Keep ordering of variable declaration like this to compile without errors.
     integer(int_pre), dimension(1:ndim) :: grid_offset
     real(dp), dimension(grid_offset(1) - 2: , &
                         grid_offset(2) - 2: , &
                         grid_offset(3) - 2: , &
                         1 : ), intent(inout) :: rho_tmp
     integer :: grid_level, patch_size
     
     ! Take a temporary regular grid patch of a given size at a given
     ! level with a given grid offset and add the content of it to the
     ! permanent rho/phi in memory.
     
     integer(int_pre), dimension(0:ndim) :: hash_key
     integer(int_pre), dimension(1:ndim) :: ix
     integer :: grid_index
     
     call patch_to_AMR(grid_offset, patch_size, grid_level, dump_rho_tmp_callback)
     
   contains
     
     subroutine dump_rho_tmp_callback(ix, grid_index)
       use amr_parameters,  only: ndim, int_pre, ngridmax
       use amr_commons,     only: ind_table2, ncoarse
       use poisson_commons, only: rho, phi
       implicit none
       integer(int_pre), dimension(1:ndim) :: ix, ixg
       integer                             :: grid_index, cell_index
       
       integer :: icell
       if (grid_index > 0) then
          do icell = 0, 7
             ixg(1:ndim) = ix(1:ndim) + ind_table2(1:ndim, icell)
             cell_index = ncoarse + icell * ngridmax + grid_index
             rho(cell_index) = rho(ncoarse + icell * ngridmax + grid_index) + rho_tmp(ixg(1), ixg(2), ixg(3), 1)
             phi(cell_index) = phi(cell_index) + rho_tmp(ixg(1), ixg(2), ixg(3), 2)
          end do
       end if
     end subroutine dump_rho_tmp_callback
   end subroutine dump_rho_tmp

   subroutine load_f_tmp(f_tmp, grid_offset, patch_size, grid_level)
     use amr_parameters,  only: ndim, dp, int_pre, MASK_VALUE
     use amr_commons,     only: ind_table2, grid_dict, ncoarse, ngridmax
     use poisson_commons, only: f
     use pm_utils,        only: patch_to_AMR
     implicit none
     ! Keep ordering of variable declaration like this to compile without errors.
     integer(int_pre), dimension(1:ndim) :: grid_offset
     real(dp), dimension(grid_offset(1) - 2: , &
                         grid_offset(2) - 2: , &
                         grid_offset(3) - 2: , &
                         1 : ), intent(inout) :: f_tmp
     integer :: grid_level, patch_size

     call patch_to_AMR(grid_offset, patch_size, grid_level, load_f_tmp_callback)
     
   contains
     
     subroutine load_f_tmp_callback(ix, grid_index)
       use amr_parameters,  only: ndim, int_pre, ngridmax
       use amr_commons,     only: ind_table2, ncoarse
       use poisson_commons, only: f
       implicit none
       integer(int_pre), dimension(1:ndim)      :: ix
       integer                                  :: grid_index

       integer(int_pre), dimension(1:ndim, 0:7) :: ixg       
       integer :: icell

       do icell = 0, 7
          ixg(1:3, icell) = ind_table2(1:3, icell) + ix(1:3)
       end do

       if (grid_index > 0) then
          do icell = 0, 7
             f_tmp(ixg(1, icell), ixg(2, icell), ixg(3, icell), 1:ndim) = f(ncoarse + icell * ngridmax + grid_index, 1:ndim)
          end do
       else
          do icell = 0, 7
             f_tmp(ixg(1, icell), ixg(2, icell), ixg(3, icell), 1:ndim) = MASK_VALUE
          end do
       end if
     end subroutine load_f_tmp_callback
   end subroutine load_f_tmp
   
 end module cic_stuff
