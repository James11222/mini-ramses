module particle_interpolation
    
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

     ! Subroutine to compute the fractional CIC cloud volumes and the
     ! integer coordinates of the respective cells for nvector  particles.
    
     real(dp) :: one_over_dx
     integer :: icloud, idim, ip
     integer, dimension(1:ndim) :: ind
     ! Auxiliary nvector workspace arrays
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

   
   subroutine ngp_nvector(xp, ix, np, dx)
     use amr_parameters, only: dp, ndim, int_pre, nvector
     use amr_commons,    only: ind_table2
     implicit none
     
     real(dp),         dimension(:, :), intent(in)    :: xp
     integer(int_pre), dimension(:, :), intent(inout) :: ix
     integer,                           intent(in)    :: np
     real(dp),                          intent(in)    :: dx

     ! Subroutine to compute the NGP cell integer coordinates
     ! for nvector particles

     real(dp) :: one_over_dx
     integer :: idim, ip
     
     one_over_dx = 1.0D0 / dx
     
     ! Scale to grid-spacing coordinates
     do idim = 1, ndim
        do ip = 1, np
           ix(ip, idim) = floor(xp(ip, idim) * one_over_dx, kind=int_pre)
        end do
     end do
   end subroutine ngp_nvector
   
 end module particle_interpolation
