recursive subroutine cic(xpart, array_size, cell_index, vol, offset, np, cic_level, level_boundary_case)
   use amr_parameters,  only: static, dp, twotondim, int_pre
   use amr_commons,     only: boxlen, nvector, ndim, ind_table2
   use hilbert,         only: hilbert3d
   implicit none

   integer,  intent(in)                                              :: offset, np, array_size, cic_level, level_boundary_case
   integer(kind=4), intent(inout), dimension(1:nvector, 0:twotondim - 1) :: cell_index
   real(dp),        intent(inout), dimension(1:nvector, 0:twotondim - 1) :: vol
   real(dp), intent(in), dimension(1:array_size, 1:ndim)             :: xpart

   ! Subroutine to do the Cloud-in-Cell interpolation for nvector particle positions at level cic_level.

   ! Input variables
   ! - xpart: particle coordinates
   ! - cic_level: grid level at which interpolation takes place)
   ! - np: number of particles
   ! - level_boundary_case: 1 (default) if CIC is performed only at cic_level or 2 if
   !   interpolation should be repeated at level cic_level - 1 whenever parts of the
   !   CIC cloud fall outside the boundaries of cic_level
                            
   ! Output variables:
   ! - cell_index: indices of the cells touched by the CIC clouds
   ! - vol: Fractional volume of of a particle that falls in a given cell

   integer(int_pre), dimension(1:nvector, 1:ndim) :: ix
   integer(kind=4), dimension(1:nvector)         :: cell_level
   real(dp),   dimension(1:nvector, 0:1, 1:ndim) :: cloud_boundary
   real(dp),        dimension(1:nvector, 1:ndim) :: xpart_grid
   logical,         dimension(1:nvector)         :: repeat_coarser
   integer,         dimension(1:ndim)            :: ind
   real(dp),        dimension(1:ndim)            :: delta
   integer  :: idim, ind_cloud, ip
   real(dp) :: part_to_grid
   integer(int_pre) :: grid_size, grid_size_one

   grid_size = 2_int_pre**cic_level
   grid_size_one = grid_size - 1
   
   if (level_boundary_case==2) repeat_coarser = .false.

   ! Convert particle coordinates (0 to boxlen)
   ! into grid-coordinates (0 to 2.**grid_level)
   part_to_grid = 2.0**cic_level / boxlen
   xpart_grid(1:np, 1:ndim) = xpart(offset + 1 : offset + np, 1:ndim) * part_to_grid

   ! Compute distances of cloud boundary from nearest "integer coordinate"
   do idim = 1, ndim       

      ! upper/right/front boundary
      do ip=1,np
         cloud_boundary(ip,1,idim) = xpart_grid(ip, idim) + 0.5D0
      end do

      ! upper/rigt/front boundary rel to nearest integer
      do ip=1,np
         cloud_boundary(ip,1,idim) = cloud_boundary(ip,1,idim) - floor(cloud_boundary(ip,1,idim), kind=dp)
      end do

      ! lower/left/back boundary rel to nearest integer
      do ip=1,np
         cloud_boundary(ip,0,idim) = 1.0D0 - cloud_boundary(ip,1,idim)
      end do
   end do


   ! Loop cloud/cell intersections
   do ind_cloud = 0, twotondim - 1
      ind(1:ndim) = ind_table2(1:ndim, ind_cloud) 

      ! Compute cloud volume
      do ip=1,np
#if NDIM==1
         vol(ip, ind_cloud) = cloud_boundary(ip,ind(1),1) 
#endif
#if NDIM==2
         vol(ip, ind_cloud) = cloud_boundary(ip,ind(1),1) * &
              cloud_boundary(ip,ind(2),2) 
#endif
#if NDIM==3
         vol(ip, ind_cloud) = cloud_boundary(ip,ind(1),1) * &
              cloud_boundary(ip,ind(2),2) * &
              cloud_boundary(ip,ind(3),3) 
#endif
      end do
      
         
      ! Compute cloud corner offset from cloud center
      delta(1:ndim) = ind(1:ndim) - 0.5D0       
      

      ! TODO: Non-periodic boundaries...
      do idim = 1, ndim
         do ip = 1, np
            ix(ip, idim) = floor(xpart_grid(ip,idim) + delta(idim), kind=int_pre)
         end do
      end do
      do idim = 1, ndim
         do ip = 1, np
            ix(ip, idim) = IAND(ix(ip, idim) + grid_size, grid_size_one)
         end do
      end do
      

      ! Get cell indices and levels where the cloud corners fall into
      call get_cell_index_from_cartesian_hash(cell_index(1:np, ind_cloud), cell_level(1:np), &
            ix, cic_level, np)  
      
      ! Exclude cloud fraction which lies in coarser level
      if (level_boundary_case == 1)then
         do ip = 1, np        
            if(cell_level(ip) < cic_level)then
               vol(ip, ind_cloud) = 0.d0
            end if
         end do
      end if
      if (level_boundary_case == 2)then
         do ip = 1, np        
            if(cell_level(ip) < cic_level)then
               repeat_coarser(ip) = .true.
            end if
         end do
      end if
   end do ! end loop over cloud/cell intersections

   if (level_boundary_case == 2)then
      do ip = 1, np        
         if (repeat_coarser(ip)) then
            call cic(xpart, array_size, cell_index(ip, 0), vol(ip, 0), offset + ip - 1, 1, cic_level - 1, 1)
         end if
      end do
   end if
      
end subroutine cic

