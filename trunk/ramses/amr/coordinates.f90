module coordinates
  
  use amr_parameters, only: ndim, nvector, dp, int_pre, ngridmax
  use amr_commons, only: ncoarse, twotondim, ind_table, xg
  implicit none

contains
  
  function grid_to_integer_nvector(xgrid, ilevel, n)
    real(dp),        intent(in), dimension(1:nvector, 1:ndim)  :: xgrid
    integer, intent(in) :: ilevel, n
    integer(int_pre), dimension(1:nvector, 1:ndim) :: grid_to_integer_nvector

    integer  :: idim, i
    real(dp) :: fact

    fact = 2.0_dp ** ilevel
    
    do idim = 1, ndim
       do i = 1, n
          grid_to_integer_nvector(i, idim) = floor(fact * xgrid(i, idim), kind=dp)
       end do
    end do

  end function grid_to_integer_nvector

  function grid_to_integer(xgrid, ilevel)
    real(dp),        intent(in), dimension(1:ndim)  :: xgrid
    integer, intent(in) :: ilevel
    integer(int_pre), dimension(1:ndim) :: grid_to_integer
    
    integer  :: idim       
    do idim = 1, ndim
          grid_to_integer(idim) = floor(2.0_dp ** ilevel * xgrid(idim), kind=dp)
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

    do idim = 1, ndim
       do i = 1, n
          ix(i, idim) = floor(xg(grid_index(i), idim) * fact - 0.5_dp, int_pre) + ind_table(ind(i), idim)
       end do
    end do
    
  end subroutine get_cell_cartesian_key
  

  
end module coordinates
