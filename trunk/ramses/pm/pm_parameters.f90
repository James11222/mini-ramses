module pm_parameters
  use amr_parameters, ONLY: dp, id_pre
  integer  :: npartmax=0                  ! Maximum number of particles
  integer(id_pre)  :: npart=0                     ! Actual number of particles
  real(dp) :: n_dump_parts_direct = 100000000   ! If a leaf cell contains more than n_dump_parts_direct
                                          ! the histogrammed quantities are used for to compute rho

end module pm_parameters
