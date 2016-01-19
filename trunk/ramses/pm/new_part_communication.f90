module particle_communication
#ifndef WITHOUTMPI
  
  ! Operator overloading for the communication routines
  interface part_data_to_domain
     module procedure part_data_to_domain_i4, part_data_to_domain_i8, part_data_to_domain_dp
  end interface part_data_to_domain
  
  interface domain_data_to_part
     module procedure domain_data_to_part_i4, domain_data_to_part_dp
  end interface domain_data_to_part
  
contains


  subroutine build_communicator(communicator, recv_tot, ndata, local_data, local_data_oft, hkey, ilevel)

    use amr_parameters, only: nlevelmax
    use amr_commons,    only: ncpu, myid, bound_key_level, bound_key
    implicit none

    include 'mpif.h'
    integer, intent(in) ::  ilevel, ndata
    integer, intent(inout) ::  recv_tot, local_data, local_data_oft

    ! Hilbert keys
    integer(kind=8), dimension(:,:), intent(in) :: hkey

    ! Wrap send/receive counters/offsets in one array to make the passing them
    ! around a little more convenient
    integer, dimension(1:ncpu, 1:4), intent(inout) :: communicator

    !----------------------------------------------------------------------------
    ! This routine sets up the communication structure for any kind 3-integer 
    ! quantity (particles, bins, etc.). The input hilbert keys (keys) are assumed
    ! to be sorted.
    !----------------------------------------------------------------------------

    integer :: receive_cpu, idata, info, icpu, idest, isource
    integer :: countrecv, countsend  
    integer, dimension(MPI_STATUS_SIZE,2*ncpu) :: statuses
    integer, dimension(2*ncpu)                 :: reqsend, reqrecv

    if (nlevelmax>20)then 
       print*, 'problem here with precision'
       stop
    end if

    communicator = 0
    recv_tot = 0

    ! count number of bins that need to be sent to every other process
    receive_cpu = 1; local_data = 0
    do idata = 1, ndata
       ! TO DO:
       ! REPLACE this by do while( gt_3_keys(keys..., bound_key_level )):
       ! do while (gt_3keys_individual_input(keys(idata,2),keys(idata,1),keys(idata,0), &
       ! bound_key_level(idata,2),bound_key_level(idata,1),bound_key_level(idata,0)))
       do while (hkey(idata, 1) > bound_key_level(receive_cpu, ilevel) ) 
          receive_cpu = receive_cpu + 1
       end do
       communicator(receive_cpu, 1) = communicator(receive_cpu, 1) + 1
    end do

    ! send 1 integer to every other process (including itself - a bit silly, but who cares...)
    ! results in  receive counter
    countrecv = 0; countsend = 0  
    do isource = 1, ncpu        
       countrecv = countrecv + 1
       call MPI_IRECV(communicator(isource, 3), 1, MPI_INTEGER, isource - 1, &
            1234, MPI_COMM_WORLD, reqrecv(countrecv), info)
    end do
    do idest = 1, ncpu
       countsend = countsend + 1
       call MPI_ISEND(communicator(idest, 1), 1, MPI_INTEGER, idest - 1, 1234, &
            MPI_COMM_WORLD, reqsend(countsend), info)
    end do

    call MPI_WAITALL(ncpu, reqrecv, statuses, info)
    call MPI_WAITALL(ncpu, reqsend, statuses, info)

    ! "Prefix sum" to compute send offsets
    do icpu = 1, ncpu - 1
       communicator(icpu + 1, 2) = communicator(icpu, 1) + communicator(icpu, 2)
    end do

    ! No information to is sent to itself, store the number of local data
    local_data = communicator(myid, 1)
    local_data_oft = communicator(myid, 2)
    communicator(myid, 1) = 0
    communicator(myid, 3) = 0

    do icpu = 1, ncpu - 1
       communicator(icpu + 1, 4) = communicator(icpu, 3) + communicator(icpu, 4)
    end do

    recv_tot=sum(communicator(:, 3))

  end subroutine build_communicator
  !################################################################

  
  !################################################################
  subroutine part_data_to_domain_i4(communicator, send_data, recv_data)
    use amr_commons,   only: ncpu
    implicit none
    include 'mpif.h'
    
    integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
    integer, dimension(:), intent(in) :: send_data
    integer, dimension(:), intent(inout) :: recv_data
    
    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)
    
    
    call MPI_IALLTOALLV(send_data, communicator(:,1), communicator(:,2), MPI_INTEGER, &
         &              recv_data, communicator(:,3), communicator(:,4), MPI_INTEGER, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine part_data_to_domain_i4
  !################################################################
  !################################################################
  subroutine part_data_to_domain_i8(communicator, send_data, recv_data)
    use amr_commons, only: ncpu
    implicit none
    include 'mpif.h'

    integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
    integer(kind = 8), dimension(:), intent(in) :: send_data
    integer(kind = 8), dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, communicator(:,1), communicator(:,2), MPI_INTEGER8, &
         &              recv_data, communicator(:,3), communicator(:,4), MPI_INTEGER8, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine part_data_to_domain_i8
  !################################################################
  !################################################################
  subroutine part_data_to_domain_dp(communicator, send_data, recv_data)
    use amr_commons,   only: ncpu, dp
    implicit none
    include 'mpif.h'

    integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
    real(dp), dimension(:), intent(in) :: send_data
    real(dp), dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, communicator(:,1), communicator(:,2), MPI_DOUBLE_PRECISION, &
         &              recv_data, communicator(:,3), communicator(:,4), MPI_DOUBLE_PRECISION, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine part_data_to_domain_dp
  !################################################################
  !################################################################
  subroutine domain_data_to_part_i4(communicator, send_data, recv_data)
    use amr_commons,   only: ncpu
    implicit none
    include 'mpif.h'

    integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
    integer, dimension(:), intent(in) :: send_data
    integer, dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, communicator(:,3), communicator(:,4), MPI_INTEGER, &
         &              recv_data, communicator(:,1), communicator(:,2), MPI_INTEGER, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine domain_data_to_part_i4
  !################################################################
  !################################################################
  subroutine domain_data_to_part_dp(communicator, send_data, recv_data)
    use amr_commons,   only: ncpu, dp
    implicit none
    include 'mpif.h'

    integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
    real(dp), dimension(:), intent(in) :: send_data
    real(dp), dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, communicator(:,3), communicator(:,4), MPI_DOUBLE_PRECISION, &
         &              recv_data, communicator(:,1), communicator(:,2), MPI_DOUBLE_PRECISION, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine domain_data_to_part_dp
  !################################################################
  !################################################################
  !################################################################
  !################################################################
  subroutine communicate_refinements(communicator, ndata_remote, ndata, ndata_local, ndata_local_oft, &
       refined, hkeys, ilevel)
    use amr_parameters, only: nhilbert, nvector
    use amr_commons,   only: ncpu, myid, bound_key_level, son, nvector, nlevelmax
    implicit none
    include 'mpif.h'
    integer, intent(in) ::  ilevel, ndata
    integer, intent(in) :: ndata_remote, ndata_local, ndata_local_oft
    integer(kind=8), dimension(:,:), intent(inout) :: hkeys
    integer, dimension(1:ncpu, 1:4), intent(in) :: communicator
    integer, dimension(1:ndata), intent(inout) :: refined

    ! This routine sorts particles between ilevel and ilevel + 1
    ! (formerly known as kill_tree_fine).

    ! The routine can only run after:
    ! call compute_particle_histogram(ilevel)
    ! call build_communicator(...)

    ! The "bins" here are the histogram bins and correspond to 
    ! actual grid cells. "local_... variables" denote properties in the
    ! MPI process which hosts the particles, while "remote_... variables" 
    ! are used for properties in the MPI process which hosts the corresponding
    ! leaf-cell.


    integer  :: idata, ioft, nd, i
    integer, dimension(1:nvector) :: dummy_int

    integer,         allocatable, dimension(:  ) :: remote_refined
    integer(kind=8), allocatable, dimension(:,:) :: remote_keys

    allocate(remote_refined(1:ndata_remote))
    allocate(remote_keys(1:ndata_remote, 1:nhilbert))
    do i = 1, nhilbert
       call part_data_to_domain(communicator, hkeys(:,i), remote_keys(:,i))
    end do


    ! Probe local cells for refinement (abuse refined to store cell index)
    do ioft = ndata_local_oft, ndata_local_oft + ndata_local - 1, nvector
       nd = min(ndata_local_oft + ndata_local - ioft, nvector) 
       call get_cell_index_from_hilbertkey(refined(ioft + 1: ioft + nd), &
            dummy_int(1: nd), hkeys(ioft + 1: ioft + nd, 1:nhilbert), nd, ilevel)
    end do

    ! Mark data corresponding to refined cells
    do idata = ndata_local_oft + 1, ndata_local_oft + ndata_local
       if (son(refined(idata)) > 0)then
          refined(idata) = 1
       else
          refined(idata) = 0
       end if
    end do

    ! Probe remote cells for refinement (abuse remote_refined to store cell index) 
    do ioft = 0, ndata_remote - 1 , nvector
       nd = min(ndata_remote - ioft, nvector) 
       call get_cell_index_from_hilbertkey(remote_refined(ioft + 1:ioft + nd), &
            dummy_int(1:nd), remote_keys(ioft + 1:ioft + nd, 1:nhilbert), nd, ilevel)
    end do

    ! Mark bins corresponding to refined cells
    do idata = 1, ndata_remote
       if (son(remote_refined(idata))>0)then
          remote_refined(idata) = 1
       else
          remote_refined(idata) = 0
       end if
    end do

    ! Send refinement information back
    call domain_data_to_part(communicator, remote_refined, refined)

    deallocate(remote_refined)
    deallocate(remote_keys)

  end subroutine communicate_refinements
!################################################################
!################################################################
!################################################################
!################################################################
#endif
end module particle_communication
