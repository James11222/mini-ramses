module particle_communication
  type hilbert_comm
     integer               :: ndata, nlocal, nrecv, local_oft
     integer, dimension(:), allocatable :: send_count, recv_count, send_oft, recv_oft
  end type hilbert_comm
  

#ifndef WITHOUTMPI  
  ! Operator overloading for the communication routines
  interface part_data_to_domain
     module procedure part_data_to_domain_i4, part_data_to_domain_i8, part_data_to_domain_dp
  end interface part_data_to_domain
  
  interface domain_data_to_part
     module procedure domain_data_to_part_i4, domain_data_to_part_dp
  end interface domain_data_to_part
#endif

contains

  subroutine build_communicator(comm, hkey, boundary_keys)
    use amr_commons,    only: ncpu, myid, nlevelmax
    implicit none
#ifndef WITHOUTMPI    
    include 'mpif.h'
#endif
    type(hilbert_comm), intent(inout) :: comm
    integer(kind=8), dimension(0:ncpu), intent(in) :: boundary_keys
    integer(kind=8), dimension(:,:), intent(in) :: hkey
    
    ! This routine sets up the communication structure for any 3-integer hilbert key
    ! quantity (particles, bins, etc.). The input hilbert keys (keys) are assumed
    ! to be sorted. The domain boundaries must be given.

    ! Note: A dummy communicator is constructed for non-MPI case to make the
    ! rest of the code contain less preprocessor directives.
    
    integer :: receive_cpu, idata, info, icpu, idest, isource
    integer :: countrecv, countsend  
    integer, dimension(MPI_STATUS_SIZE,2*ncpu) :: statuses
    integer, dimension(2*ncpu)                 :: reqsend, reqrecv

    allocate(comm%send_count(1:ncpu))
    allocate(comm%recv_count(1:ncpu))
    allocate(comm%send_oft(1:ncpu))
    allocate(comm%recv_oft(1:ncpu))

    if (nlevelmax > 20)then 
       print*, 'problem here with precision'
       stop
    end if

    comm%ndata = size(hkey, 1)
    comm%send_count = 0; comm%recv_count = 0

    ! count number of bins that need to be sent to every other process
    receive_cpu = 1
    do idata = 1, comm%ndata
       ! TO DO:
       ! REPLACE this by do while( gt_3_keys(keys..., bound_key_level )):
       ! do while (gt_3keys_individual_input(keys(idata,2),keys(idata,1),keys(idata,0), &
       ! bound_key_level(idata,2),bound_key_level(idata,1),bound_key_level(idata,0)))
       do while (hkey(idata, 1) > boundary_keys(receive_cpu) ) 
          receive_cpu = receive_cpu + 1
       end do
       comm%send_count(receive_cpu) = comm%send_count(receive_cpu) + 1
    end do

#ifndef WITHOUTMPI
    ! Send 1 integer to every other process (including itself - a bit silly, but who cares...)
    ! Store results in  receive counter
    countrecv = 0; countsend = 0  
    do isource = 1, ncpu        
       countrecv = countrecv + 1
       call MPI_IRECV(comm%recv_count(isource), 1, MPI_INTEGER, isource - 1, &
            1234, MPI_COMM_WORLD, reqrecv(countrecv), info)
    end do
    do idest = 1, ncpu
       countsend = countsend + 1
       call MPI_ISEND(comm%send_count(idest), 1, MPI_INTEGER, idest - 1, 1234, &
            MPI_COMM_WORLD, reqsend(countsend), info)
    end do

    call MPI_WAITALL(ncpu, reqrecv, statuses, info)
    call MPI_WAITALL(ncpu, reqsend, statuses, info)
#endif
    
    ! "Prefix sum" to compute send offsets
    comm%send_oft(1) = 0
    do icpu = 1, ncpu - 1
       comm%send_oft(icpu + 1) = comm%send_oft(icpu) + comm%send_count(icpu)
    end do

    ! No information is sent to itself, store the number of local data
    comm%nlocal = comm%send_count(myid)
    comm%local_oft = comm%send_oft(myid)
    comm%send_count(myid) = 0
    comm%recv_count(myid) = 0

    ! "Prefix sum" to compute recv offsets
    comm%recv_oft(1) = 0
    do icpu = 1, ncpu - 1
       comm%recv_oft(icpu + 1) = comm%recv_oft(icpu) + comm%recv_count(icpu)
    end do

    ! Total number of received data
    comm%nrecv = sum(comm%recv_count(:))

  end subroutine build_communicator
  !################################################################

#ifndef WITHOUTMPI      
  !################################################################
  subroutine part_data_to_domain_i4(comm, send_data, recv_data)
    implicit none
    include 'mpif.h'
    
    type(hilbert_comm), intent(in) :: comm
    integer, dimension(:), intent(in) :: send_data
    integer, dimension(:), intent(inout) :: recv_data
    
    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)
    
    call MPI_IALLTOALLV(send_data, comm%send_count, comm%send_oft, MPI_INTEGER, &
         &              recv_data, comm%recv_count, comm%recv_oft, MPI_INTEGER, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine part_data_to_domain_i4
  !################################################################
  !################################################################
  subroutine part_data_to_domain_i8(comm, send_data, recv_data)
    implicit none
    include 'mpif.h'
    
    type(hilbert_comm), intent(in) :: comm
    integer(kind = 8), dimension(:), intent(in) :: send_data
    integer(kind = 8), dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, comm%send_count, comm%send_oft, MPI_INTEGER8, &
         &              recv_data, comm%recv_count, comm%recv_oft, MPI_INTEGER8, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine part_data_to_domain_i8
  !################################################################
  !################################################################
  subroutine part_data_to_domain_dp(comm, send_data, recv_data)
    use amr_commons,   only: dp
    implicit none
    include 'mpif.h'

    type(hilbert_comm), intent(in) :: comm
    real(dp), dimension(:), intent(in) :: send_data
    real(dp), dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, comm%send_count, comm%send_oft, MPI_DOUBLE_PRECISION, &
         &              recv_data, comm%recv_count, comm%recv_oft, MPI_DOUBLE_PRECISION, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine part_data_to_domain_dp
  !################################################################
  !################################################################
  subroutine domain_data_to_part_i4(comm, send_data, recv_data)
    implicit none
    include 'mpif.h'

    
    type(hilbert_comm), intent(in) :: comm
    integer, dimension(:), intent(in) :: send_data
    integer, dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, comm%recv_count, comm%recv_oft, MPI_INTEGER, &
         &              recv_data, comm%send_count, comm%send_oft, MPI_INTEGER, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine domain_data_to_part_i4
  !################################################################
  !################################################################
  subroutine domain_data_to_part_dp(comm, send_data, recv_data)
    use amr_commons,   only: dp
    implicit none
    include 'mpif.h'
    
    type(hilbert_comm), intent(in) :: comm
    real(dp), dimension(:), intent(in) :: send_data
    real(dp), dimension(:), intent(inout) :: recv_data

    integer  :: info, request
    integer  :: status(MPI_STATUS_SIZE)

    call MPI_IALLTOALLV(send_data, comm%recv_count, comm%recv_oft, MPI_DOUBLE_PRECISION, &
         &              recv_data, comm%send_count, comm%send_oft, MPI_DOUBLE_PRECISION, &
         &              MPI_COMM_WORLD, request, info)

    ! Finish communication (CAN BE MOVED OUTSIDE BY PASSING
    ! REQUEST HANDLE OUT OF SUBROUTINE)
    call MPI_WAIT(request, status, info)

  end subroutine domain_data_to_part_dp
!################################################################
!################################################################
!################################################################
!################################################################
#endif
end module particle_communication
