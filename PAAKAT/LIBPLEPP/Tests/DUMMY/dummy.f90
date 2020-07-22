program main
  implicit none
  include 'mpif.h'

  ! 0.0. 
  integer        :: error  = -1
  integer        :: gcomm, grank = -1, gsize = -1
  integer        :: lcomm, lrank = -1, lsize = -1

  integer        :: plepp_dim, dummy, timeSteps, itime, ii  
  character(64)  :: app_type, namei, namej

  real(8) :: send(2), recv(2) 

  logical :: test = .False. 

  ! 1.0. mpi init 
  lcomm   = mpi_comm_null
  gcomm   = MPI_COMM_WORLD
  call MPI_Init(error)
  call MPI_Comm_rank(gcomm, grank, error);
  call MPI_Comm_size(gcomm, gsize, error);

  ! 2.0. init_plepp
  plepp_dim  = 0   
  app_type   = "PLEPP"
  namei      = "ALYA"
  namej      = "HEMELB"
  lcomm      = gcomm  

  call plepp_create(trim(app_type), len_trim(app_type), &
                    trim(namei),    len_trim(namei),    &
                    trim(namej),    len_trim(namej),    &
                    plepp_dim, lcomm)
  call plepp_test() 

  ! 3.0. test_plepp     
  dummy     = 3
  timeSteps = -1  
  call plepp_sendrecv_int(dummy, 1, timeSteps, 1, 0) 

  recv = huge(1.0) 
  send = (/2.0,1.0/)
  do itime = 1,timeSteps 
    do ii = 1,size(send)  
      call plepp_exchange( send(ii), recv(ii), 0)
    enddo  
  enddo 

  test = all( recv==(/-1.0,-2.0/) )
  if(.not.test) then
    print *, "Test ERROR : recv==(/-1.0,-2.0/) !!" 
    call exit(0)
  endif 

  call MPI_Comm_rank(lcomm, lrank, error);
  if(lrank==0)  print *, "       ["//trim(namei)//"] recv:", recv

  call MPI_Finalize(error)

end program
