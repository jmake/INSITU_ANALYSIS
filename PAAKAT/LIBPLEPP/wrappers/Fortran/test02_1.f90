program main
  use iso_c_binding
  implicit none
  include 'mpif.h'

  interface
    function commdom_locator_get_pointer2() bind(c)
      import :: c_ptr
      implicit none
      type(c_ptr) :: commdom_locator_get_pointer2 
    end function 
  end interface

  real(c_double), pointer :: f_p, f_p2(:)
  class(*      ), pointer :: ptr_class  ! point to any type?! 
  type(c_ptr)             :: ptr_plepp1  ! -> void**   
!  type(c_ptr), value      :: ptr_plepp2  ! -> void*   :(  



  integer        :: error  = -1
  integer        :: gcomm, grank = -1, gsize = -1
  integer        :: lcomm, lrank = -1, lsize = -1

  integer(4)     :: commij
  integer        :: plepp_type, plepp_dim
  character(64)  :: app_type, namei, namej
  real(8)        :: dummy ! :) 8 or c_double  

  ! 1.0. mpi init 
  commij  = mpi_comm_null
  lcomm   = mpi_comm_null
  gcomm   = MPI_COMM_WORLD
  call MPI_Init(error)
  call MPI_Comm_rank(gcomm, grank, error);
  call MPI_Comm_size(gcomm, gsize, error);

  ! 2.0. plepp init
  plepp_type = 4  
  plepp_dim  = 2  
  app_type   = "INSITU"
  namei      = "NEK5K"
  namej      = "EXPLORER"
  lcomm      = gcomm  

  call plepp_create(trim(app_type), len_trim(app_type), &
                    trim(namei),    len_trim(namei),    &
                    trim(namej),    len_trim(namej),    &
                    plepp_dim, lcomm)

!!call plepp_getcomij(commij)  

  call MPI_Comm_rank(lcomm, lrank, error) 
  !print *, lrank, grank   

  call MPI_Finalize(error)
  !if(grank==0) 
  !print *, "gsize:", gsize, ", OK!!"

end program
