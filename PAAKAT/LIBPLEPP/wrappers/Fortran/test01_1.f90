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



  integer           error
  integer           world_comm, world_rank, world_size
  integer(4)     :: commi, commij
  integer        :: plepp_type, plepp_dim
  character(64)  :: app_type, namei, namej
  real(8)        :: dummy ! :) 8 or c_double  

  ! 1.0. mpi init 
  world_comm = MPI_COMM_WORLD
  call MPI_Init(error)
  call MPI_Comm_rank(world_comm, world_rank, error);
  call MPI_Comm_size(world_comm, world_size, error);

  ! 2.0. plepp init
  plepp_type = 4  
  plepp_dim  = 2  
  app_type   = "INSITU"
  namei      = "NEK5K"
  namej      = "EXPLORER"
  call commdom_create()
  call commdom_set_names(trim(app_type), len_trim(app_type), trim(namei), len_trim(namei))

  commi      = mpi_comm_null
  commij     = mpi_comm_null
  call commdom_create_commij(world_comm, commi)
  call commdom_get_commij(trim(namej), len_trim(namej), commij)


  !call commdom_delete() 
  call MPI_Finalize(error)
  !if(world_rank==0) 
  print *, "world_size:", world_size, ", OK!!"

end program
