program main
  use iso_c_binding
  !use mpi 
  implicit none
  include 'mpif.h'

  integer      error
  integer      world_comm, world_rank, world_size

  world_comm = MPI_COMM_WORLD

  call MPI_Init(error)
  call MPI_Comm_rank(world_comm, world_rank, error);
  call MPI_Comm_size(world_comm, world_size, error);

  call nekcatalystinitialize(world_comm)
  call nekcatalystfinalize()

  call MPI_Finalize(error)
  if(world_rank==0) print *, "world_size:", world_size, ", OK!!"

end program
