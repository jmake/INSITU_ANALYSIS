      subroutine paakat_init
        implicit none
        include 'SIZE'
        include 'PARALLEL'
        character(len=200) :: arg
        integer :: ilen
        !
        !! SEE : core/comm_mpi.f , core/drive1.f
        call nekcatalystinitialize(iglobalcomm) !nekcomm)
      end


      subroutine paakat_end()
        implicit none
        call nekcatalystfinalize()
      end


      subroutine paakat_process()
        implicit none
        include 'SIZE'
        !include 'TOTAL'
        include 'GEOM'
        include 'INPUT'
        include 'TSTEP'
        include 'SOLN'
        include 'PARALLEL'
        integer flag, dim, npr, nt, nvx
        !
        npr = product(shape(pr))
        nt  = product(shape(t))
        nvx = product(shape(vx))
       !print*, "rank,size:", nid_global, np_global
       !print*, "rank,size:", nid_global, npr, nt, nvx, nel
       !print*, "rank,size:", nid_global, nx1, ny1, nz1, nelv, nelt
       !print*, "rank,size:", nid_global, shape(vx), nt
        !
        dim = 2
        if (IF3D) dim = 3
        call buildvtkgrid(xm1, ym1, zm1, lx1, ly1, lz1, lelt, dim)
        call add_scalar_field(t , "TEMPE"//char(0), nt)
        call add_scalar_field(vx, "VELOX"//char(0), nvx)
        call add_scalar_field(vy, "VELOY"//char(0), nvx)
        call add_scalar_field(vz, "VELOZ"//char(0), nvx)
      !!call add_scalar_field(pr, "PRESS"//char(0), npr) !! assert(ndata==nPts) :(
      !!call add_vector_field(vx, vy, vz, dim, "velocity"//char(0),nvx) !! :(
        call nekcatalystcoprocess(time, istep)
       !!call exit()
      end
