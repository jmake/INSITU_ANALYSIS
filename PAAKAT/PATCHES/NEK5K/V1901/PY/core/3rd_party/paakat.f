      subroutine paakat_plepp_init
      use mod_plepp
      include 'SIZE'
      include 'PARALLEL'
      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal
      character(64)  :: namei
      !
      call plepp_get_key(namei)
      call paakat_initialize(nekcomm, 5, namei)   
      !
      end


      subroutine paakat_plepp_end
      call paakat_finalize() 
      end


      subroutine paakat_plepp_process()
      use mod_plepp
      implicit none
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
      include 'SOLN'
      !
      integer dim
      character(64)  :: namei
      !
      dim = 2
      if (IF3D) dim = 3
      !
      call paakat_buildvtkgrid(xm1,ym1,zm1,lx1,ly1,lz1,lelt,dim)
      call paakat_add_scalar(t, "temperature"//char(0), size(t) )
      call paakat_coprocess(time, istep)
      !
      call plepp_get_key(namei)  
      call paakat_savemesh(namei)   
      ! 
      end
