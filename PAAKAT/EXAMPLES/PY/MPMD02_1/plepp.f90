module mod_plepp 
  implicit none
  character(64)  :: token, key

  public :: plepp_simple 
  public :: plepp_get_key
  public :: init_argvs 
  public :: init_plepp 

  !=============================================================| contains |===!
  contains

!-------------------------------------------------------------------------||---!


      subroutine plepp_simple(lcomm, show)
      implicit none
      integer(4),    intent(inout)  :: lcomm
      integer(4),    intent(in   )  :: show
      !
      integer(4)     :: iarg
      integer(4)     :: dummy = -1, lrank = -1, error = -1 
      ! 
      call plepp_init_argvs()
      do iarg = 1, command_argument_count()
        token = ""
        call get_command_argument(iarg, token)
        call plepp_set_argvs(trim(token), len_trim(token))
      enddo
      ! 
      key   = "??"
      token = "--namei"
      call plepp_analyse_argvs(trim(token), len_trim(token),0)!1=print
      call plepp_get_argvs(key, show); 
      !print *, "[init_argvs] --name:",  "'"//trim(key)//"'"
      !
      call plepp_split(lcomm, dummy, 0) ! 1 = print   
      call MPI_Comm_rank(dummy, lrank, error)  
      lcomm = dummy  
      !  
      end subroutine


      subroutine plepp_get_key( string )
      implicit none
      character(64),  intent(inout)  :: string 

      string = key 

      end subroutine


      subroutine init_plepp(typeij, namei, namej, dimij, lcomm)
      !
      ! SEE :
      !      /Users/poderozita/z2020_2/REPOSITORY/ALYAs/JUL08_3/Sources
      ! 
      implicit none
      !
      character(64), intent(in   )   :: typeij, namei, namej
      integer(4),    intent(in   )   :: dimij
      integer(4),    intent(inout)   :: lcomm
      integer        :: lrank = -1, lsize = -1, error = -1, dummy  
      !
      call plepp_split(lcomm, dummy)
      !
!#ifdef PLEPP 
       call plepp_create(trim(typeij), len_trim(typeij), &
                         trim(namei ), len_trim(namei ), &
                         trim(namej ), len_trim(namej ), &
                         dimij, lcomm)
      !
      call plepp_test() !> iam_hemelb.x : EXAMPLES/DUMMY/dummy.cxx -> hemelbCoupling->Test(...);    
      ! 
      call MPI_Comm_rank(lcomm, lrank, error)
!#endif
      ! 
      end subroutine


      subroutine init_argvs(namei, namej, show)
      implicit none
      character(64), intent(inout)  :: namei, namej
      integer(4),    intent(in   )  :: show
      !
!      character(64)  :: token, key
      integer(4)     :: iarg
      ! 
      call plepp_init_argvs()
      do iarg = 1, command_argument_count()
        token = ""
        call get_command_argument(iarg, token)
        call plepp_set_argvs(trim(token), len_trim(token))
        !
      enddo
      ! 
      key   = "??"
      token = "--namei"
      call plepp_analyse_argvs(trim(token), len_trim(token))
      call plepp_get_argvs(key, 0); ! 1=print 
      !print *, "[init_argvs] --name:",  "'"//trim(key)//"'"
      !!  
      if(key(1:5)=="HEART") then
        namei  = "ALYA"
        namej  = "HEMELB"
      endif
      if(key(1:4)=="BODY") then
        namej  = "ALYA"
        namei  = "HEMELB"
      endif
      !
      end subroutine

!-------------------------------------------------------------------------||---!
end module
