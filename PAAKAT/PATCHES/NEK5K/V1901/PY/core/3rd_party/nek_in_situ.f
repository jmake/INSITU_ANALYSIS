      subroutine in_situ_init()
      call paakat_plepp_init()
      end 

      subroutine in_situ_check()
      call paakat_plepp_process()
      end

      subroutine in_situ_end()
      call paakat_plepp_end()
      end
