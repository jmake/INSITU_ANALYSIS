      subroutine in_situ_init()
      call paakat_init()
      end

      subroutine in_situ_check()
      call paakat_process()
      end

      subroutine in_situ_end()
      call paakat_end()
      end
