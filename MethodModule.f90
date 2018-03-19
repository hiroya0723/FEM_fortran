module MethodModule
    implicit none
contains
    subroutine printArray2(array)
        real(8) array(:,:)
        integer i, n
        n = size(array, 1)
        do i = 1, n
            print *, array(i, :)
        enddo
    end subroutine printArray2

    subroutine printArray3(array)
        real(8) array(:,:,:)
        integer i, j, n, m
        n = size(array, 1)
        m = size(array, 2)
        do i = 1, n
            print *, i
            do j = 1, m
                print *, array(i, j, :)
            enddo
        enddo
    end subroutine printArray3
    
    function readCsv(filename, m) result(output)
        integer m, n 
        integer i
        character(*) filename
        real(8), allocatable :: output(:,:)
        open (17, file=filename, status='old')
        ! === レコード数を調べる ===
        n = 0
        read (17, '()')
        do
          read (17, *, end=100)! ファイル終端ならば999に飛ぶ
          n = n + 1
        end do
      100 continue
        allocate(output(n,m))
        print *, n
        rewind (17)  ! ファイルの最初に戻る
        read (17, '()')
        do i = 1, n
          read (17, *) output(i,:)
        end do
        close (17)
    end function readCsv
end module MethodModule
