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
        rewind (17)  ! ファイルの最初に戻る
        read (17, '()')
        do i = 1, n
          read (17, *) output(i,:)
        end do
        close (17)
    end function readCsv
    
    subroutine solve(a,b,x)
        real(8) a(:,:)
	real(8) b(:), x(:)
	real(8) temp, alpha, reserve
	integer i, j, k, L, dim, maxline

        x = 0.0d0
        dim = size(a,1)
	do k = 1, dim - 1
            maxline = k
            do i = k + 1, dim
                if (abs(a(i,k)) .gt. abs(a(maxline,k))) then
                        maxline = i
                end if

                do L = 1, dim
                        reserve = a(k,L)
                        a(k,L) = a(maxline,L)
                        a(maxline,L) = reserve
                end do

                reserve = b(k)
                b(k) = b(maxline)
                b(maxline) = reserve
            end do

            do i = k + 1, dim
                alpha = a(i,k) / a(k,k)
                do j = k , dim
                    a(i,j) = a(i,j) - alpha * a(k,j)
                end do
                b(i) = b(i) - alpha * b(k)
            end do
	end do 

	x(dim) = b(dim) / a(dim,dim)

	do k = dim-1, 1, -1
            temp = 0.0
            do j = k + 1, dim
                temp = temp + a(k,j) * x(j)
            end do
            x(k) = (b(k) - temp) / a(k,k)
	end do
        print *, x
    end subroutine solve
end module MethodModule
