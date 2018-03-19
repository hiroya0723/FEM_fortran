module ElementModule
    use MethodModule
    implicit none
    real(8) :: nu = 1/3.0d0
    real(8) :: E  = 100
    real(8) :: h  = 1.0d0
contains
    function makeElements(points) result(elements)
        integer, allocatable :: index(:,:)
        real(8), allocatable :: elements(:,:,:)
        real(8) points(:,:)
        integer i, j, n
        
        index = readCsv('element.csv', 3)
        n = size(index, 1)
        allocate(elements(n,3,3))
        do i = 1, n
            do j = 1, 3
                elements(i,j,:) = points(index(i,j),:)
            enddo
        enddo
    end function makeElements
        
    function getL(point) result(l)
        real(8) point(3, 3)
        integer i
        real(8) l(3)
        
        do i = 1, 3
            l(i) = (point((mod(i,3)+1), 2)-point((mod(i+1,3)+1), 2))**2 + &
                   (point((mod(i,3)+1), 3)-point((mod(i+1,3)+1), 3))**2
            l(i) = sqrt(l(i))
        enddo
    end function getL
    
    function makeA(l) result(A)
        real(8) l(3)
        real(8) A
        real(8) s
        integer i
        
        s = 0.0d0
        do i = 1, 3
            s = s + l(i)
        enddo
        s = s/2.0d0
        A = s
        do i = 1, 3
            A = A * (s-l(i))
        enddo
        A = sqrt(A)
    end function makeA
    
    function makeB(point) result(B)
        real(8) point(3, 3)
        real(8) A
        real(8) l(3)
        real(8) B(3, 6)
        integer i, j
        
        l = getL(point)
        A = makeA(l)
        B(1,1) = point(2,3)-point(3,3)
        B(1,2) = 0.0d0
        B(1,3) = point(3,3)-point(1,3)
        B(1,4) = 0.0d0
        B(1,5) = point(1,3)-point(2,3)
        B(1,6) = 0.0d0
        B(2,1) = 0.0d0
        B(2,2) = point(3,2)-point(2,2)
        B(2,3) = 0.0d0
        B(2,4) = point(1,2)-point(3,2)
        B(2,5) = 0.0d0
        B(2,6) = point(2,2)-point(1,2)
        B(3,1) = point(3,2)-point(2,2)
        B(3,2) = point(2,3)-point(3,3)
        B(3,3) = point(1,2)-point(3,2)
        B(3,4) = point(3,3)-point(1,3)
        B(3,5) = point(2,2)-point(1,2)
        B(3,6) = point(1,3)-point(2,3)
        B = B/(2.0d0*A)
    end function makeB
    
    function makeD() result(D)
        real(8) D(3, 3)
        
        D(1, 1) = 1.0d0
        D(1, 2) = nu
        D(1, 3) = 0.0d0
        D(2, 1) = nu
        D(2, 2) = 1.0d0
        D(2, 3) = 0.0d0
        D(3, 1) = 0.0d0
        D(3, 2) = 0.0d0
        D(3, 3) = (1-nu)/2.0d0
        D = E*D/(1-nu**2)
    end function makeD
    
    function makeK(point, num) result(K)
        real(8) point(3, 3)
        real(8) A 
        real(8) B(3, 6), B_t(6, 3)
        real(8) D(3, 3)
        real(8) K_element(6, 6)
        real(8), allocatable :: K(:,:)
        real(8) l(3)
        integer i, j
        integer num
        integer index(6)
        
        allocate(K(num, num))
        
        do i = 1, num
            do j = 1, num
                K(i,j) = 0.0d0
            enddo
        enddo
        l = getL(point)
        A = makeA(l)
        B = makeB(point)
        B_t = transpose(B)
        D = makeD()
        
        K_element = matmul(matmul(B_t, D), B)
        K_element = A * h * K_element
        
        do i = 1, 3
            index(2*(i-1)+1) = 2*point(i,1)-1
            index(2*i)       = 2*point(i,1)
        enddo
        
        do i = 1, 6
            do j = 1, 6
                K(index(i), index(j)) = K_element(i,j)
            enddo
        enddo
    end function makeK
end module ElementModule