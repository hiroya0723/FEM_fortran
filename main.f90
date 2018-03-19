program main
    use ElementModule
    use MethodModule
    
    implicit none
    
    real(8), allocatable :: points(:,:)
    real(8), allocatable :: elements(:,:,:)
    integer i, j
    real(8) A
    real(8) resB(3, 6)
    real(8) resD(3, 3)
    real(8) l(3)
    real(8), allocatable :: K(:,:)
    integer n, m 
    
    points = readCsv('point.csv', 3)
    elements = makeElements(points)
    n = size(points, 1) * 2
    m = size(elements, 1)
    allocate(K(n, n))
    do i = 1, n
        do j = 1, n
            K(i,j) = 0.0d0
        enddo
    enddo
    do i = 1, m
        K = K + makeK(elements(i,:,:), n)
    enddo
    call printArray2(K)
end program main

