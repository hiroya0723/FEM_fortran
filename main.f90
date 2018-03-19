program main
    use ElementModule
    use MethodModule
    
    implicit none
    
    real(8), allocatable :: points(:,:)
    real(8), allocatable :: elements(:,:,:)
    integer, allocatable :: condition(:,:)
    real(8), allocatable :: f(:,:)
    integer i, j
    real(8) A
    real(8) resB(3, 6)
    real(8) resD(3, 3)
    real(8) l(3)
    real(8), allocatable :: K_all(:,:)
    integer n, m 
    real(8), allocatable :: d(:)
    
    points = readCsv('point.csv', 3)
    condition = readCsv('condition.csv', 3)
    f = readCsv('force.csv', 3)
    elements = makeElements(points)
    n = size(points, 1) * 2
    m = size(elements, 1)
    allocate(K_all(n, n))
    
    d = initD(condition)
    do i = 1, n
        do j = 1, n
            K_all(i,j) = 0.0d0
        enddo
    enddo
    do i = 1, m
        K_all = K_all + makeK(elements(i,:,:), n)
    enddo
    call calcD(K_all, d, f)
    !call printArray2(K_all)
    !print *, d
end program main

