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
    real(8), allocatable :: K(:,:,:)
    integer n 
    points = readCsv('point.csv', 3)
    elements = makeElements(points)
    n = size(elements, 1)
    allocate(K(n,6,6))
    do i = 1, n
        K(i,:,:) = makeK(elements(i,:,:))
    enddo
    !call printArray2(points)
    !call printArray3(elements)
    call printArray3(K)
end program main

