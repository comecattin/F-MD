module initialization_module
    implicit none
    
contains
    subroutine initialize(pos, vel, n)

        real(8), dimension(n,3), intent(out) :: pos, vel
        integer, intent(in) :: n
        integer :: i, j

        ! Initialize the positions and velocities between -1 and 1
        do i = 1, n
            do j = 1, 3
                pos(i,j) = 2.0 * (rand() - 0.5)
                vel(i,j) = 2.0 * (rand() - 0.5)
            end do
        end do

    end subroutine initialize
    
end module initialization_module