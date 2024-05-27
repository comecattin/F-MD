module initialization_module
    implicit none
    
contains
    subroutine initialize(pos, vel, n, box_length)

        real(8), dimension(n,3), intent(out) :: pos, vel
        integer, intent(in) :: n
        real(8), intent(in) :: box_length
        integer :: i, j

        call random_seed() ! Initialize the random number generator
        call random_number(pos) ! Initialize the positions between 0 and 1
        pos = pos * box_length

        call random_number(vel)
        vel = vel - 0.5 ! Shift the velocities to be between -0.5 and 0.5

    end subroutine initialize
    
end module initialization_module