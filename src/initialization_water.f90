module initialization_water
    implicit none
    contains
    subroutine initialize(pos, vel, charges, n, box_length)
        real(8), dimension(n, 3), intent(out) :: pos, vel
        real(8), dimension(n), intent(out) :: charges
        integer, intent(in) :: n
        real(8), intent(in) :: box_length
        integer :: i, j

        call random_seed()
        call random_number(pos)
        pos = pos * box_length

        call random_number(vel)
        vel = vel - 0.5d0

        ! Initialize charges for TIP3P water model
        do i = 1, n, 3
            charges(i) = -0.834  ! Oxygen charge
            charges(i+1) = 0.417  ! Hydrogen charge
            charges(i+2) = 0.417  ! Hydrogen charge
        end do
    end subroutine initialize
end module initialization_water
