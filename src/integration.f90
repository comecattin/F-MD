module integration_module
    implicit none
    
contains
    subroutine integrate(pos, vel, frc, dt, n)
        real(8), dimension(n, 3), intent(inout) :: pos, vel
        real(8), dimension(n, 3), intent(in) :: frc
        real(8), intent(in) :: dt
        integer, intent(in) :: n
        integer :: i

        ! Integrate using velocity verlet
        do i = 1, n
            pos(i, :) = pos(i, :) + vel(i, :) * dt + 0.5 * frc(i, :) * dt**2
            vel(i, :) = vel(i, :) + 0.5 * frc(i, :) * dt
        end do

    end subroutine integrate
    
end module integration_module