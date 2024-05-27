module integration_module
    implicit none
    
contains
    subroutine integrate(pos, vel, frc, dt, n, box_length)
        real(8), dimension(n, 3), intent(inout) :: pos, vel
        real(8), dimension(n, 3), intent(in) :: frc
        real(8), intent(in) :: dt, box_length
        integer, intent(in) :: n
        integer :: i, j

        do i = 1, n
            do j = 1, 3
                pos(i, j) = pos(i, j) + vel(i, j) * dt + 0.5 * frc(i, j) * dt**2
                vel(i, j) = vel(i, j) + 0.5 * frc(i, j) * dt
                
                ! Apply periodic boundary conditions
                if (pos(i, j) < 0.0d0) then
                    pos(i, j) = pos(i, j) + box_length
                else if (pos(i, j) >= box_length) then
                    pos(i, j) = pos(i, j) - box_length
                end if

            end do
        end do

    end subroutine integrate
    
end module integration_module