module integration_module
    use contraints_module, only: shake
    use forces_module, only: compute_forces_water
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

    subroutine integrate_constraints(pos, vel, forces, charges, n, dt, box_length, tolerance, max_iter)
        implicit none
        integer, intent(in) :: n, max_iter
        real(8), intent(in) :: dt, box_length, tolerance
        real(8), dimension(n, 3), intent(inout) :: pos, vel, forces
        real(8), dimension(n, 3) :: old_pos
        real(8), dimension(n) :: charges
        integer :: i

        ! Save old positions
        old_pos = pos

        ! Velocity Verlet integration
        do i = 1, n
            vel(i, :) = vel(i, :) + 0.5d0 * dt * forces(i, :)
            pos(i, :) = pos(i, :) + dt * vel(i, :)
        end do

        ! Apply periodic boundary conditions
        pos = pos - box_length * floor(pos / box_length)

        ! Apply constraints using SHAKE
        call shake(pos, old_pos, vel, box_length, n, dt, tolerance, max_iter)

        ! Compute forces at new positions (not shown here, but call your force calculation subroutine)
        call compute_forces_water(pos, charges, forces, n, box_length)

        ! Update velocities
        do i = 1, n
            vel(i, :) = vel(i, :) + 0.5d0 * dt * forces(i, :)
        end do
    end subroutine integrate_constraints
    
end module integration_module