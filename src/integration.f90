module integration_module
    use contraints_module, only: shake
    use forces_module, only: compute_forces_water
    implicit none
    
contains
    subroutine integrate(pos, vel, frc, dt, n, box_length, mass, charges)
        real(8), dimension(n, 3), intent(inout) :: pos, vel
        real(8), dimension(n, 3), intent(in) :: frc
        real(8), intent(in) :: dt, box_length
        integer, intent(in) :: n
        real(8), dimension(n), intent(in) :: mass, charges
        integer :: i, j
        real(8), dimension(n, 3):: a, a_new, frc_new

        a = frc / spread(mass, 2, size(frc, 2))
        pos = pos + vel * dt + 0.5 * a * dt*dt
        ! Periodic bondary conditions
        pos = pos - box_length * floor(pos / box_length)
        ! New acceleration
        call compute_forces_water(pos, charges, frc_new, n, box_length)
        a_new = frc_new / spread(mass, 2, size(frc_new,2))
        vel = vel + 0.5 * (a + a_new) * dt

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
        !call shake(pos, old_pos, vel, box_length, n, dt, tolerance, max_iter)

        ! Compute forces at new positions (not shown here, but call your force calculation subroutine)
        call compute_forces_water(pos, charges, forces, n, box_length)

        ! Update velocities
        do i = 1, n
            vel(i, :) = vel(i, :) + 0.5d0 * dt * forces(i, :)
        end do
    end subroutine integrate_constraints
    
end module integration_module