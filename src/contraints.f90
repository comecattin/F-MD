module contraints_module
    implicit none
    
    contains

    subroutine shake(position, old_position, vel, box_length, n, dt, tolerance, max_iter)

        real(8), dimension(n, 3), intent(inout) :: position, vel
        real(8), dimension(n, 3), intent(in) :: old_position
        real(8), intent(in) :: box_length, dt, tolerance
        integer, intent(in) :: n, max_iter
        integer :: i, j, iter
        real(8) :: dx, dy, dz, r, r_old, diff, corr
        real(8) :: oh_distance

        oh_distance = 0.9572d0

        ! Iterate to apply constraints
        do iter = 1, max_iter
            do i = 1, n, 3
                ! Contraints for O-H bonds
                do j = i+1, i+2
                    ! Current bond
                    dx = position(i, 1) - position(j, 1)
                    dy = position(i, 2) - position(j, 2)
                    dz = position(i, 3) - position(j, 3)
                    ! Minimum image convention
                    dx = dx - box_length * nint(dx / box_length)
                    dy = dy - box_length * nint(dy / box_length)
                    dz = dz - box_length * nint(dz / box_length)

                    r = sqrt(dx**2 + dy**2 + dz**2)

                    ! Old bond
                    dx = old_position(i, 1) - old_position(j, 1)
                    dy = old_position(i, 2) - old_position(j, 2)
                    dz = old_position(i, 3) - old_position(j, 3)
                    ! Minimum image convention
                    dx = dx - box_length * nint(dx / box_length)
                    dy = dy - box_length * nint(dy / box_length)
                    dz = dz - box_length * nint(dz / box_length)

                    r_old = sqrt(dx**2 + dy**2 + dz**2)

                    ! Correction factor
                    diff = r - oh_distance
                    corr = diff / r

                    ! Apply correction
                    position(i, 1) = position(i, 1) - 0.5 * corr * dx
                    position(i, 2) = position(i, 2) - 0.5 * corr * dy
                    position(i, 3) = position(i, 3) - 0.5 * corr * dz
                    position(j, 1) = position(j, 1) + 0.5 * corr * dx
                    position(j, 2) = position(j, 2) + 0.5 * corr * dy
                    position(j, 3) = position(j, 3) + 0.5 * corr * dz

                    ! Update velocities
                    vel(i, 1) = (position(i, 1) - old_position(i, 1)) / dt
                    vel(i, 2) = (position(i, 2) - old_position(i, 2)) / dt
                    vel(i, 3) = (position(i, 3) - old_position(i, 3)) / dt
                    vel(j, 1) = (position(j, 1) - old_position(j, 1)) / dt
                    vel(j, 2) = (position(j, 2) - old_position(j, 2)) / dt
                    vel(j, 3) = (position(j, 3) - old_position(j, 3)) / dt

                    ! Check convergence
                    if (abs(diff) < tolerance) exit
                end do
            end do
        end do
    end subroutine shake
    
end module contraints_module