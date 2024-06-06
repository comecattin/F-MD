module contraints_module
    implicit none
    
    contains

    subroutine shake(position, old_position, vel, box_length, n, dt, tolerance, max_iter)
        real(8), dimension(n, 3), intent(inout) :: position, vel
        real(8), dimension(n, 3), intent(in) :: old_position
        real(8), intent(in) :: box_length, dt, tolerance
        integer, intent(in) :: n, max_iter
        integer :: i, j, iter
        real(8) :: dx, dy, dz, r, diff, corr, cos_theta, current_cos_theta
        real(8) :: hx1, hy1, hz1, hx2, hy2, hz2, length1, length2
        real(8) :: oh_distance, hoh_angle, cos_hoh_angle

        oh_distance = 0.9572d0
        hoh_angle = 104.52d0
        cos_hoh_angle = cos(hoh_angle * acos(-1.0d0) / 180.0d0)

        ! Iterate to apply constraints
        do iter = 1, max_iter
            corr = 0.0d0
            do i = 1, n, 3
                ! Constraints for O-H bonds
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

                    ! Correction factor
                    diff = r - oh_distance
                    corr = corr + abs(diff)
                    dx = dx / r
                    dy = dy / r
                    dz = dz / r
                    position(i, 1) = position(i, 1) - 0.5 * diff * dx
                    position(i, 2) = position(i, 2) - 0.5 * diff * dy
                    position(i, 3) = position(i, 3) - 0.5 * diff * dz
                    position(j, 1) = position(j, 1) + 0.5 * diff * dx
                    position(j, 2) = position(j, 2) + 0.5 * diff * dy
                    position(j, 3) = position(j, 3) + 0.5 * diff * dz

                    ! Update velocities
                    vel(i, 1) = (position(i, 1) - old_position(i, 1)) / dt
                    vel(i, 2) = (position(i, 2) - old_position(i, 2)) / dt
                    vel(i, 3) = (position(i, 3) - old_position(i, 3)) / dt
                    vel(j, 1) = (position(j, 1) - old_position(j, 1)) / dt
                    vel(j, 2) = (position(j, 2) - old_position(j, 2)) / dt
                    vel(j, 3) = (position(j, 3) - old_position(j, 3)) / dt
                end do

                ! Constraints for H-O-H angle
                hx1 = position(i+1, 1) - position(i, 1)
                hy1 = position(i+1, 2) - position(i, 2)
                hz1 = position(i+1, 3) - position(i, 3)
                hx2 = position(i+2, 1) - position(i, 1)
                hy2 = position(i+2, 2) - position(i, 2)
                hz2 = position(i+2, 3) - position(i, 3)

                ! Minimum image convention
                hx1 = hx1 - box_length * nint(hx1 / box_length)
                hy1 = hy1 - box_length * nint(hy1 / box_length)
                hz1 = hz1 - box_length * nint(hz1 / box_length)
                hx2 = hx2 - box_length * nint(hx2 / box_length)
                hy2 = hy2 - box_length * nint(hy2 / box_length)
                hz2 = hz2 - box_length * nint(hz2 / box_length)

                length1 = sqrt(hx1*hx1 + hy1*hy1 + hz1*hz1)
                length2 = sqrt(hx2*hx2 + hy2*hy2 + hz2*hz2)

                current_cos_theta = (hx1*hx2 + hy1*hy2 + hz1*hz2) / (length1 * length2)
                diff = current_cos_theta - cos_hoh_angle

                ! Apply corrections to satisfy the angle constraint
                corr = corr + abs(diff)
                hx1 = hx1 / length1
                hy1 = hy1 / length1
                hz1 = hz1 / length1
                hx2 = hx2 / length2
                hy2 = hy2 / length2
                hz2 = hz2 / length2
                position(i+1, :) = position(i+1, :) - 0.5 * diff * (/hx2, hy2, hz2/)
                position(i+2, :) = position(i+2, :) - 0.5 * diff * (/hx1, hy1, hz1/)
                position(i, :) = position(i, :) + 0.5 * diff * (/hx1 + hx2, hy1 + hy2, hz1 + hz2/)
            end do

            ! Check for convergence
            if (corr < tolerance) exit
        end do
    end subroutine shake
end module contraints_module
