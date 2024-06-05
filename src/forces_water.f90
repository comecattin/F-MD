module water_forces_module
    implicit none
    real(8), parameter :: epsilon = 0.636  ! Lennard-Jones epsilon for TIP3P
    real(8), parameter :: sigma = 3.1507   ! Lennard-Jones sigma for TIP3P
    real(8), parameter :: r_cutoff = 10.0  ! Cutoff distance
    
contains

    subroutine compute_forces(pos, charges, forces, n, box_length)
        integer, intent(in) :: n
        real(8), intent(in) :: pos(n,3), charges(n)
        real(8), intent(out) :: forces(n,3)
        real(8), intent(in) :: box_length

        real(8) r2_cutoff
        real(8) :: dx, dy, dz, r2, r6, r12, lj_potential, coulomb_potential
        real(8) :: fx, fy, fz, f_lj, f_coulomb, f_total
        integer :: i, j

        r2_cutoff = r_cutoff**2

        forces = 0.0

        do i = 1, n-1
            do j = i+1, n
                dx = pos(i,1) - pos(j,1)
                dy = pos(i,2) - pos(j,2)
                dz = pos(i,3) - pos(j,3)

                ! Minimum image convention
                dx = dx - box_length * nint(dx/box_length)
                dy = dy - box_length * nint(dy/box_length)
                dz = dz - box_length * nint(dz/box_length)

                r2 = dx**2 + dy**2 + dz**2

                if r2 < r2_cutoff then
                    r6 = (sigma**2 / r2)**3
                    r12 = r6**2
                    lj_potential = 4.0 * epsilon * (r12 - r6)

                    coulomb_potential = charges(i) * charges(j) / sqrt(r2)

                    f_lj = 48.0 * epsilon * (r12 - 0.5 * r6) / r2
                    f_coulomb = coulomb_potential / r2
                    f_total = f_lj + f_coulomb

                    fx = f_total * dx
                    fy = f_total * dy
                    fz = f_total * dz

                    forces(i,1) = forces(i,1) + fx
                    forces(i,2) = forces(i,2) + fy
                    forces(i,3) = forces(i,3) + fz

                    forces(j,1) = forces(j,1) - fx
                    forces(j,2) = forces(j,2) - fy
                    forces(j,3) = forces(j,3) - fz
                end if
            end do
        end do
    end subroutine compute_forces
    
end module water_forces_module