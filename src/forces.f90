module forces_module
    implicit none
    real(8), parameter :: epsilon = 0.1554253 ! LJ O-O SPC/Fw
    real(8), parameter :: sigma = 3.165492    ! LJ O-O SPC/Fw
    real(8), parameter :: pi = 3.14159265358979323846
    
contains

    subroutine compute_forces(pos, frc, n, box_length)
        real(8), dimension(n, 3), intent(in) :: pos
        real(8), dimension(n, 3), intent(out) :: frc
        integer, intent(in) :: n
        real(8), intent(in) :: box_length
        integer :: i, j
        real(8) :: r, r2, r6, lj_force
        real(8), dimension(3) :: dist

        ! Initialize forces to zero
        frc = 0.0

        ! Compute Leenard-Jones forces
        do i = 1, n-1
            do j = i+1, n
                dist = pos(i, :) - pos(j, :)
                dist = dist - nint(dist / box_length) * box_length
                r = sqrt(sum(dist**2))
                r2 = r**2
                r6 = r2**3

                lj_force = 24.0d0 * epsilon * (2.0d0 * sigma**6 / r6**2 - sigma**12 / r6)

                ! Add forces to both atoms
                frc(i, :) = frc(i, :) + lj_force * dist / r
                frc(j, :) = frc(j, :) - lj_force * dist / r
            end do
        end do

    end subroutine compute_forces

    subroutine compute_forces_water(position, charges, forces, n, box_length)
        
        real(8), dimension(n, 3), intent(in) :: position
        real(8), dimension(n), intent(in) :: charges
        real(8), dimension(n, 3), intent(out) :: forces
        integer, intent(in) :: n
        real(8), intent(in) :: box_length
        integer :: i, j
        real(8), dimension(3) :: dist, vec1, vec2
        real(8) :: r, r2, f_yukawa, force_lj, force_harmonic_bond, force_harmonic_angle
        integer :: bond_type
        real(8) :: cos_theta, theta

        real(8) :: alpha, g
        real(8) :: ka, kb, r_eq, angle_eq

        logical :: is_OO, is_OH, is_HH

        alpha = 1.0
        g = 1.0
        r_eq = 1.012                     ! Angstrom
        angle_eq = 113.24                ! Deg
        angle_eq = angle_eq * pi / 180.0 ! Rad
        ka = 75.90
        kb = 1059.162

        forces = 0.0

        do i = 1, n-1

            ! Harmonic potential on HOH angle
            if (mod(i, 3) == 1) then
                vec1 = position(i+1, :) - position(i, :)
                vec2 = position(i+2, :) - position(i, :)
                
                cos_theta = dot_product(vec1, vec2) / (norm2(vec1) * norm2(vec2))
                theta = acos(cos_theta)
                
                call harmonic(theta, angle_eq, ka, force_harmonic_angle)

                forces(i, :) = forces(i, :) - force_harmonic_angle * (vec1/norm2(vec1) + vec2/norm2(vec2))
                forces(i + 1, :) = forces(i + 1, :) - force_harmonic_angle * vec1/norm2(vec1)
                forces(i + 2, :) = forces(i + 2, :) - force_harmonic_angle * vec2/norm2(vec2)
            end if

            do j = i+1, n

                bond_type = mod(i, 3) + mod(j, 3)
                bond_type = mod(bond_type, 3)

                is_OH = (bond_type == 1 .or. bond_type == 0)
                is_OO = (bond_type == 2 .and. mod(i, 3) == 1)
                is_HH = (bond_type == 2 .and. mod(i, 3) /= 1)
                
                dist = position(i, :) - position(j, :)
                r2 = sum(dist**2)
                r = sqrt(r2)

                ! Oxygen Oxygen interaction
                if (is_OO) then
                    ! Lennard Jones forces
                    call lennard_jones(r, force_lj)
                    force_lj = 0.0
                    forces(i, :) = forces(i, :) + force_lj * dist/r
                    forces(j, :) = forces(j, :) - force_lj * dist/r

                ! Oxygen Hydrogen interaction
                else if (is_OH) then
                    ! Harmonic bond forces
                    call harmonic(r, r_eq, kb, force_harmonic_bond)
                    forces(i, :) = forces(i, :) + force_harmonic_bond * dist/r
                    forces(j, :) = forces(j, :) - force_harmonic_bond * dist/r
                
                end if

                ! Coulomb force for all pair
                call coulomb_yukawa(r, 1.0_8, charges(i), charges(j), alpha, g, f_yukawa)
                forces(i, :) = forces(i, :) + f_yukawa * dist/r
                forces(j, :) = forces(j, :) - f_yukawa * dist/r

            end do
        end do

    end subroutine compute_forces_water

    subroutine coulomb_yukawa(r, mass, q1, q2, alpha, g, f_yukawa)

        real(8), intent(in) :: r
        real(8), intent(in) :: q1, q2
        real(8), intent(in) :: mass, alpha, g
        real(8), intent(out) :: f_yukawa

        f_yukawa = q1*q2 * exp(-alpha * mass * r) * (alpha * mass * r + 1)/r**2

    end subroutine coulomb_yukawa

    subroutine lennard_jones(r, force_lj)
    
        real(8), intent(in) :: r
        real(8) :: r6
        real(8), intent(out) :: force_lj

        r6 = r**6
        
        force_lj = 24.0d0 * epsilon * (2.0d0 * sigma**6 / r6**2 - sigma**12 / r6)

    end subroutine lennard_jones

    subroutine harmonic(r, r_eq, kb, force_harmonic)
    
        real(8), intent(in) :: r
        real(8), intent(in) :: r_eq, kb
        real(8), intent(out) :: force_harmonic

        force_harmonic = kb * (r - r_eq)

    end subroutine harmonic

end module forces_module
