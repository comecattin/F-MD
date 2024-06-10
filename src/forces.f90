module forces_module
    implicit none
    real(8), parameter :: pi = 3.14159265358979323846
    real(8), parameter :: epsilon = 0.1554253          ! kcal/mol
    real(8), parameter :: sigma = 3.165492             ! Angstrom 
    real(8), parameter :: ka_HOH = 75.90               ! kcal.mol-1.rad-2
    real(8), parameter :: kb_OH = 1059.162             ! kcal.mol-1.angstrom-2
    real(8), parameter :: alpha = 1.0
    real(8), parameter :: r_eq_OH = 1.012              ! Angstrom
    real(8), parameter :: angle_eq = 113.24            ! deg
    
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
        real(8) :: angle_eq_rad
        logical :: is_OO, is_OH, is_HH

        
        angle_eq_rad = angle_eq * pi / 180.0 ! Rad
        

        forces = 0.0

        ! Intra molecular interactions
        do i = 1, n-1, 3
            ! Harmonic potential on HOH angle
            vec1 = position(i+1, :) - position(i, :)
            vec2 = position(i+2, :) - position(i, :)

            ! Minimum image convention
            vec1 = vec1 - nint(vec1 / box_length) * box_length
            vec2 = vec2 - nint(vec2 / box_length) * box_length
            
            cos_theta = dot_product(vec1, vec2) / (norm2(vec1) * norm2(vec2))
            theta = acos(cos_theta)
            
            call harmonic(theta, angle_eq_rad, ka_HOH, force_harmonic_angle)

            forces(i, :) = forces(i, :) - force_harmonic_angle * (vec1/norm2(vec1) + vec2/norm2(vec2))
            forces(i + 1, :) = forces(i + 1, :) - force_harmonic_angle * vec1/norm2(vec1)
            forces(i + 2, :) = forces(i + 2, :) - force_harmonic_angle * vec2/norm2(vec2)


            ! Harmonic potential on O-H bond
            !    O-H1 bond
            call harmonic(norm2(vec1), r_eq_OH, kb_OH, force_harmonic_bond)
            forces(i, :) = forces(i, :) + force_harmonic_bond * vec1/norm2(vec1)
            forces(i + 1, :) = forces(i + 1, :) - force_harmonic_bond * vec1/norm2(vec1)
            
            !    O-H2 bond
            call harmonic(norm2(vec2), r_eq_OH, kb_OH, force_harmonic_bond)
            forces(i, :) = forces(i, :) + force_harmonic_bond * vec2/norm2(vec2)
            forces(i + 2, :) = forces(i + 2, :) - force_harmonic_bond * vec2/norm2(vec2)
            
        end do

        ! Van der Waals interactions
        do i = 1, n-4, 3
            do j = i+3, n, 3
                dist = position(i, :) - position(j, :)
                ! Minimum image convention
                dist = dist - nint(dist / box_length) * box_length
                r = norm2(dist)
                call lennard_jones(r, force_lj)
                forces(i, :) = forces(i, :) + force_lj * dist/r
                forces(j, :) = forces(j, :) - force_lj * dist/r
            end do
        end do

        ! Coulomb interactions
        do i = 1, n-1
            do j = i+1, n
                dist = position(i, :) - position(j, :)
                ! Minimum image convention
                dist = dist - nint(dist / box_length) * box_length
                r = norm2(dist)
                call coulomb_yukawa(r, 1.0_8, charges(i), charges(j), alpha, f_yukawa)
                forces(i, :) = forces(i, :) + f_yukawa * dist/r
                forces(j, :) = forces(j, :) - f_yukawa * dist/r
            end do
        end do

    end subroutine compute_forces_water

    subroutine coulomb_yukawa(r, mass, q1, q2, alpha_param, f_yukawa)

        real(8), intent(in) :: r
        real(8), intent(in) :: q1, q2
        real(8), intent(in) :: mass, alpha_param
        real(8), intent(out) :: f_yukawa

        f_yukawa = q1*q2 * exp(-alpha_param * mass * r) * (alpha_param * mass * r + 1)/r**2

    end subroutine coulomb_yukawa

    subroutine lennard_jones(r, force_lj)
    
        real(8), intent(in) :: r
        real(8), intent(out) :: force_lj

        force_lj = 48*epsilon/sigma * ((sigma / r)**13 - 0.5 * (sigma / r)**7)

    end subroutine lennard_jones

    subroutine harmonic(r, r_eq, kb, force_harmonic)
    
        real(8), intent(in) :: r
        real(8), intent(in) :: r_eq, kb
        real(8), intent(out) :: force_harmonic

        force_harmonic = kb * (r - r_eq)

    end subroutine harmonic

end module forces_module
