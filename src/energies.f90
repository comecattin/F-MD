module energies_module
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
    subroutine compute_energies(pos, vel, charges, n, ke, pe, te)
        real(8), dimension(n, 3), intent(in) :: pos, vel
        real(8), dimension(n), intent(in) :: charges
        integer, intent(in) :: n
        real(8), intent(out) :: ke, pe, te

        call compute_kinetic_energy(vel, n, ke)
        ! call compute_potential_energy(pos, n, pe)
        call compute_potential_energy_water(pos, n, charges, pe)
        te = ke + pe

        print *, "Kinetic energy: ", ke
        print *, "Potential energy: ", pe
        print *, "Total energy: ", te
    end subroutine compute_energies

    subroutine compute_kinetic_energy(vel, n, ke)
        real(8), dimension(n, 3), intent(in) :: vel
        integer, intent(in) :: n
        real(8), intent(out) :: ke
        integer :: i

        ke = 0.0
        do i = 1, n
            ke = ke + 0.5 * sum(vel(i, :) ** 2)
        end do
    end subroutine compute_kinetic_energy

    subroutine compute_potential_energy(pos, n, pe)
        real(8), dimension(n, 3), intent(in) :: pos
        integer, intent(in) :: n
        real(8), intent(out) :: pe
        real(8) :: harmonic_p_bond, harmonic_p_angle
        integer :: i, j
        real(8) :: r

        pe = 0.0

        ! Intra molecular potential
        do i = 1, n-1, 3
            ! Harmonic potential on HOH angle
            vec1 = position(i+1, :) - position(i, :)
            vec2 = position(i+2, :) - position(i, :)
            ! Minimum image convention
            vec1 = vec1 - nint(vec1 / box_length) * box_length
            vec2 = vec2 - nint(vec2 / box_length) * box_length
            cos_theta = dot_product(vec1, vec2) / (norm2(vec1) * norm2(vec2))
            theta = acos(cos_theta)
            call harmonic_potential(theta, angle_eq_rad, ka_HOH, harmonic_p_angle)
            p = p + harmonic_p_angle

            ! Harmonic potential on O-H bond
            !    O-H1 bond
            call harmonic_potential(norm2(vec1), r_eq_OH, kb_OH, harmonic_p_bond)
            p = p + harmonic_p_bond
            !    O-H2 bond
            call harmonic_potential(norm2(vec2), r_eq_OH, kb_OH, harmonic_p_bond)
            p = p + harmonic_p_bond
        end do

        ! Lennard-Jones potential
        do i = 1, n
            do j = i + 1, n
                r = sqrt(sum((pos(i, :) - pos(j, :)) ** 2))
                pe = pe + 4.0 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
            end do
        end do
    end subroutine compute_potential_energy

    subroutine yukawa_potential(r, mass, q1, q2, alpha_param, yukawa_p)
        
        real(8), intent(in) :: r
        real(8), intent(in) :: q1, q2
        real(8), intent(in) :: mass, alpha_param
        real(8), intent(out) :: yukawa_p

        yukawa_p = - q1*q2 * exp(-alpha_param * mass * r) / r

    end subroutine yukawa_potential

    subroutine lennard_jones_potential(r, lj_p)
    
        real(8), intent(in) :: r
        real(8), intent(out) :: lj_p

        lj_p = 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)

    end subroutine lennard_jones_potential

    subroutine harmonic_potential(r, r_eq, kb, harmonic_p)
    
        real(8), intent(in) :: r
        real(8), intent(in) :: r_eq, kb
        real(8), intent(out) :: harmonic_p

        harmonic_p = 0.5 * kb * (r - r_eq)**2

    end subroutine harmonic_potential

    
end module energies_module