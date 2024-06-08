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
        real(8), dimension(3) :: dist
        real(8) :: r, r2, coulomb_force

        call compute_forces(position, forces, n, box_length)

        ! Add Coulomb forces
        do i = 1, n-1
            do j = i+1, n
                
                dist = position(i, :) - position(j, :)
                r2 = sum(dist**2)
                r = sqrt(r2)

                coulomb_force = charges(i) * charges(j) / r2

                forces(i, :) = forces(i, :) + coulomb_force * dist / r
                forces(j, :) = forces(j, :) - coulomb_force * dist / r  

            end do
        end do

    end subroutine compute_forces_water

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
