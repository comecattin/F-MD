module energies_module
    implicit none
    real(8), parameter :: epsilon = 1.0 ! LJ potential well depth
    real(8), parameter :: sigma = 1.0 ! LJ distance where potential is zero
    
contains 
    subroutine compute_energies(pos, vel, n, ke, pe, te)
        real(8), dimension(n, 3), intent(in) :: pos, vel
        integer, intent(in) :: n
        real(8), intent(out) :: ke, pe, te

        call compute_kinetic_energy(vel, n, ke)
        call compute_potential_energy(pos, n, pe)
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
        integer :: i, j
        real(8) :: r

        pe = 0.0
        do i = 1, n
            do j = i + 1, n
                r = sqrt(sum((pos(i, :) - pos(j, :)) ** 2))
                pe = pe + 4.0 * epsilon * ((sigma / r) ** 12 - (sigma / r) ** 6)
            end do
        end do
    end subroutine compute_potential_energy   
end module energies_module