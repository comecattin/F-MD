module forces_module
    implicit none
    real(8), parameter :: epsilon = 1.0 ! LJ potential well depth
    real(8), parameter :: sigma = 1.0 ! LJ distance where potential is zero
    
contains

    subroutine compute_forces(pos, frc, n)
        real(8), dimension(n, 3), intent(in) :: pos
        real(8), dimension(n, 3), intent(out) :: frc
        integer, intent(in) :: n
        integer :: i, j
        real(8) :: r, r2, r6, lj_force

        ! Initialize forces to zero
        frc = 0.0

        ! Compute Leenard-Jones forces
        do i = 1, n-1
            do j = i+1, n
                r = sqrt(sum((pos(i, :) - pos(j, :))**2))
                r2 = r**2
                r6 = r2**3

                lj_force = 24.0d0 * epsilon * (2.0d0 * sigma**6 / r6**2 - sigma**12 / r6)

                ! Add forces to both atoms
                frc(i, :) = frc(i, :) + lj_force * (pos(i, :) - pos(j, :)) / r
                frc(j, :) = frc(j, :) - lj_force * (pos(i, :) - pos(j, :)) / r
            end do
        end do

    end subroutine compute_forces

end module forces_module
