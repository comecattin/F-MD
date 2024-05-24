module forces_module
    implicit none
    
contains

    subroutine forces(pos, frc, n)
        real(8), dimension(n, 3), intent(in) :: pos
        real(8), dimension(n, 3), intent(out) :: frc
        integer, intent(in) :: n
        integer :: i

        ! Compute the forces
        do i = 1, n
            frc(i, :) = 0.0
        end do
    end subroutine forces

end module forces_module
