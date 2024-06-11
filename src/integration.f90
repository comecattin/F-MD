module integration_module
    use forces_module, only: compute_forces
    implicit none
    
contains
    subroutine integrate(pos, vel, frc, dt, n, box_length, mass, charges)
        real(8), dimension(n, 3), intent(inout) :: pos, vel
        real(8), dimension(n, 3), intent(in) :: frc
        real(8), intent(in) :: dt, box_length
        integer, intent(in) :: n
        real(8), dimension(n), intent(in) :: mass, charges
        integer :: i, j
        real(8), dimension(n, 3):: a, a_new, frc_new

        a = frc / spread(mass, 2, size(frc, 2))
        vel = vel + 0.5 * dt * a
        pos = pos + dt * vel
        pos = pos - box_length * floor(pos / box_length)
        call compute_forces(pos, charges, frc_new, n, box_length)
        a_new = frc_new / spread(mass, 2, size(frc_new, 2))
        vel = vel + 0.5 * dt * a_new

    end subroutine integrate
    
end module integration_module