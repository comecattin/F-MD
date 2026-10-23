program test_physics
    use forces_module, only: compute_forces
    use energies_module, only: compute_potential_energy, compute_energies
    use integration_module, only: integrate
    implicit none
    integer, parameter :: n = 6
    real(8), parameter :: box = 12.0d0, h = 1.0d-5
    real(8) :: pos(n,3), shifted(n,3), vel(n,3), forces(n,3), charges(n), mass(n)
    real(8) :: plus, minus, numerical, max_error, ke, pe, total, initial, drift(2), dt
    integer :: i, j, scenario, run, step

    ! Distorted, non-planar waters; the first straddles the periodic boundary.
    pos(1,:) = [11.7d0, 2.0d0, 3.0d0]
    pos(2,:) = [0.8d0, 2.1d0, 3.2d0]
    pos(3,:) = [11.3d0, 2.9d0, 2.8d0]
    pos(4,:) = [3.5d0, 3.0d0, 4.0d0]
    pos(5,:) = [4.4d0, 3.3d0, 4.1d0]
    pos(6,:) = [3.2d0, 3.9d0, 4.3d0]
    mass = [15.999d0, 1.00784d0, 1.00784d0, 15.999d0, 1.00784d0, 1.00784d0]

    do scenario = 1, 2
        charges = 0.0d0
        if (scenario == 2) charges = [-0.82d0, 0.41d0, 0.41d0, -0.82d0, 0.41d0, 0.41d0]
        call compute_forces(pos, charges, forces, n, box)
        max_error = 0.0d0
        do i = 1, n
            do j = 1, 3
                shifted = pos
                shifted(i,j) = pos(i,j) + h
                call compute_potential_energy(shifted, charges, n, box, plus)
                shifted(i,j) = pos(i,j) - h
                call compute_potential_energy(shifted, charges, n, box, minus)
                numerical = -(plus-minus)/(2.0d0*h)
                max_error = max(max_error, abs(forces(i,j)-numerical))
            end do
        end do
        print *, 'Gradient scenario, maximum absolute error:', scenario, max_error
        if (.not. (max_error < 2.0d-6)) error stop 'Force does not match energy gradient'
        if (maxval(abs(sum(forces, dim=1))) > 1.0d-10) error stop 'Net internal force is nonzero'
        ! Translate and wrap all atoms to exercise minimum-image invariance.
        shifted = modulo(pos + 5.0d0, box)
        call compute_potential_energy(pos, charges, n, box, plus)
        call compute_potential_energy(shifted, charges, n, box, minus)
        if (abs(plus-minus) > 1.0d-10) error stop 'Energy changes on translation'
        call compute_forces(shifted, charges, vel, n, box)
        if (maxval(abs(vel-forces)) > 1.0d-10) error stop 'Force changes on translation'
    end do

    ! Same physical duration at two timesteps: Verlet energy error should scale as dt**2.
    do run = 1, 2
        shifted = pos
        vel = 0.0d0
        dt = 5.0d-6 / run ! ps
        call compute_energies(shifted, vel, mass, charges, n, box, ke, pe, initial)
        drift(run) = 0.0d0
        do step = 1, 1000 * run
            call compute_forces(shifted, charges, forces, n, box)
            call integrate(shifted, vel, forces, dt, n, box, mass, charges)
            call compute_energies(shifted, vel, mass, charges, n, box, ke, pe, total)
            drift(run) = max(drift(run), abs(total-initial))
        end do
        print *, 'Timestep, maximum energy error:', dt, drift(run)
        if (.not. (drift(run)/abs(initial) < 1.0d-4)) error stop 'Energy drift too large'
    end do
    if (drift(2) > 0.3d0*drift(1)) error stop 'Energy error does not converge quadratically'
    print *, 'Physics regression checks passed'
end program test_physics
