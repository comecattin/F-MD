program test_units
    use forces_module, only: compute_forces, coulomb_yukawa
    use energies_module, only: compute_kinetic_energy, yukawa_potential
    use integration_module, only: integrate
    implicit none
    real(8) :: pos(3,3), vel(3,3), charges(3), mass(3), forces(3,3), new_forces(3,3)
    real(8) :: expected_pos(3,3), expected_vel(3,3), energy, force, dt

    ! 1 g/mol at 1 A/ps: 0.5 * 1e-3 kg/mol * (100 m/s)**2 / 4184 J/kcal.
    mass = 1.0d0
    vel = 0.0d0
    vel(1,1) = 1.0d0
    call compute_kinetic_energy(vel, mass, 3, energy)
    call check(energy, 5.0d0/4184.0d0, 1.0d-14, 'Kinetic energy SI reference')

    ! Two unit charges 1 A apart, without screening: 332.063713 kcal/mol.
    call yukawa_potential(1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0, energy)
    call coulomb_yukawa(1.0d0, 1.0d0, 1.0d0, 1.0d0, 0.0d0, force)
    call check(energy, 332.06371307417d0, 1.0d-9, 'Coulomb energy reference')
    call check(force, 332.06371307417d0, 1.0d-9, 'Coulomb force reference')
    call yukawa_potential(2.0d0, 1.0d0, -1.0d0, 1.0d0, 0.5d0, energy)
    call coulomb_yukawa(2.0d0, 1.0d0, -1.0d0, 1.0d0, 0.5d0, force)
    call check(energy, -332.06371307417d0*exp(-1.0d0)/2.0d0, 1.0d-9, 'Yukawa energy reference')
    call check(force, -332.06371307417d0*exp(-1.0d0)/2.0d0, 1.0d-9, 'Yukawa force reference')

    ! Verify both velocity half-steps and the displacement in ps/A units.
    pos(1,:) = [3.0d0, 3.0d0, 3.0d0]
    pos(2,:) = [4.1d0, 3.1d0, 3.0d0]
    pos(3,:) = [2.6d0, 3.9d0, 3.2d0]
    mass = [15.999d0, 1.00784d0, 1.00784d0]
    charges = 0.0d0
    vel = 0.1d0
    dt = 1.0d-5
    call compute_forces(pos, charges, forces, 3, 12.0d0)
    ! 1 kcal/mol/A divided by 1 g/mol = 418.4 A/ps**2, derived from SI.
    expected_vel = vel + 0.5d0 * dt * 418.4d0 * forces / spread(mass, 2, 3)
    expected_pos = pos + dt * expected_vel
    call compute_forces(expected_pos, charges, new_forces, 3, 12.0d0)
    expected_vel = expected_vel + 0.5d0 * dt * 418.4d0 * new_forces / spread(mass, 2, 3)
    call integrate(pos, vel, forces, dt, 3, 12.0d0, mass, charges)
    call check(maxval(abs(pos-expected_pos)), 0.0d0, 1.0d-12, 'Acceleration displacement reference')
    call check(maxval(abs(vel-expected_vel)), 0.0d0, 1.0d-12, 'Acceleration velocity reference')
    print *, 'Physical unit checks passed'
contains
    subroutine check(actual, expected, tolerance, label)
        real(8), intent(in) :: actual, expected, tolerance
        character(len=*), intent(in) :: label
        if (.not. (abs(actual-expected) <= tolerance)) then
            print *, label, ': actual, expected = ', actual, expected
            error stop 'Physical unit check failed'
        end if
    end subroutine check
end program test_units
