program md_simulation
    use initialization, only : initialize
    implicit none

    ! Define parameters
    integer, parameter :: n_atoms = 100
    integer, parameter :: n_steps = 1000
    real(8), parameter :: dt = 0.001

    ! Define position, velocity and forces arrays
    real(8), dimension(n_atoms, 3) :: positions, velocities, forces

    ! Initialize the positions and velocities
    call initialize(positions, velocities, n_atoms)

    print *, 'Initial positions:'
    print *, positions(1:5, :)

    print *, 'Initial velocities:'
    print *, velocities(1:5, :)

    
    
end program md_simulation