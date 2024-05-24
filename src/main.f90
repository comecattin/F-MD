program md_simulation
    use initialization_module, only : initialize
    use forces_module, only : compute_forces
    implicit none

    ! Define parameters
    integer, parameter :: n_atoms = 100
    integer, parameter :: n_steps = 1000
    real(8), parameter :: dt = 0.001

    ! Define position, velocity and forces arrays
    real(8), dimension(n_atoms, 3) :: positions, velocities, forces

    ! Other variables
    integer :: step

    ! Initialize the positions and velocities
    call initialize(positions, velocities, n_atoms)

    ! Perform the molecular dynamics simulation
    do step = 1, n_steps
        ! Compute the forces
        call compute_forces(positions, forces, n_atoms)

        ! Integrate positions and velocities
    end do

    print *, 'Simulation finished'
    

    
    
end program md_simulation