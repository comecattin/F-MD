program md_simulation
    use initialization_module, only : initialize
    use forces_module, only : compute_forces
    use integration_module, only: integrate
    use output_module, only: output_positions
    implicit none

    ! Define parameters
    integer, parameter :: n_atoms = 100
    integer, parameter :: n_steps = 1000
    real(8), parameter :: dt = 0.001
    real(8), parameter :: box_length = 10.0
    character(len=100) :: output_file

    ! Define position, velocity and forces arrays
    real(8), dimension(n_atoms, 3) :: positions, velocities, forces

    ! Other variables
    integer :: step

    ! Initialize the positions and velocities
    call initialize(positions, velocities, n_atoms, box_length)

    ! Open the output file
    output_file = 'trajectories.dat'

    ! Perform the molecular dynamics simulation
    do step = 1, n_steps
        print *, 'Step: ', step

        ! Compute the forces
        call compute_forces(positions, forces, n_atoms, box_length)

        ! Integrate positions and velocities
        call integrate(positions, velocities, forces, dt, n_atoms, box_length)
        
        ! Output the positions
        call output_positions(step, positions, n_atoms, output_file)

    end do

    print *, 'Simulation finished'
    

    
    
end program md_simulation