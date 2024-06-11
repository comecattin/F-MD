program md_simulation
    use parser_module, only : read_config
    use initialization_module, only : initialize_water
    use forces_module, only : compute_forces
    use energies_module
    use integration_module, only: integrate
    use output_module, only: output_positions, output_energies, output_positions_water
    implicit none

    ! Define parameters
    integer :: n_atoms, n_steps, step, num_args, stride
    real(8) :: dt, box_length
    character(len=100) :: input_file, output_file, output_file_energies, command
    logical :: file_exists

    ! Define position, velocity and forces arrays
    real(8), allocatable:: positions(:,:), velocities(:,:), forces(:,:), charges(:), mass(:)

    ! Kinetic, potential and total energies
    real(8) :: ke, pe, te

    call read_config(input_file, n_atoms, n_steps, dt, box_length, stride)

    ! Allocate the arrays
    allocate(positions(n_atoms, 3))
    allocate(velocities(n_atoms, 3))
    allocate(forces(n_atoms, 3))
    allocate(charges(n_atoms))
    allocate(mass(n_atoms))


    ! Initialize the positions and velocities
    call initialize_water(positions, velocities, charges, n_atoms, box_length, mass)

    ! Open the output file
    output_file = 'trajectories.dat'
    output_file_energies = 'energies.dat'
    ! Check if the file already exists
    inquire(file=output_file, exist=file_exists)
    if (file_exists) then
        print *, 'The file ', output_file, ' already exists. Removing it and the .arc file.'
        call system('rm ' //trim(output_file))
        call system('rm ' //trim(output_file(1:index(output_file, '.')-1)) // '.arc')
    end if
    inquire(file=output_file_energies, exist=file_exists)
    if (file_exists) then
        print *, 'The file ', output_file_energies, ' already exists. Removing it'
        call system('rm ' //trim(output_file_energies))
    end if

    call output_positions_water(0, positions, n_atoms, output_file)

    ! Perform the molecular dynamics simulation
    do step = 1, n_steps

        ! Compute the forces
        call compute_forces(positions, charges, forces, n_atoms, box_length)

        ! Integrate positions and velocities
        call integrate(positions, velocities, forces, dt, n_atoms, box_length, mass, charges)

        ! Compute the energies
        call compute_energies(positions, velocities, mass, charges, n_atoms, box_length, ke, pe, te)
        
        if (mod(step, stride) .eq. 0) then
            print *, ''
            print *, 'Step: ', step
            ! Output the positions and the energies
            call output_positions_water(step, positions, n_atoms, output_file)
            call output_energies(step, ke, pe, te, output_file_energies)
        end if

    end do

    print *, 'Simulation finished'

    ! Deallocate the arrays
    deallocate(positions, velocities, forces, charges, mass)
    

    
    
end program md_simulation