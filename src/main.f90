program md_simulation
    use initialization_module, only : initialize
    use forces_module, only : compute_forces
    use integration_module, only: integrate
    use output_module, only: output_positions
    implicit none

    ! Define parameters
    integer :: n_atoms, n_steps, step, num_args
    real(8) :: dt, box_length
    character(len=100) :: input_file, output_file
    logical :: file_exists

    ! Define position, velocity and forces arrays
    real(8), allocatable:: positions(:,:), velocities(:,:), forces(:,:)

    num_args = command_argument_count()
    if (num_args /= 1) then
        print *, 'Usage: ./md_simulation input_file'
        stop
    else
        call get_command_argument(1, input_file)
        print *, 'Input file: ', input_file
    end if

    ! Read the input parameters
    open(unit=20, file='input.dat')
    read(20, *) n_atoms
    read(20, *) n_steps
    read(20, *) dt
    read(20, *) box_length
    close(20)

    ! Print the input parameters
    print *, 'Number of atoms: ', n_atoms
    print *, 'Number of steps: ', n_steps
    print *, 'Time step: ', dt
    print *, 'Box length: ', box_length

    ! Allocate the arrays
    allocate(positions(3, n_atoms))
    allocate(velocities(3, n_atoms))
    allocate(forces(3, n_atoms))


    ! Initialize the positions and velocities
    call initialize(positions, velocities, n_atoms, box_length)

    ! Open the output file
    output_file = 'trajectories.dat'
    ! Check if the file already exists
    inquire(file=output_file, exist=file_exists)
    if (file_exists) then
        print *, 'The file ', output_file, ' already exists. Removing it and the .arc file.'
        call system('rm ' //trim(output_file))
        call system('rm ' //trim(output_file(1:index(output_file, '.')-1)) // '.arc')
    end if

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

    ! Deallocate the arrays
    deallocate(positions, velocities, forces)
    

    
    
end program md_simulation