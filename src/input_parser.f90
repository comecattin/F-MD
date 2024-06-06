module parser_module
    implicit none
    
contains

    subroutine read_config(input_file, n_atoms, n_steps, dt, box_length, tolerance, max_iter)
        !
        ! Routine to read the input configuration file
        !
        character(len=100):: input_file
        
        integer :: num_args
        character(len=100) :: command

        integer, intent(out) :: n_atoms, n_steps, max_iter
        real(8), intent(out) :: dt, box_length, tolerance

        ! Check that CLI parameters is correct
        num_args = command_argument_count()
        if (num_args /= 1) then
            print *, 'Usage: ./md_simulation input_file'
            stop
        else
            call get_command_argument(1, input_file)
            print *, 'Input file: ', input_file
        end if

        ! Create the python command
        write(command, '(A, A)') 'python input_parser.py ', trim(input_file)
        call system(command)

        ! Read the temporary file created by python 
        input_file = 'parsed_config.tmp'
        open(unit=20, file=input_file)
        read(20, *) n_atoms
        read(20, *) n_steps
        read(20, *) dt
        read(20, *) box_length
        read(20, *) tolerance
        read(20, *) max_iter
        close(20)

        print *, 'Number of atoms: ', n_atoms
        print *, 'Number of steps: ', n_steps
        print *, 'Time step: ', dt
        print *, 'Box length: ', box_length
        print *, 'SHAKE Tolerance: ', tolerance
        print *, 'SHAKE Max iterations: ', max_iter

        ! Remove the temporary file
        call system('rm parsed_config.tmp')


    end subroutine read_config
end module parser_module