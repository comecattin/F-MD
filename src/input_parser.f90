module parser_module
    implicit none
contains
    subroutine read_config(input_file, n_atoms, n_steps, dt, box_length, stride)
        character(len=*), intent(out) :: input_file
        integer, intent(out) :: n_atoms, n_steps, stride
        real(8), intent(out) :: dt, box_length
        character(len=4096) :: parser, python, executable
        character(len=:), allocatable :: command
        integer :: status, command_status, slash, config_unit

        if (command_argument_count() /= 1) then
            print *, 'Usage: md_simulation input_file'
            stop 1
        end if
        call get_command_argument(1, input_file, status=status)
        if (status /= 0) error stop 'Input path is too long'
        call get_environment_variable('FMD_INPUT_PARSER', parser, status=status)
        if (status == -1) error stop 'Parser path is too long'
        if (len_trim(parser) == 0) then
            call get_command_argument(0, executable)
            slash = scan(trim(executable), '/', back=.true.)
            parser = executable(:slash) // 'input_parser.py'
        end if
        call get_environment_variable('FMD_PYTHON', python, status=status)
        if (status == -1) error stop 'Python path is too long'
        if (len_trim(python) == 0) python = 'python3'
        command = shell_quote(trim(python)) // ' ' // shell_quote(trim(parser)) // ' ' // &
            shell_quote(trim(input_file))
        call execute_command_line(command, exitstat=status, cmdstat=command_status)
        if (command_status /= 0) error stop 'Could not launch input parser'
        if (status /= 0) error stop 'Invalid simulation configuration'

        open(newunit=config_unit, file='parsed_config.tmp', status='old', action='read', iostat=status)
        if (status /= 0) error stop 'Input parser did not create a configuration'
        read(config_unit, *) n_atoms, n_steps, dt, box_length, stride
        close(config_unit, status='delete')
        print *, 'Input file: ', trim(input_file)
        print *, 'Number of atoms: ', n_atoms
        print *, 'Number of steps: ', n_steps
        print *, 'Time step (ps): ', dt
        print *, 'Box length (Angstrom): ', box_length
        print *, 'Stride: ', stride
    end subroutine read_config

    function shell_quote(value) result(quoted)
        character(len=*), intent(in) :: value
        character(len=:), allocatable :: quoted
        integer :: i
        quoted = achar(39)
        do i = 1, len(value)
            if (value(i:i) == achar(39)) then
                quoted = quoted // achar(39) // achar(34) // achar(39) // achar(34) // achar(39)
            else
                quoted = quoted // value(i:i)
            end if
        end do
        quoted = quoted // achar(39)
    end function shell_quote
end module parser_module
