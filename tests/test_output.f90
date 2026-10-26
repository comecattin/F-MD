program test_output
    use output_module, only: output_energies
    implicit none
    character(len=100) :: label
    real(8) :: energy
    integer :: step

    call output_energies(100000, 123456.789d0, -987654.321d0, -864197.532d0, 'energies.dat')
    open(unit=20, file='energies.dat', status='old', action='read')
    read(20, *) label, step
    if (step /= 100000) error stop 'Output step is incorrect'
    read(20, *) label, label, energy
    if (abs(energy - 123456.789d0) > 1.0d-8) error stop 'Kinetic energy output is incorrect'
    read(20, *) label, label, energy
    if (abs(energy + 987654.321d0) > 1.0d-8) error stop 'Potential energy output is incorrect'
    read(20, *) label, label, energy
    if (abs(energy + 864197.532d0) > 1.0d-8) error stop 'Total energy output is incorrect'
    close(20)
    print *, 'Energy output checks passed'
end program test_output
