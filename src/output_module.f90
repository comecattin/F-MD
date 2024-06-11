module output_module
    implicit none
    
contains
    subroutine output_positions(step, pos, n, filename)
        integer, intent(in) :: step, n
        real(8), dimension(n,3), intent(in) :: pos
        character(len=*), intent(in) :: filename
        integer :: i
        character (len=100) :: arc_file

        arc_file = filename(1:index(filename, '.')-1) // '.arc'
        
        open(unit=10, file=filename, status='unknown', action='write', position='append')
        open(unit=11, file=arc_file, status='unknown', action='write', position='append')
        write(10, '(A, I6)') 'Step: ', step
        write(11, '(A, I0)') '', n
        do i = 1, n
            write(10, '(3F10.5)') pos(i, 1), pos(i, 2), pos(i, 3)
            write(11, *) i, 'H', pos(i, 1), pos(i, 2), pos(i, 3)
        end do
        close(10)
        close(11)

    end subroutine output_positions

    subroutine output_positions_water(step, pos, n, filename)
        integer, intent(in) :: step, n
        real(8), dimension(n,3), intent(in) :: pos
        character(len=*), intent(in) :: filename
        integer :: i
        character (len=100) :: arc_file

        arc_file = filename(1:index(filename, '.')-1) // '.arc'
        
        open(unit=10, file=filename, status='unknown', action='write', position='append')
        open(unit=11, file=arc_file, status='unknown', action='write', position='append')
        write(10, '(A, I6)') 'Step: ', step
        write(11, '(A, I0)') '', n
        do i = 1, n, 3
            write(10, '(3F10.5)') pos(i, 1), pos(i, 2), pos(i, 3)
            write(10, '(3F10.5)') pos(i+1, 1), pos(i+1, 2), pos(i+1, 3)
            write(10, '(3F10.5)') pos(i+2, 1), pos(i+2, 2), pos(i+2, 3)
            write(11, *) i, 'O', pos(i, 1), pos(i, 2), pos(i, 3)
            write(11, *) i+1, 'H', pos(i+1, 1), pos(i+1, 2), pos(i+1, 3)
            write(11, *) i+2, 'H', pos(i+2, 1), pos(i+2, 2), pos(i+2, 3)
        end do
        close(10)
        close(11)

    end subroutine output_positions_water

    subroutine output_energies(step, ek, ep, et, filename)
        integer, intent(in) :: step
        real(8), intent(in) :: ek, ep, et
        character(len=*), intent(in) :: filename

        open(unit=10, file=filename, status='unknown', action='write', position='append')
        write(10, '(A, I6)') 'Step: ', step
        write(10, '(A, F10.5)') 'Kinetic Energy: ', ek
        write(10, '(A, F10.5)') 'Potential Energy: ', ep
        write(10, '(A, F10.5)') 'Total Energy: ', et
        close(10)

        print *, "Kinetic energy: ", ek
        print *, "Potential energy: ", ep
        print *, "Total energy: ", et

    end subroutine output_energies
    
end module output_module