module output_module
    implicit none
    
contains
    subroutine output_positions(step, pos, n, filename)
        integer, intent(in) :: step, n
        real(8), dimension(n,3), intent(in) :: pos
        character(len=*), intent(in) :: filename
        integer :: i
        
        open(unit=10, file=filename, status='unknown', action='write', position='append')
        write(10, '(A, I6)') 'Step: ', step
        do i = 1, n
            write(10, '(3F10.5)') pos(i, 1), pos(i, 2), pos(i, 3)
        end do
        close(10)

    end subroutine output_positions
    
end module output_module