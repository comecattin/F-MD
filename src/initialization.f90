module initialization_module
    implicit none
    
contains
    subroutine initialize(pos, vel, n, box_length)

        real(8), dimension(n,3), intent(out) :: pos, vel
        integer, intent(in) :: n
        real(8), intent(in) :: box_length
        integer :: i, j

        call random_seed() ! Initialize the random number generator
        call random_number(pos) ! Initialize the positions between 0 and 1
        pos = pos * box_length

        call random_number(vel)
        vel = vel - 0.5 ! Shift the velocities to be between -0.5 and 0.5

    end subroutine initialize

    subroutine initialize_water(pos, vel, charges, n, box_length)
        real(8), dimension(n,3), intent(out) :: pos, vel
        real(8), dimension(n), intent(out) :: charges
        real(8) :: oh_distance, hoh_angle, o_charge, h_charge
        real(8) :: pi
        real(8) :: x_oxygen, y_oxygen, z_oxygen
        real(8) :: x_h1, y_h1, z_h1, x_h2, y_h2, z_h2
        integer, intent(in) :: n
        real(8), intent(in) :: box_length
        integer :: i, j, overlap_i
        logical :: overlap
        real(8) :: cos_theta, theta
        real(8), dimension(3):: vec1, vec2

        if (mod(n, 3) /= 0) then
            print *, "Error: The number of atoms must be a multiple of 3"
            stop
        end if

        call random_seed()

        ! Velocities
        call random_number(vel)
        vel = vel - 0.5

        pi = 3.141592653589793d0

        ! TIP3P water model
        oh_distance = 1.012                 ! Angstroms
        hoh_angle = 113.24 * pi / 180.0     ! Radians
        o_charge = -0.834d0
        h_charge = 0.417d0



            
            ! Hydrogen atoms relative to the oxygen atom
            pos(i+1 ,1) = pos(i, 1) + oh_distance
            pos(i+1, 2) = pos(i, 2)
            pos(i+1, 3) = pos(i, 3)

            pos(i+2 ,1) = pos(i, 1) + oh_distance * cos(hoh_angle)
            pos(i+2, 2) = pos(i, 2) + oh_distance * sin(hoh_angle)
            pos(i+2, 3) = pos(i, 3)

            ! Periodic boundary conditions
            pos(i+1, :) = mod(pos(i+1, :), box_length)
            pos(i+2, :) = mod(pos(i+2, :), box_length)


            ! Assign charges to the atoms
            charges(i) = o_charge
            charges(i+1) = h_charge
            charges(i+2) = h_charge

        end do
    end subroutine initialize_water
    
end module initialization_module