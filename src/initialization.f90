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
        integer :: i

        if (mod(n, 3) /= 0) then
            print *, "Error: The number of atoms must be a multiple of 3"
            stop
        end if

        call random_seed()

        pi = 3.141592653589793d0

        ! TIP3P water model
        oh_distance = 0.9572d0              ! Angstroms
        hoh_angle = 104.52d0 * pi / 180.0d0 ! Radians
        o_charge = -0.834d0
        h_charge = 0.417d0


        ! Initialize the positions at random for every atom
        call random_number(pos)
        pos = pos * box_length
        call random_number(vel)
        vel = vel - 0.5

        ! Move hydrogen atoms and assign charges
        do i = 1, n, 3
            ! Oxygen atom
            x_oxygen = pos(i,1)
            y_oxygen = pos(i,2)
            z_oxygen = pos(i,3)

            ! Hydrogen atoms relative to the oxygen atom
            x_h1 = x_oxygen + oh_distance * cos(hoh_angle / 2.0)
            y_h1 = y_oxygen + oh_distance * sin(hoh_angle / 2.0)
            z_h1 = z_oxygen
            x_h2 = x_oxygen - oh_distance * cos(hoh_angle / 2.0)
            y_h2 = y_oxygen + oh_distance * sin(hoh_angle / 2.0)
            z_h2 = z_oxygen
            
            ! Perioidic boundary conditions
            x_h1 = mod(x_h1, box_length)
            y_h1 = mod(y_h1, box_length)
            z_h1 = mod(z_h1, box_length)

            x_h2 = mod(x_h2, box_length)
            y_h2 = mod(y_h2, box_length)
            z_h2 = mod(z_h2, box_length)

            ! Assign positions to the hydrogen atoms
            pos(i+1,1) = x_h1
            pos(i+1,2) = y_h1
            pos(i+1,3) = z_h1

            pos(i+2,1) = x_h2
            pos(i+2,2) = y_h2
            pos(i+2,3) = z_h2

            ! Assign charges to the atoms
            charges(i) = o_charge
            charges(i+1) = h_charge
            charges(i+2) = h_charge
        
        end do
    end subroutine initialize_water
    
end module initialization_module