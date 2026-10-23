module forces_module
    implicit none
    real(8), parameter :: pi = 3.14159265358979323846
    real(8), parameter :: epsilon = 0.1554253          ! kcal/mol
    real(8), parameter :: sigma = 3.165492             ! Angstrom 
    real(8), parameter :: ka_HOH = 75.90               ! kcal.mol-1.rad-2
    real(8), parameter :: kb_OH = 1059.162             ! kcal.mol-1.angstrom-2
    real(8), parameter :: alpha = 1.0
    real(8), parameter :: r_eq_OH = 1.012              ! Angstrom
    real(8), parameter :: angle_eq = 113.24            ! deg
    
contains

    subroutine compute_forces(position, charges, forces, n, box_length)
        
        real(8), dimension(n, 3), intent(in) :: position
        real(8), dimension(n), intent(in) :: charges
        real(8), dimension(n, 3), intent(out) :: forces
        integer, intent(in) :: n
        real(8), intent(in) :: box_length
        integer :: i, j
        real(8), dimension(3) :: dist, vec1, vec2, angle_force1, angle_force2
        real(8) :: r, r2, f_yukawa, force_lj, force_harmonic_bond, force_harmonic_angle
        integer :: bond_type
        real(8) :: cos_theta, theta, length1, length2, sin_theta
        real(8) :: angle_eq_rad
        logical :: is_OO, is_OH, is_HH

        
        angle_eq_rad = angle_eq * pi / 180.0 ! Rad
        

        forces = 0.0

        ! Intra molecular interactions
        do i = 1, n-1, 3
            ! Harmonic potential on HOH angle
            vec1 = position(i+1, :) - position(i, :)
            vec2 = position(i+2, :) - position(i, :)

            ! Minimum image convention
            vec1 = vec1 - nint(vec1 / box_length) * box_length
            vec2 = vec2 - nint(vec2 / box_length) * box_length
            
            length1 = norm2(vec1)
            length2 = norm2(vec2)
            if (min(length1, length2) <= tiny(1.0d0)) error stop 'Zero-length O-H bond'
            cos_theta = dot_product(vec1, vec2) / (length1 * length2)
            cos_theta = max(-1.0d0, min(1.0d0, cos_theta))
            theta = acos(cos_theta)
            sin_theta = sqrt(max(0.0d0, 1.0d0 - cos_theta**2))
            if (sin_theta <= 1.0d-12) error stop 'Undefined force for collinear H-O-H angle'
            
            call harmonic(theta, angle_eq_rad, ka_HOH, force_harmonic_angle)

            ! F = -dV/dtheta * dtheta/dr for V = ka/2 * (theta-theta_eq)**2.
            angle_force1 = force_harmonic_angle / (length1 * sin_theta) * &
                (vec2 / length2 - cos_theta * vec1 / length1)
            angle_force2 = force_harmonic_angle / (length2 * sin_theta) * &
                (vec1 / length1 - cos_theta * vec2 / length2)
            forces(i, :) = forces(i, :) - angle_force1 - angle_force2
            forces(i + 1, :) = forces(i + 1, :) + angle_force1
            forces(i + 2, :) = forces(i + 2, :) + angle_force2


            ! Harmonic potential on O-H bond
            !    O-H1 bond
            call harmonic(norm2(vec1), r_eq_OH, kb_OH, force_harmonic_bond)
            forces(i, :) = forces(i, :) + force_harmonic_bond * vec1/norm2(vec1)
            forces(i + 1, :) = forces(i + 1, :) - force_harmonic_bond * vec1/norm2(vec1)
            
            !    O-H2 bond
            call harmonic(norm2(vec2), r_eq_OH, kb_OH, force_harmonic_bond)
            forces(i, :) = forces(i, :) + force_harmonic_bond * vec2/norm2(vec2)
            forces(i + 2, :) = forces(i + 2, :) - force_harmonic_bond * vec2/norm2(vec2)
            
        end do

        ! Van der Waals interactions
        do i = 1, n-4, 3
            do j = i+3, n, 3
                dist = position(i, :) - position(j, :)
                ! Minimum image convention
                dist = dist - nint(dist / box_length) * box_length
                r = norm2(dist)
                call lennard_jones(r, force_lj)
                forces(i, :) = forces(i, :) + force_lj * dist/r
                forces(j, :) = forces(j, :) - force_lj * dist/r
            end do
        end do

        ! Coulomb interactions
        do i = 1, n-1
            do j = i+1, n
                dist = position(i, :) - position(j, :)
                ! Minimum image convention
                dist = dist - nint(dist / box_length) * box_length
                r = norm2(dist)
                call coulomb_yukawa(r, 1.0_8, charges(i), charges(j), alpha, f_yukawa)
                forces(i, :) = forces(i, :) + f_yukawa * dist/r
                forces(j, :) = forces(j, :) - f_yukawa * dist/r
            end do
        end do

    end subroutine compute_forces

    subroutine coulomb_yukawa(r, mass, q1, q2, alpha_param, f_yukawa)

        real(8), intent(in) :: r
        real(8), intent(in) :: q1, q2
        real(8), intent(in) :: mass, alpha_param
        real(8), intent(out) :: f_yukawa

        f_yukawa = q1*q2 * exp(-alpha_param * mass * r) * (alpha_param * mass * r + 1)/r**2

    end subroutine coulomb_yukawa

    subroutine lennard_jones(r, force_lj)
    
        real(8), intent(in) :: r
        real(8), intent(out) :: force_lj

        force_lj = 48*epsilon/sigma * ((sigma / r)**13 - 0.5 * (sigma / r)**7)

    end subroutine lennard_jones

    subroutine harmonic(r, r_eq, kb, force_harmonic)
    
        real(8), intent(in) :: r
        real(8), intent(in) :: r_eq, kb
        real(8), intent(out) :: force_harmonic

        force_harmonic = kb * (r - r_eq)

    end subroutine harmonic

    subroutine harmonic_angle(dea, n, position, ideal, force)
        real(8), dimension(n, 3), intent(out) :: dea
        integer, intent(in) :: n
        real(8), dimension(n, 3), intent(in) :: position
        real(8), intent(in) :: ideal
        real(8), intent(in) :: force
        real(8) :: xia, yia, zia, xib, yib, zib, xic, yic, zic
        real(8) :: xab, yab, zab, xcb, ycb, zcb
        real(8) :: rab2, rcb2
        real(8) :: xp, yp, zp, rp
        real(8) :: dot, cosine, angle
        integer :: i
        real(8) :: dt, dt2, deddt
        real(8) :: terma, termc, dedxia, dedyia, dedzia, dedxic, dedyic, dedzic, dedxib, dedyib, dedzib

        do i = 1, n
            dea(i, 1) = 0.0d0
            dea(i, 2) = 0.0d0
            dea(i, 3) = 0.0d0
        end do

        do i = 1, n-1, 3

            xia = position(i+1, 1)
            yia = position(i+1, 2)
            zia = position(i+1, 3)
            xib = position(i, 1)
            yib = position(i, 2)
            zib = position(i, 3)
            xic = position(i+2, 1)
            yic = position(i+2, 2)
            zic = position(i+2, 3)

            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib

            rab2 = xab*xab + yab*yab + zab*zab
            rcb2 = xcb*xcb + ycb*ycb + zcb*zcb

            xp = ycb*zab - zcb*yab
            yp = zcb*xab - xcb*zab
            zp = xcb*yab - ycb*xab
            rp = sqrt(xp*xp + yp*yp + zp*zp)

            dot = xab*xcb + yab*ycb + zab*zcb
            cosine = dot / sqrt(rab2 * rcb2)
            cosine = min(1.0d0, max(-1.0d0, cosine))
            angle = acos(cosine)

            dt = angle - ideal
            dt2 = dt * dt
            deddt = force * dt * 2.0d0

            terma = -deddt / (rab2 * rp)
            termc =  deddt / (rcb2 * rp)
            dedxia = terma * (yab*zp-zab*yp)
            dedyia = terma * (zab*xp-xab*zp)
            dedzia = terma * (xab*yp-yab*xp)
            dedxic = termc * (ycb*zp-zcb*yp)
            dedyic = termc * (zcb*xp-xcb*zp)
            dedzic = termc * (xcb*yp-ycb*xp)
            dedxib = -dedxia - dedxic
            dedyib = -dedyia - dedyic
            dedzib = -dedzia - dedzic

            dea(i+1, 1) = dea(i+1, 1) + dedxia
            dea(i+1, 2) = dea(i+1, 2) + dedyia
            dea(i+1, 3) = dea(i+1, 3) + dedzia
            dea(i, 1) = dea(i, 1) + dedxib
            dea(i, 2) = dea(i, 2) + dedyib
            dea(i, 3) = dea(i, 3) + dedzib
            dea(i+2, 1) = dea(i+2, 1) + dedxic
            dea(i+2, 2) = dea(i+2, 2) + dedyic
            dea(i+2, 3) = dea(i+2, 3) + dedzic

        end do

    end subroutine harmonic_angle




end module forces_module
