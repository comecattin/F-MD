! Physical units: Angstrom, picosecond, g/mol, elementary charge, kcal/mol.
module units_module
    implicit none
    private
    public :: force_to_acceleration, kinetic_to_kcal, coulomb_constant

    ! F [(kcal/mol)/A] / M [g/mol] -> acceleration [A/ps**2].
    ! 4184 J/kcal / (1e-3 kg/g * (1e-10 m/A)**2 / (1e-12 s/ps)**2).
    real(8), parameter :: force_to_acceleration = 418.4d0
    ! M [g/mol] * v**2 [(A/ps)**2] -> energy [kcal/mol].
    real(8), parameter :: kinetic_to_kcal = 1.0d0 / force_to_acceleration

    ! SI constants from CODATA 2022: https://physics.nist.gov/cuu/pdf/all.pdf
    real(8), parameter :: avogadro = 6.02214076d23
    real(8), parameter :: elementary_charge = 1.602176634d-19
    real(8), parameter :: vacuum_permittivity = 8.8541878188d-12
    ! N_A * e**2 / (4*pi*epsilon_0), converted to kcal*A/(mol*e**2).
    real(8), parameter :: coulomb_constant = avogadro * elementary_charge**2 / &
        (4.0d0 * acos(-1.0d0) * vacuum_permittivity * 1.0d-10 * 4184.0d0)
end module units_module
