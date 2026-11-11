!> Result-file form of the turbulence output contract
!> (docs/drift-turbulence.md): one set of variable names for every provider.
!> Arrays are written in the normalised units of the model; the scalars
!> `unit_*` convert them to Gaussian cgs units.
module milca_turbulence_output
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_result_io, only: result_file_t
    use milca_turbulence_statistics, only: turbulence_spectrum_t, turbulence_scalar_t
    implicit none
    private
    public :: turbulence_units_t, write_turbulence_spectrum, write_turbulence_scalar

    !> Names of the normalised units and their size in Gaussian cgs units.
    type :: turbulence_units_t
        character(len=32) :: length = 'rho_s', time = '1/Omega_i', potential = 'T_e/e'
        real(dp) :: length_cm = 1, time_s = 1, potential_statvolt = 1
    end type

    character(len=*), parameter :: interval = 'candidate; statistical standard error'
contains
    subroutine write_turbulence_scalar(result, name, scalar, units, long_name)
        type(result_file_t), intent(inout) :: result
        character(len=*), intent(in) :: name, units, long_name
        type(turbulence_scalar_t), intent(in) :: scalar

        call result%write(name, scalar%value, units, long_name, status=interval)
        call result%write(name//'_error', scalar%error, units, &
            'one standard error of '//name//'; negative: not available', status=interval)
    end subroutine

    subroutine write_turbulence_spectrum(result, spectrum, units, model, provider)
        type(result_file_t), intent(inout) :: result
        type(turbulence_spectrum_t), intent(in) :: spectrum
        type(turbulence_units_t), intent(in) :: units
        character(len=*), intent(in) :: model, provider
        character(len=:), allocatable :: l, t, p, density

        l = trim(units%length); t = trim(units%time); p = trim(units%potential)
        density = '('//p//')^2 ('//l//')^2 '//t
        call result%attribute('turbulence_model', model)
        call result%attribute('turbulence_provider', provider)
        call result%attribute('turbulence_frame', &
            'zero mean ExB velocity; laboratory frequency is omega + ky V_E')
        call result%attribute('turbulence_time_convention', 'exp[i(kx x + ky y - omega t)]')
        call result%attribute('turbulence_parallel_model', 'one constant parallel wave number')
        call result%write('unit_length', units%length_cm, 'cm', 'size of '//l)
        call result%write('unit_time', units%time_s, 's', 'size of '//t)
        call result%write('unit_potential', units%potential_statvolt, 'statvolt', 'size of '//p)
        call result%write('kx', spectrum%kx, '1/('//l//')', 'radial wave number', dims='kx')
        call result%write('ky', spectrum%ky, '1/('//l//')', 'binormal wave number', dims='ky')
        call result%write('omega', spectrum%omega, '1/('//t//')', 'frequency', dims='omega')
        call result%write('S_phi_phi', spectrum%s_phi_phi, density, &
            'power spectral density of the potential; its sum times dkx dky domega '// &
            'is the variance', dims='kx,ky,omega', status='candidate; Welch estimate')
        call result%write('S_n_phi_re', spectrum%s_n_phi_re, density, &
            'real part of the cross spectral density <n conjg(phi)>, n in units of e phi/T_e', &
            dims='kx,ky,omega', status='candidate; Welch estimate')
        call result%write('S_n_phi_im', spectrum%s_n_phi_im, density, &
            'imaginary part of the cross spectral density <n conjg(phi)>', &
            dims='kx,ky,omega', status='candidate; Welch estimate')
        call result%write('lag_t', spectrum%lag_t, t, 'time lag', dims='lag_t')
        call result%write('lag_x', spectrum%lag_x, l, 'radial lag', dims='lag_x')
        call result%write('lag_y', spectrum%lag_y, l, 'binormal lag', dims='lag_y')
        call result%write('C_phi_t', spectrum%c_t, '('//p//')^2', &
            'autocorrelation of the potential at one point', dims='lag_t', status='candidate')
        call result%write('C_phi_x', spectrum%c_x, '('//p//')^2', &
            'radial correlation of the potential', dims='lag_x', status='candidate')
        call result%write('C_phi_y', spectrum%c_y, '('//p//')^2', &
            'binormal correlation of the potential', dims='lag_y', status='candidate')
        call write_turbulence_scalar(result, 'phi_variance', spectrum%variance, &
            '('//p//')^2', 'variance of the potential in the retained modes')
        call write_turbulence_scalar(result, 'tau_c_euler', spectrum%tau_euler, t, &
            'first 1/e crossing of the single-point autocorrelation')
        call write_turbulence_scalar(result, 'tau_c_mode', spectrum%tau_mode, t, &
            'first 1/e crossing of the summed mode-amplitude autocorrelation')
        call write_turbulence_scalar(result, 'ell_x', spectrum%ell_x, l, &
            'first 1/e crossing of the radial correlation')
        call write_turbulence_scalar(result, 'ell_y', spectrum%ell_y, l, &
            'first 1/e crossing of the binormal correlation')
        call write_turbulence_scalar(result, 'kubo', spectrum%kubo, '1', &
            'tau_c_euler times the rms radial ExB velocity over ell_x')
        call result%write('spectrum_segments', real(spectrum%segments, dp), '1', &
            'number of half-overlapping Hann segments')
        call result%write('spectrum_batches', real(spectrum%batches, dp), '1', &
            'number of complete batches behind the standard errors')
    end subroutine
end module milca_turbulence_output
