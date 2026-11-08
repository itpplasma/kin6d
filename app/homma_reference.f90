program homma_reference_main
    ! Benchmark outputs B01-B03 from the Rosenbluth oracle: Maxwellian A, D (B01),
    ! distortion A^(T), D^(T) (B02), and thermal force F_T (B03).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use milca_rosenbluth_reference, only: homma_pair_t, homma_coeff_t, homma_vt, &
        homma_point, homma_force_maxwellian, homma_force_shell, homma_force_beam, &
        homma_closure
    use milca_result_io, only: run_t, result_file_t, runtime_git
    implicit none
    integer, parameter :: n_b01 = 7, n_b02 = 5, n_grad2 = 6, n_shell = 4, n_beam = 2
    integer, parameter :: n_charge = 4, n_mass = 3, n_dens = 3, n_temp = 3, n_grad3 = 4
    integer, parameter :: n_dir = 4, n_dist = 1 + n_shell + n_shell*n_beam
    integer, parameter :: n_b03 = n_dist*n_charge*n_mass*n_dens*n_temp*n_grad3
    real(dp), parameter :: e_charge = 1.602176634e-19_dp, amu = 1.66053906660e-27_dp
    real(dp), parameter :: unit_dir = 0.57735026918962576_dp
    real(dp) :: ppu = 8.0_dp, xmax = 7.0_dp
    real(dp) :: speeds_b01(n_b01) = [0.0_dp, 0.1_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp, 10.0_dp]
    real(dp) :: speeds_b02(n_b02) = [0.1_dp, 0.5_dp, 1.0_dp, 2.0_dp, 5.0_dp]
    real(dp) :: dirs(3, n_dir) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, unit_dir, unit_dir, unit_dir], [3, n_dir])
    real(dp) :: grad_b02(3, n_grad2) = reshape([0.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, -5.0_dp, &
        100.0_dp, 0.0_dp, 0.0_dp, 300.0_dp, 0.0_dp, 0.0_dp, -100.0_dp, 0.0_dp, 0.0_dp, &
        300.0_dp, 0.0_dp, 5.0_dp], [3, n_grad2])
    real(dp) :: shell_x(n_shell) = [0.05_dp, 0.1_dp, 0.5_dp, 1.0_dp]
    real(dp) :: beam_dirs(3, n_beam) = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], &
        [3, n_beam])
    real(dp) :: grad_b03(3, n_grad3) = reshape([0.0_dp, 0.0_dp, 5.0_dp, 100.0_dp, 0.0_dp, 0.0_dp, &
        300.0_dp, 0.0_dp, 0.0_dp, 300.0_dp, 0.0_dp, 5.0_dp], [3, n_grad3])
    real(dp) :: charges(n_charge) = [1.0_dp, 3.0_dp, 10.0_dp, 30.0_dp]
    real(dp) :: masses_amu(n_mass) = [183.84_dp, 12.011_dp, 39.948_dp]
    real(dp) :: densities(n_dens) = [1.0e19_dp, 1.0e20_dp, 1.0e21_dp]
    real(dp) :: temps_ev(n_temp) = [10.0_dp, 50.0_dp, 200.0_dp]
    character(len=256) :: case_files(3) = 'not supplied', case_sha256(3) = 'not supplied'
    namelist /homma_reference/ ppu, xmax, speeds_b01, speeds_b02, dirs, &
        grad_b02, shell_x, beam_dirs, grad_b03, charges, masses_amu, densities, temps_ev, &
        case_files, case_sha256

    type(run_t) :: run
    type(result_file_t) :: result
    type(homma_pair_t) :: p
    type(homma_coeff_t) :: co
    character(len=256) :: msg
    integer :: unit, parse_status, i, j, k, ig, ic, im, id, it, dist, row
    integer :: n_b1, n_b2
    real(dp) :: vt, v3(3), speed, b_z(3), grad_si(3), force(3), ferr
    real(dp), allocatable :: a_b01(:, :), d_b01(:, :, :), e_b01(:, :), s_b01(:), dir_b01(:, :)
    real(dp), allocatable :: a_b02(:, :), d_b02(:, :, :), e_b02(:, :), s_b02(:), dir_b02(:, :)
    real(dp), allocatable :: f_b03(:, :), c_b03(:), m_b03(:), n_b03_v(:), t_b03(:), g_b02(:, :)
    real(dp), allocatable :: e_b03(:)
    character(len=64) :: label, param
    character(len=:), allocatable :: commit, tree

    call run%start('homma_reference')
    parse_status = 0; msg = ''
    call run%open_config(unit)
    if (run%has_config()) read (unit, nml=homma_reference, iostat=parse_status, iomsg=msg)
    call run%close_config(unit, 'homma_reference', parse_status, msg)
    if (ppu <= 0.0_dp .or. xmax <= 0.0_dp) error stop 'ppu and xmax must be positive'

    call run%open_resolved(unit)
    write (unit, nml=homma_reference)
    call run%close_resolved(unit)

    b_z = [0.0_dp, 0.0_dp, 1.0_dp]
    ! Background H+ of B01-B03: n=1e20 m^-3, T=50 eV, lnL=15, B along z; W3+ test.
    p%n_b = 1.0e20_dp; p%t_b = 50.0_dp*e_charge; p%m_b = 1.67262192369e-27_dp
    p%q_b = e_charge; p%lnlambda = 15.0_dp
    p%m_a = 183.84_dp*amu; p%q_a = 3.0_dp*e_charge
    vt = homma_vt(p)
    call execute_command_line('mkdir -p output/homma')
    call runtime_git(commit, tree, 'homma_reference')

    ! B01: Maxwellian coefficients for each listed speed and direction.
    p%heat_flux = 0.0_dp
    n_b1 = n_b01*n_dir
    allocate (a_b01(3, n_b1), d_b01(3, 3, n_b1), e_b01(2, n_b1), s_b01(n_b1), dir_b01(3, n_b1))
    row = 0
    do i = 1, n_b01
        do j = 1, n_dir
            row = row + 1
            speed = speeds_b01(i)
            v3 = speed*vt*dirs(:, j)
            call homma_point(p, v3, ppu, xmax, co)
            a_b01(:, row) = co%a_m
            d_b01(:, :, row) = co%d_m
            e_b01(:, row) = [co%err_a_m, co%err_d_m]
            s_b01(row) = speeds_b01(i)
            dir_b01(:, row) = dirs(:, j)
        end do
    end do

    ! B02: distortion corrections A^(T), D^(T) for each heat flux from grad T.
    n_b2 = n_grad2*n_b02*n_dir
    allocate (a_b02(3, n_b2), d_b02(3, 3, n_b2), e_b02(2, n_b2), s_b02(n_b2), dir_b02(3, n_b2))
    allocate (g_b02(3, n_b2))
    row = 0
    do ig = 1, n_grad2
        grad_si = grad_b02(:, ig)*e_charge
        call set_heat_flux(p, b_z, grad_si)
        do i = 1, n_b02
            do j = 1, n_dir
                row = row + 1
                v3 = speeds_b02(i)*vt*dirs(:, j)
                call homma_point(p, v3, ppu, xmax, co)
                a_b02(:, row) = co%a_t
                d_b02(:, :, row) = co%d_t
                e_b02(:, row) = [co%err_a_t, co%err_d_t]
                s_b02(row) = speeds_b02(i)
                dir_b02(:, row) = dirs(:, j)
                g_b02(:, row) = grad_b02(:, ig)
            end do
        end do
    end do

    ! B03: thermal force per test particle for each distribution and parameter set.
    allocate (f_b03(3, n_b03), c_b03(n_b03), m_b03(n_b03), n_b03_v(n_b03), t_b03(n_b03))
    allocate (e_b03(n_b03))
    open (newunit=unit, file='output/homma/B03.csv', status='replace', action='write')
    call write_header(unit, 3, 'F_T = m_a <A^(T)> per test particle; test temperature '// &
        'equals the background temperature; beams use the shell speeds', &
        'FT N; error_estimate N (largest component change under refinement)')
    write (unit, '(a)') 'distribution,parameter,charge_number,test_mass_amu,density_m3,'// &
        'temperature_ev,grad_x_ev_m,grad_y_ev_m,grad_z_ev_m,FT_x,FT_y,FT_z,error_estimate,status'
    row = 0
    do ig = 1, n_grad3
        do it = 1, n_temp
            do id = 1, n_dens
                do im = 1, n_mass
                    do ic = 1, n_charge
                        do dist = 1, n_dist
                            row = row + 1
                            p%n_b = densities(id)
                            p%t_b = temps_ev(it)*e_charge
                            p%m_a = masses_amu(im)*amu
                            p%q_a = charges(ic)*e_charge
                            call set_heat_flux(p, b_z, grad_b03(:, ig)*e_charge)
                            call thermal_force(p, dist, force, ferr, label, param)
                            f_b03(:, row) = force
                            e_b03(row) = ferr
                            c_b03(row) = charges(ic)
                            m_b03(row) = masses_amu(im)
                            n_b03_v(row) = densities(id)
                            t_b03(row) = temps_ev(it)
                            write (unit, '(a,",",a,11(",",es24.16e3),",candidate")') &
                                trim(label), trim(param), charges(ic), masses_amu(im), &
                                densities(id), temps_ev(it), grad_b03(:, ig), force, ferr
                        end do
                    end do
                end do
            end do
        end do
    end do
    close (unit)

    call write_b01_csv()
    call write_b02_csv()

    call execute_command_line('mkdir -p output')
    call result%create(run)
    call result%write('b01_speed_norm', s_b01, '1', 'normalised speed |v|/v_t (B01)', &
        dims='b01_row')
    call result%write('b01_A', a_b01, 'm s^-2', 'Maxwellian drift A (B01)', &
        dims='component,b01_row')
    call result%write('b01_D', d_b01, 'm^2 s^-3', 'Maxwellian diffusion D (B01)', &
        dims='component,component2,b01_row')
    call result%write('b01_err_A', e_b01(1, :), 'm s^-2', 'error estimate of A (B01)', &
        dims='b01_row')
    call result%write('b01_err_D', e_b01(2, :), 'm^2 s^-3', 'error estimate of D (B01)', &
        dims='b01_row')
    call result%write('b02_speed_norm', s_b02, '1', 'normalised speed (B02)', dims='b02_row')
    call result%write('b02_A_T', a_b02, 'm s^-2', 'distortion drift A^(T) (B02)', &
        dims='component,b02_row')
    call result%write('b02_D_T', d_b02, 'm^2 s^-3', 'distortion diffusion D^(T) (B02)', &
        dims='component,component2,b02_row')
    call result%write('b02_err_A_T', e_b02(1, :), 'm s^-2', 'error estimate of A^(T) (B02)', &
        dims='b02_row')
    call result%write('b02_err_D_T', e_b02(2, :), 'm^2 s^-3', 'error estimate of D^(T) (B02)', &
        dims='b02_row')
    call result%write('b03_F_T', f_b03, 'N', 'thermal force per test particle (B03)', &
        dims='component,b03_row')
    call result%write('b03_err', e_b03, 'N', 'error estimate of F_T (B03)', dims='b03_row')
    call result%write('b03_charge_number', c_b03, '1', 'test charge number (B03)', &
        dims='b03_row')
    call result%write('b03_mass_amu', m_b03, 'amu', 'test mass (B03)', dims='b03_row')
    call result%write('b03_density_m3', n_b03_v, 'm^-3', 'background density (B03)', &
        dims='b03_row')
    call result%write('b03_temperature_ev', t_b03, 'eV', 'background temperature (B03)', &
        dims='b03_row')
    call result%close()
    print '(a)', 'homma_reference wrote output/homma_reference.nc and output/homma/B0{1,2,3}.csv'

contains

    !> Heat flux of the background from grad T (eV/m given as SI J/m) with B along b.
    subroutine set_heat_flux(pair, b, grad_joule_m)
        type(homma_pair_t), intent(inout) :: pair
        real(dp), intent(in) :: b(3), grad_joule_m(3)
        real(dp) :: q(3), bmag
        bmag = sqrt(dot_product(b, b))
        call homma_closure(pair, b/bmag, bmag, grad_joule_m, q)
        pair%heat_flux = q
    end subroutine set_heat_flux

    !> Thermal force of distribution dist: 1 Maxwellian, then shells, then beams
    !> (each shell speed along each beam direction).
    subroutine thermal_force(pair, dist, force_out, err_out, label_out, param_out)
        type(homma_pair_t), intent(in) :: pair
        integer, intent(in) :: dist
        real(dp), intent(out) :: force_out(3), err_out
        character(len=64), intent(out) :: label_out, param_out
        integer :: ish, ib, k
        if (dist == 1) then
            label_out = 'maxwellian'
            param_out = '-'
            call homma_force_maxwellian(pair, pair%t_b, ppu, xmax, force_out, err_out)
        else if (dist <= 1 + n_shell) then
            ish = dist - 1
            label_out = 'shell'
            write (param_out, '(f0.4)') shell_x(ish)
            call homma_force_shell(pair, shell_x(ish), ppu, xmax, force_out, err_out)
        else
            k = dist - 2 - n_shell
            ish = 1 + mod(k, n_shell)
            ib = 1 + k/n_shell
            label_out = 'beam'
            write (param_out, '(a,f0.4,a,f0.4,",",f0.4,",",f0.4,a)') '"speed=', shell_x(ish), &
                ';dir=', beam_dirs(:, ib), '"'
            call homma_force_beam(pair, shell_x(ish), beam_dirs(:, ib), ppu, xmax, force_out, &
                err_out)
        end if
    end subroutine thermal_force

    !> Header of a result file (thermal-force-benchmark conventions, Result files).
    subroutine write_header(u, case_index, deviation, units)
        integer, intent(in) :: u, case_index
        character(len=*), intent(in) :: deviation, units
        write (u, '(a)') '# code,milca homma_reference (numerical Rosenbluth integrals)'
        write (u, '(a)') '# code_revision,'//commit//' ('//tree//')'
        write (u, '(a)') '# case_file,'//trim(case_files(case_index))
        write (u, '(a)') '# case_file_sha256,'//trim(case_sha256(case_index))
        write (u, '(a)') '# convention_deviations,written to output/homma/ (not results/<code>/); '// &
            deviation
        write (u, '(a)') '# precision,double (IEEE binary64)'
        write (u, '(a,es10.3,a,es10.3,a)') '# resolution,Gauss-Legendre 16 per panel; ', ppu, &
            ' panels per v_t; tail cut at ', xmax, ' v_t'
        write (u, '(a)') '# units,'//units
    end subroutine write_header

    subroutine write_b01_csv()
        integer :: u, r
        open (newunit=u, file='output/homma/B01.csv', status='replace', action='write')
        call write_header(u, 1, 'none else', 'speed_x = |v|/sqrt(2 T_b/m_b); A m s^-2; '// &
            'D m^2 s^-3; error_A, error_D: largest component change under refinement, same units')
        write (u, '(a)') 'speed_x,dir_x,dir_y,dir_z,A_x,A_y,A_z,D_xx,D_yy,D_zz,D_xy,D_xz,D_yz,'// &
            'error_A,error_D,status'
        do r = 1, n_b1
            write (u, '(15(es24.16e3,","),a)') s_b01(r), dir_b01(:, r), a_b01(:, r), &
                d_b01(1, 1, r), d_b01(2, 2, r), d_b01(3, 3, r), d_b01(1, 2, r), &
                d_b01(1, 3, r), d_b01(2, 3, r), e_b01(:, r), 'candidate'
        end do
        close (u)
    end subroutine write_b01_csv

    subroutine write_b02_csv()
        integer :: u, r
        open (newunit=u, file='output/homma/B02.csv', status='replace', action='write')
        call write_header(u, 2, 'none else', 'grad eV m^-1; speed_x = |v|/sqrt(2 T_b/m_b); '// &
            'AT m s^-2; DT m^2 s^-3; error_A, error_D: largest component change under '// &
            'refinement, same units')
        write (u, '(a)') 'grad_x_ev_m,grad_y_ev_m,grad_z_ev_m,speed_x,dir_x,dir_y,dir_z,'// &
            'AT_x,AT_y,AT_z,DT_xx,DT_yy,DT_zz,DT_xy,DT_xz,DT_yz,error_A,error_D,status'
        do r = 1, n_b2
            write (u, '(18(es24.16e3,","),a)') g_b02(:, r), s_b02(r), dir_b02(:, r), &
                a_b02(:, r), d_b02(1, 1, r), d_b02(2, 2, r), d_b02(3, 3, r), d_b02(1, 2, r), &
                d_b02(1, 3, r), d_b02(2, 3, r), e_b02(:, r), 'candidate'
        end do
        close (u)
    end subroutine write_b02_csv

end program homma_reference_main
