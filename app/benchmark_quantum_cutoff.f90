program benchmark_quantum_cutoff
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, i64 => int64
    use fortnum_interval, only: interval_t, interval, isub
    use kin6d_quantum_fock, only: assemble_jellium_pair_sector
    use kin6d_quantum_evolution, only: quantum_evolution_t, &
         initialize_quantum_evolution, solve_quantum_evolution, sample_quantum_evolution
    use kin6d_quantum_reference_gallery, only: quantum_initial_norm, &
          quantum_matrix_defect, quantum_physical_radius, quantum_observable_packet_t, &
                                 prepare_quantum_pair_sector, sample_quantum_observables
    implicit none
    integer, parameter :: times(4) = [0, 1, 4, 16], sizes(3) = [27, 125, 343]
    real(dp), parameter :: length = 6.2831853071795864769_dp
    real(dp) :: bins(16, 4, 3), radius(4, 3), refined(16, 4), unused(4)
    real(dp) :: difference(4, 2), step_difference(4), seconds(3)
    integer :: extent, t, unit, bin_unit
    open (newunit=bin_unit, file='quantum-large-cutoff-bins.csv', status='replace')
  write (bin_unit, '(a)') '# candidate arithmetic intervals do not certify trajectories'
    write (bin_unit, '(a)') 'modes,steps,time,bin,value,lower,upper,status'
    do extent = 1, 3
        call compare_sector(extent, 1024, extent == 2, bins(:, :, extent), &
                            radius(:, extent), seconds(extent))
    end do
    call compare_sector(3, 2048, .false., refined, unused, seconds(3))
    close (bin_unit)
    step_difference = maxval(abs(refined - bins(:, :, 3)), dim=1)
    difference(:, 1) = maxval(abs(bins(:, :, 1) - bins(:, :, 2)), dim=1)
    difference(:, 2) = maxval(abs(bins(:, :, 2) - bins(:, :, 3)), dim=1)
    call require(maxval(step_difference) < 1.e-5_dp, &
                 '343-mode candidate time-step refinement resolves model differences')
    call require(maxval(abs(bins(:, 1, :) - 1._dp)) < 1.e-12_dp, &
                 'larger sectors share the identical uniform pair initial state')
    open (newunit=unit, file='quantum-large-cutoff-evidence.csv', status='replace')
    write (unit, '(a)') '# cube halfwidth1/2/3 -> 27/125/343 modes; g=2,L=binary64(2pi)'
    write (unit, '(a)') '# 125 numerical/physical Hilbert certificate when radius>=0;'
  write (unit, '(a)') '# radius=-1 means unavailable; cutoff differences are indicators'
    write (unit, '(a)') 'time,max_bin_27v125,max_bin_125v343,'// &
        '343_step_refinement,max_125_bin_error_bound'
    do t = 1, 4
        write (unit, '(i0,4(",",es24.16e3))') times(t), difference(t, :), &
            step_difference(t), radius(t, 2)
        print '(a,i2,a,4es12.4)', 'time=', times(t), ' larger-cutoff indicators=', &
            difference(t, :), step_difference(t), radius(t, 2)
    end do
    close (unit)
    print *, 'PASS: physical sectors27/125/343; model/step differences are indicators'
contains
    subroutine compare_sector(extent, steps, certify, values, bounds, elapsed)
        integer, intent(in) :: extent, steps
        logical, intent(in) :: certify
        real(dp), intent(out) :: values(16, 4), bounds(4), elapsed
        integer :: count, shell, ix, iy, iz, k, j, t, cell, b
        integer, allocatable :: modes(:, :)
        real(dp), allocatable :: hr(:, :), hi(:, :), re(:), im(:)
        type(interval_t), allocatable :: hrb(:, :), hib(:, :)
        type(interval_t) :: norm0, rounding_box
        type(quantum_evolution_t) :: state
        type(quantum_observable_packet_t) :: packet
        complex(qp), allocatable :: candidate(:)
        real(dp) :: start, finish, h, radius, physical_radius, defect
        real(dp) :: result(6), lower(6), upper(6)
        real(qp) :: x, norm, ell, pi, wave, reference, voxel
        logical :: ok, bounded
        character(160) :: status
        count = sizes(extent); h = 16._dp/real(steps, dp)
        allocate (modes(3, count), hr(count, count), hi(count, count), &
           hrb(count, count), hib(count, count), re(count), im(count), candidate(count))
        k = 0
        do shell = 0, 3*extent**2
            do iz = -extent, extent
                do iy = -extent, extent
                    do ix = -extent, extent
                        if (ix*ix + iy*iy + iz*iz /= shell) cycle
                        k = k + 1; modes(:, k) = [ix, iy, iz]
                    end do
                end do
            end do
        end do
        call require(k == count, 'three-dimensional cube cardinality')
        call assemble_jellium_pair_sector(modes, length, 1._dp, 1._dp, 2._dp, &
                                          100000000_i64, hr, hi, hrb, hib, ok, status)
        call require(ok, 'strict physical symmetry constructor beyond mask capacity')
        pi = acos(-1._qp); ell = real(length, qp); wave = 2*pi/ell
        do k = 1, count
            do j = 1, count
                if (j == k) then
                    reference = wave**2*sum(modes(:, k)**2)
                else
                    reference = 8*pi/(ell**3*wave**2* &
                                      sum((modes(:, k) - modes(:, j))**2))
                end if
                call require(real(hrb(k, j)%lo, qp) <= reference .and. &
                             real(hrb(k, j)%hi, qp) >= reference, &
                         'independent real128 large-sector physical matrix containment')
            end do
        end do
        call quantum_matrix_defect(hr, hi, hrb, hib, defect, ok)
        call require(ok, 'physical coefficient packet defect')
        re = 0; im = 0; re(1) = 1
        call quantum_initial_norm(re, im, norm0, ok)
        call require(ok, 'large-sector initial norm')
        call initialize_quantum_evolution(state, hr, hi, re, im, steps, h, 1._dp, ok)
        call require(ok, 'same general shared evolution')
        call cpu_time(start)
        call solve_quantum_evolution(state, ok, certify=certify)
        call cpu_time(finish); elapsed = finish - start
        call require(ok, 'large-sector candidate solve')
        print '(a,i0,a,i0,a,l1,a,f9.3)', 'modes=', count, ' steps=', steps, &
            ' requested_certificate=', certify, ' solve_plus_requested_checker_cpu_seconds=', elapsed
        bounds = -1._dp
        do t = 1, 4
            cell = max(1, times(t)*(steps/16))
            call sample_quantum_evolution(state, cell, merge(0._dp, 1._dp, t == 1), &
                                          re, im, radius, ok, bounded)
            call require(ok, 'large-sector reconstruction at named physical time')
            candidate = cmplx(real(re, qp), real(im, qp), qp)
            norm = sqrt(sum(abs(candidate)**2))
            call require(abs(norm - 1._qp) < 1.e-5_qp, 'candidate norm diagnostic')
            if (bounded) then
                call quantum_physical_radius(norm0, radius, defect, &
                                         real(times(t), dp), 1._dp, physical_radius, ok)
                call require(ok, 'large-sector physical coefficient error ledger')
                bounds(t) = 0._dp
            else
                physical_radius = 0._dp
            end if
            call prepare_quantum_pair_sector(packet, modes, re, im, norm0, &
                                   physical_radius, 0._dp, [1._dp/16, 1._dp, 1._dp], ok)
            call require(ok, 'large sector shared voxel observable packet')
            do b = 1, 16
                x = -.5_qp + (real(b, qp) - .5_qp)/16
                voxel = pair_voxel(modes, candidate, x)
                call sample_quantum_observables(packet, real(x, dp), result, &
                                                lower, upper, ok)
                call require(ok, 'large sector shared voxel sample')
                call require(real(lower(1), qp) <= voxel .and. &
                         real(upper(1), qp) >= voxel, 'strict candidate Fourier oracle')
                values(b, t) = result(1)
                if (bounded) then
                    rounding_box = isub(interval(result(1)), interval(lower(1)))
                    bounds(t) = max(bounds(t), rounding_box%hi)
                    rounding_box = isub(interval(upper(1)), interval(result(1)))
                    bounds(t) = max(bounds(t), rounding_box%hi)
                    status = 'certified-continuous finite basis'
                else
                    status = 'candidate arithmetic only; trajectory unavailable'
                end if
                write (bin_unit, '(4(i0,","),3(es24.16e3,","),a)') &
                  count, steps, times(t), b, result(1), lower(1), upper(1), trim(status)
            end do
        end do
    end subroutine
    real(qp) function pair_voxel(modes, state, x) result(value)
        integer, intent(in) :: modes(:, :)
        complex(qp), intent(in) :: state(:)
        real(qp), intent(in) :: x
        integer :: a, b, delta(3)
        real(qp) :: phase, z, sinc, pi
        complex(qp) :: phase_factor(-6:6)
        real(qp) :: filters(-6:6)
        pi = acos(-1._qp)
        do a = -6, 6
            phase = 2*pi*a*x; z = pi*a/16; sinc = 1
            if (a /= 0) sinc = sin(z)/z
            phase_factor(a) = cmplx(cos(phase), sin(phase), qp)
            filters(a) = sinc
        end do
        value = 0
        do a = 1, size(state)
            do b = 1, size(state)
                delta = modes(:, a) - modes(:, b)
                if (any(delta(2:3) /= 0)) cycle
                value = value + real(state(a)*conjg(state(b))* &
                                     phase_factor(delta(1)), qp)*filters(delta(1))
            end do
        end do
    end function
    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: '//message
            error stop 1
        end if
    end subroutine
end program
