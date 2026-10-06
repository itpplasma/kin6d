!> Figures of the H1a correlation map (docs/figure-style.md). Display only:
!> every plotted value and band comes from the CSV/NetCDF data.
module correlation_map_plots
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortplot, only: figure, subplot, plot, fill_between, title, suptitle, xlabel, &
        ylabel, xlim, ylim, legend, tight_layout, savefig
    use kin6d_demo_style, only: blue, vermillion, green, orange, purple, sky, black, &
        grey, num
    use kin6d_ocp_hnc, only: hnc_t
    use kin6d_correlation_hierarchy, only: ocp_energy_literature
    use correlation_map_md, only: md_result_t
    use correlation_map_types
    implicit none
    private
    public :: plot_all
contains
    subroutine plot_all(chains, ng, hnc, md)
        type(chain_t), intent(in) :: chains(:)
        integer, intent(in) :: ng
        type(hnc_t), intent(in) :: hnc(:)
        type(md_result_t), intent(in) :: md
        call plot_pairs(chains(1:ng))
        call plot_triplets(chains(1:ng))
        call plot_sizes(chains, ng)
        call plot_ybg(chains(1:ng))
        call plot_energy(chains, ng, md, hnc)
    end subroutine

    character(16) function gamma_label(gamma)
        real(dp), intent(in) :: gamma
        gamma_label = '$\Gamma$ = '//num(gamma, 1)
    end function

    subroutine plot_pairs(c)
        type(chain_t), intent(in) :: c(:)
        integer :: i, m
        call figure(figsize=[13.2_dp, 8.8_dp])
        do i = 1, min(6, size(c))
            call subplot(2, 3, i)
            m = count(c(i)%r <= 6)
            call fill_between(c(i)%r(:m), c(i)%g(:m) - 2*c(i)%g_err(:m), &
                              c(i)%g(:m) + 2*c(i)%g_err(:m), color=blue, alpha=0.25_dp)
            call plot(c(i)%r(:m), c(i)%g(:m), color=blue, label='MC $\pm 2\sigma$')
            call plot(c(i)%r(:m), c(i)%g_hnc(:m), color=vermillion, linestyle='-.', &
                      label='HNC')
            call plot(c(i)%r(:m), max(-0.2_dp, c(i)%g_dhl(:m)), color=black, &
                      linestyle=':', label='DH linear')
            call plot(c(i)%r(:m), c(i)%g_dh(:m), color=black, linestyle='--', &
                      label='Debye-Hückel')
            call xlim(0.0_dp, 6.0_dp); call ylim(-0.2_dp, max(1.3_dp, &
                maxval(c(i)%g(:m)) + 0.1_dp))
            call title(trim(gamma_label(c(i)%gamma)))
            call xlabel('r / a'); call ylabel('$g_2(r)$')
            if (i == 1) call legend()
        end do
        call suptitle('OCP pair distribution: Metropolis MC (N = '//num(real(c(1)%n, dp), 3) &
                      //') vs reduced models')
        call tight_layout()
        call savefig('output/correlation-map-g2.png')
    end subroutine

    subroutine plot_triplets(c)
        type(chain_t), intent(in) :: c(:)
        integer :: i
        call figure(figsize=[13.2_dp, 8.8_dp])
        do i = 1, min(6, size(c))
            call subplot(2, 3, i)
            call fill_between(c(i)%eq_r, c(i)%eq_g3 - 2*c(i)%eq_g3_err, &
                              c(i)%eq_g3 + 2*c(i)%eq_g3_err, color=blue, alpha=0.25_dp)
            call plot(c(i)%eq_r, c(i)%eq_g3, color=blue, marker='o', markersize=4.0_dp, &
                      label='MC $g_3(r,r,r)$')
            call plot(c(i)%eq_r, c(i)%eq_lin, color=vermillion, linestyle=':', &
                      label='$1+\sum h$ (cumulant cut)')
            call plot(c(i)%eq_r, c(i)%eq_kirk, color=black, linestyle='--', &
                      label='Kirkwood $g g g$')
            call xlim(0.0_dp, c(i)%eq_r(size(c(i)%eq_r)) + 0.1_dp)
            call title(trim(gamma_label(c(i)%gamma)))
            call xlabel('r / a (equilateral)'); call ylabel('$g_3$')
            if (i == 1) call legend()
        end do
        call suptitle('Equilateral triplet distribution vs superposition and cumulant cut')
        call tight_layout()
        call savefig('output/correlation-map-triplet.png')
    end subroutine

    !> log10 of sqrt(squared norm) with a 1-sigma band; unresolved values
    !> (below 2 sigma) are shown as upper limits at value + 2 sigma.
    subroutine norm_series(c, idx, colour, marker, label)
        type(chain_t), intent(in) :: c(:)
        integer, intent(in) :: idx
        real(dp), intent(in) :: colour(3)
        character(*), intent(in) :: marker, label
        real(dp) :: x(size(c)), y(size(c)), e(size(c))
        integer :: i
        x = log10(c%gamma)
        do i = 1, size(c)
            y(i) = c(i)%out(idx); e(i) = c(i)%out_err(idx)
        end do
        where (y > 2*e)
            y = log10(sqrt(y))
        elsewhere
            y = log10(sqrt(max(y, 0.0_dp) + 2*e))
        end where
        call plot(x, y, color=colour, marker=marker, markersize=6.0_dp, label=label)
        do i = 1, size(c)
            if (c(i)%out(idx) > 2*e(i)) then
                call plot([x(i), x(i)], [log10(sqrt(c(i)%out(idx) - min(e(i), &
                          c(i)%out(idx)/2))), log10(sqrt(c(i)%out(idx) + e(i)))], &
                          color=colour)
            else
                call plot([x(i)], [y(i)], color=colour, linestyle='none', marker='v', &
                          markersize=8.0_dp)
            end if
        end do
    end subroutine

    subroutine plot_sizes(chains, ng)
        type(chain_t), intent(in) :: chains(:)
        integer, intent(in) :: ng
        real(dp) :: x(ng), y(ng), e(ng), p, a, pe, xf(2)
        integer :: used, i, weak
        call figure(figsize=[12.8_dp, 6.4_dp])
        call subplot(1, 2, 1)
        call norm_series(chains(1:ng), O_H2, blue, 'o', '$\|h_2\|$')
        call norm_series(chains(1:ng), O_RED, green, 's', '$\|K - P\|$ (pair products)')
        call norm_series(chains(1:ng), O_U3, orange, '^', '$\|U_3\|$ (Ursell)')
        call norm_series(chains(1:ng), O_C3, vermillion, 'D', '$\|c_3\|$ (connected)')
        call norm_series(chains(ng + 1:), O_C3, purple, 'd', '$\|c_3\|$, small N')
        x = log10(chains(1:ng)%gamma)
        do i = 1, ng
            y(i) = log10(sqrt(max(chains(i)%bin_indicator, 1.0e-300_dp)))
        end do
        call plot(x, y, color=grey, linestyle=':', label='bin indicator')
        weak = count(chains(1:ng)%gamma <= 1.0_dp)
        do i = 1, ng
            y(i) = sqrt(max(chains(i)%out(O_H2), 0.0_dp))
            e(i) = chains(i)%out_err(O_H2)/(2*max(y(i), 1.0e-300_dp))
        end do
        call power_fit(chains(1:weak)%gamma, y(1:weak), e(1:weak), p, a, pe, used)
        xf = [x(1), x(weak)]
        if (used >= 2) call plot(xf, log10(a) + p*xf, color=blue, linestyle='--', &
                                 label='fit $\Gamma^{'//num(p, 2)//'}$')
        call xlabel('$\log_{10} \Gamma$'); call ylabel('$\log_{10}$ size')
        call title('Correlation sizes (n-weighted $L^2$, sides < 3a)')
        call legend()
        call subplot(1, 2, 2)
        do i = 1, ng
            y(i) = log10(max(chains(i)%out(O_DK), 1.0e-12_dp))
        end do
        call plot(x, y, color=blue, marker='o', markersize=6.0_dp, &
                  label='rms Kirkwood defect $\delta_K$')
        do i = 1, ng
            y(i) = 0.5_dp*log10(max(chains(i)%out(O_C3), 1.0e-24_dp)/ &
                                max(chains(i)%out(O_RED), 1.0e-24_dp))
        end do
        call plot(x, y, color=vermillion, marker='D', markersize=6.0_dp, &
                  label='$\|c_3\| / \|K - P\|$')
        call plot([x(1), x(ng)], [0.0_dp, 0.0_dp], color=black, linestyle='--', &
                  label='equal size')
        call xlabel('$\log_{10} \Gamma$'); call ylabel('$\log_{10}$ ratio')
        call title('Beyond superposition')
        call legend()
        call tight_layout()
        call savefig('output/correlation-map-size.png')
    end subroutine

    subroutine plot_ybg(c)
        type(chain_t), intent(in) :: c(:)
        real(dp) :: x(size(c)), y(size(c))
        integer :: i, j, m, pick(2)
        integer, parameter :: idx(5) = [O_EEX, O_EK, O_EL, O_EHNC, O_EDH]
        integer, parameter :: uidx(4) = [O_UEX, O_UK, O_UL, 0]
        character(24), parameter :: names(5) = [character(24) :: 'YBG, exact force', &
            'YBG, Kirkwood', 'YBG, first order', 'HNC', 'Debye-Hückel']
        real(dp) :: colours(3, 5)
        character(4) :: markers(5) = ['o', 's', '^', 'D', 'v']
        colours = reshape([blue, vermillion, green, orange, black], [3, 5])
        x = log10(c%gamma)
        call figure(figsize=[13.2_dp, 8.8_dp])
        call subplot(2, 2, 1)
        do j = 1, 5
            do i = 1, size(c)
                y(i) = log10(max(c(i)%out(idx(j)), 1.0e-6_dp))
            end do
            call plot(x, y, color=colours(:, j), marker=trim(markers(j)), &
                      markersize=6.0_dp, label=trim(names(j)))
        end do
        call xlabel('$\log_{10} \Gamma$'); call ylabel('$\log_{10} \|g - g_{MC}\|$')
        call title('Pair-structure error per closure')
        call legend()
        call subplot(2, 2, 2)
        do j = 1, 5
            do i = 1, size(c)
                if (j <= 3) then
                    y(i) = abs(c(i)%out(uidx(j)) - c(i)%u)
                else if (j == 4) then
                    y(i) = abs(c(i)%u_hnc - c(i)%u)
                else
                    y(i) = abs(c(i)%u_dh - c(i)%u)
                end if
                y(i) = log10(max(y(i)/abs(c(i)%u), 1.0e-7_dp))
            end do
            call plot(x, y, color=colours(:, j), marker=trim(markers(j)), &
                      markersize=6.0_dp, label=trim(names(j)))
        end do
        do i = 1, size(c)
            y(i) = log10(max(c(i)%u_err, c(i)%u_batch_err)/abs(c(i)%u))
        end do
        call plot(x, y, color=grey, linestyle=':', label='MC 1$\sigma$')
        call xlabel('$\log_{10} \Gamma$'); call ylabel('$\log_{10} |\Delta u| / |u|$')
        call title('Excess energy error per closure')
        call legend()
        pick = [min(size(c), 3), size(c) - 1]
        do j = 1, 2
            i = pick(j)
            call subplot(2, 2, 2 + j)
            m = count(c(i)%r <= 5)
            call plot(c(i)%r(:m), c(i)%g(:m), color=blue, label='MC')
            call plot(c(i)%r(:m), c(i)%g_kirk(:m), color=vermillion, linestyle='-.', &
                      label='YBG, Kirkwood')
            call plot(c(i)%r(:m), c(i)%g_lin(:m), color=green, linestyle=':', &
                      label='YBG, first order')
            call plot(c(i)%r(:m), c(i)%g_exact(:m), color=black, linestyle='--', &
                      label='YBG, exact force')
            call xlim(0.0_dp, 5.0_dp)
            call title('$g_2$ from the first YBG equation, '//trim(gamma_label(c(i)%gamma)))
            call xlabel('r / a'); call ylabel('$g_2(r)$')
            call legend()
        end do
        call tight_layout()
        call savefig('output/correlation-map-ybg.png')
    end subroutine

    subroutine plot_energy(chains, ng, md, hnc)
        type(chain_t), intent(in) :: chains(:)
        integer, intent(in) :: ng
        type(md_result_t), intent(in) :: md
        type(hnc_t), intent(in) :: hnc(:)
        real(dp) :: x(ng), y(ng), e(ng), ul, dl, cols(3, 3)
        integer :: i, s, k, m, pk(3)
        call figure(figsize=[13.2_dp, 8.8_dp])
        call subplot(2, 2, 1)
        do s = 0, 1
            do i = 1, ng
                associate (c => chains(s*ng + i))
                    call ocp_energy_literature(c%gamma, ul, dl)
                    x(i) = log10(c%gamma) + 0.02_dp*s
                    y(i) = (c%u - ul)/abs(ul)
                    e(i) = max(c%u_err, c%u_batch_err)/abs(ul)
                end associate
                call plot([x(i), x(i)], [y(i) - e(i), y(i) + e(i)], &
                          color=merge(blue, vermillion, s == 0))
            end do
            call plot(x, y, color=merge(blue, vermillion, s == 0), linestyle='none', &
                      marker=merge('o', 's', s == 0), markersize=6.0_dp, &
                      label='MC, N = '//num(real(chains(s*ng + 1)%n, dp), 3))
        end do
        call ocp_energy_literature(md%gamma, ul, dl)
        call plot([log10(md%gamma) - 0.02_dp], [(md%u - ul)/abs(ul)], color=green, &
                  linestyle='none', marker='^', markersize=8.0_dp, label='MD, N = '// &
                  num(real(md%n, dp), 3))
        call plot([x(1), x(ng)], [0.0_dp, 0.0_dp], color=black, linestyle='--', &
                  label='Potekhin-Chabrier fit')
        call xlabel('$\log_{10} \Gamma$'); call ylabel('$(u - u_{fit}) / |u_{fit}|$')
        call title('Excess energy vs literature')
        call legend()
        call subplot(2, 2, 2)
        do i = 1, ng
            call ocp_energy_literature(chains(i)%gamma, ul, dl)
            x(i) = log10(chains(i)%gamma)
            y(i) = (chains(i)%dudg - dl)/abs(dl); e(i) = chains(i)%dudg_err/abs(dl)
            call plot([x(i), x(i)], [y(i) - e(i), y(i) + e(i)], color=blue)
        end do
        call plot(x, y, color=blue, linestyle='none', marker='o', markersize=6.0_dp, &
                  label='fluctuation $(u - Var(\beta U)/N)/\Gamma$')
        call plot([x(1), x(ng)], [0.0_dp, 0.0_dp], color=black, linestyle='--', &
                  label='fit derivative (JVP leaf)')
        call xlabel('$\log_{10} \Gamma$'); call ylabel('relative difference')
        call title('$du/d\Gamma$: fluctuation formula vs fit')
        call legend()
        call subplot(2, 2, 3)
        k = min(3, ng)
        m = count(chains(k)%r <= 6)
        call plot(chains(k)%r(:m), chains(k)%g(:m), color=blue, label='MC')
        m = count(md%r <= 6)
        call plot(md%r(:m), md%g(:m), color=vermillion, linestyle='--', label='MD (NVE)')
        call xlim(0.0_dp, 6.0_dp)
        call title('MC vs MD, '//trim(gamma_label(chains(k)%gamma)))
        call xlabel('r / a'); call ylabel('$g_2(r)$')
        call legend()
        call subplot(2, 2, 4)
        pk = [1, min(3, ng), min(5, ng)]
        cols = reshape([blue, vermillion, green], [3, 3])
        do s = 1, 3
            associate (c => chains(pk(s)))
                call plot(c%s_k, c%s_mc, color=cols(:, s), linestyle='none', &
                          marker='o', markersize=5.0_dp, &
                          label='MC '//trim(gamma_label(c%gamma)))
                call plot(hnc(pk(s))%k(:60), hnc(pk(s))%s(:60), color=cols(:, s), &
                          linestyle='-.')
                call plot(c%s_k, c%s_dh, color=black, linestyle=':')
            end associate
        end do
        call xlim(0.0_dp, 2.5_dp)
        call title('S(k): MC (points), HNC (dash-dot), DH (dotted)')
        call xlabel('k a'); call ylabel('S(k)')
        call legend()
        call tight_layout()
        call savefig('output/correlation-map-energy.png')
    end subroutine
end module correlation_map_plots
