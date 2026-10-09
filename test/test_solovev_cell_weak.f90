program test_solovev_cell_weak
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf
    use kin6d_solovev_geometry, only: solovev_geometry_t, &
        initialize_solovev_geometry, geometry_ok, geometry_invalid, &
        geometry_interface, geometry_arithmetic
    use kin6d_solovev_cell_weak, only: solovev_cell_p2_forms
    use kin6d_grad_shafranov, only: gs_state_t, initialize_gs
    use fortfem_gauss_quadrature_2d, only: gauss_quadrature_triangle_t
    implicit none
    type(solovev_geometry_t) :: map, uninitialized
    type(gs_state_t) :: affine_state
    type(gauss_quadrature_triangle_t) :: rule
    real(dp) :: polygon(2, 8), gamma, vertex(2), vertices(2, 3), saved_vertices(2, 3)
    real(dp) :: nodes(2, 2), saved_nodes(2, 2), weights(2), matrix(6, 6), load(6)
    real(dp) :: want_matrix(6, 6), want_load(6), coordinates(6, 2), want_energy
    real(dp) :: physical_points(2, 2), want_jacobians(2, 2, 2)
    real(dp) :: native_point(2), jg(2, 2), det, area, phi(6), reference(2, 6)
    real(dp) :: inverse(2, 2), gradient(2, 6), wrong_gradient(2, 6)
    real(dp) :: wrong_matrix(6, 6), wrong_load(6), polygon_load(6), no_det_load(6)
    real(dp) :: sum_load, mu0, pp, ff, source, old_source, ar(3), v(6)
    real(dp) :: affine_vertices(2, 3), qnodes(2, 7), nan, infinity
    integer :: status, passed, i, j, k, q
    logical :: certified
    passed = 0
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    infinity = ieee_value(0.0_dp, ieee_positive_inf)
    gamma = 1.44_dp*(sqrt(2.0_dp)-1)
    do i = 1, 8
        select case (i)
        case (1); vertex = [-gamma, -1.44_dp]
        case (2); vertex = [gamma, -1.44_dp]
        case (3); vertex = [1.44_dp, -gamma]
        case (4); vertex = [1.44_dp, gamma]
        case (5); vertex = [gamma, 1.44_dp]
        case (6); vertex = [-gamma, 1.44_dp]
        case (7); vertex = [-1.44_dp, gamma]
        case (8); vertex = [-1.44_dp, -gamma]
        end select
        polygon(:, i) = [3.0_dp, 0.0_dp] &
            +[4*vertex(1)-3*vertex(2), 3*vertex(1)+4*vertex(2)]/25
    end do
    call initialize_solovev_geometry(map, 3.0_dp, 5.3145_dp, 0.75_dp, polygon, status)
    call require(status == geometry_ok, 'noncircular map initialized')
    mu0 = 2; pp = -7.0_dp/5; ff = 3.0_dp/7
    weights = [0.125_dp, 0.375_dp]
    ! Independent rational implicit equations at t=4/5 and801/1000 on the
    ! boundary ray (3/10,1/5), facet ratio4/5; no emitted geometry/basis used.
    ! Constants below are rounded exact rational sums, not sampled fit data.
    nodes = reshape([ &
        2.50000000000000000e-01_dp, 1.25000000000000000e-01_dp, &
        2.56433889523809522e-01_dp, 1.33773485714285711e-01_dp ], [2, 2])
    vertices = reshape([ &
        3.19959166666666661e+00_dp, 1.34091927083333340e-01_dp, &
        3.21521666666666661e+00_dp, 1.37998177083333340e-01_dp, &
        3.21000833333333313e+00_dp, 1.45810677083333340e-01_dp ], [2, 3])
    want_matrix = reshape([ &
        1.51908192728213515e-01_dp, -1.05661068248556322e-03_dp, &
        2.32688926908400998e-02_dp, -2.93941653900748057e-02_dp, &
        -7.91324215648154666e-02_dp, -6.55938877816777677e-02_dp, &
        -1.05661068248556322e-03_dp, 1.74784089605332986e-04_dp, &
        2.81610351141248043e-03_dp, 1.58126507802039892e-02_dp, &
        -2.57913760881641960e-03_dp, -1.51677900899198186e-02_dp, &
        2.32688926908400998e-02_dp, 2.81610351141248043e-03_dp, &
        7.80588535325365818e-02_dp, 3.81221511748546682e-01_dp, &
        -8.95855561118004406e-02_dp, -3.95779805371535398e-01_dp, &
        -2.93941653900748057e-02_dp, 1.58126507802039892e-02_dp, &
        3.81221511748546682e-01_dp, 2.00392376531669258e+00_dp, &
        -3.86071232014101484e-01_dp, -1.98549253044126695e+00_dp, &
        -7.91324215648154666e-02_dp, -2.57913760881641960e-03_dp, &
        -8.95855561118004406e-02_dp, -3.86071232014101484e-01_dp, &
        1.21885760242350288e-01_dp, 4.35482587057183512e-01_dp, &
        -6.55938877816777677e-02_dp, -1.51677900899198186e-02_dp, &
        -3.95779805371535398e-01_dp, -1.98549253044126695e+00_dp, &
        4.35482587057183512e-01_dp, 2.02655142662721621e+00_dp ], [6, 6])
    want_load = reshape([ &
        -1.62734375812422107e-04_dp, 1.45764223701419688e-04_dp, &
        1.13082446626230939e-04_dp, -7.29608817641124916e-04_dp, &
        -1.56529445199780367e-04_dp, -3.76667543530382232e-04_dp ], [6])
    coordinates = reshape([ &
        0.00000000000000000e+00_dp, 2.36302600197204615e-02_dp, &
        1.80214875457270392e-02_dp, 1.14892658434416119e-02_dp, &
        2.02691018553802192e-02_dp, 8.77983601193860728e-03_dp, &
        0.00000000000000000e+00_dp, 8.03144149136813916e-03_dp, &
        1.76839686835390635e-02_dp, 3.84481329666940810e-03_dp, &
        1.25079410211097631e-02_dp, 8.66312772444035629e-03_dp ], [6, 2])
    want_energy = 4.02741509680428164e-05_dp
    physical_points = reshape([ &
        3.24000000000000021e+00_dp, 1.60000000000000003e-01_dp, &
        3.24029999999999996e+00_dp, 1.60200000000000009e-01_dp ], [2, 2])
    want_jacobians = reshape([ &
        2.29785316868832239e-02_dp, 7.68962659333881620e-03_dp, &
        1.73287642629523023e-02_dp, 1.71473988315515337e-02_dp, &
        2.29953042792552693e-02_dp, 7.69842339050486740e-03_dp, &
        1.73449711904866410e-02_dp, 1.71599523993678782e-02_dp ], [2, 2, 2])

    saved_vertices = vertices; saved_nodes = nodes
    call evaluate(map)
    call require(status == geometry_ok .and. .not. certified, &
        'uncertified forms succeed')
    call close(maxval(abs(matrix-want_matrix)), 0.0_dp, 2.0e-10_dp, &
        'two-node full stiffness rational oracle')
    call close(maxval(abs(load-want_load)), 0.0_dp, 2.0e-14_dp, &
        'two-node signed weak source rational oracle')
    call structure()
    ! These reference P2 coefficient vectors match the physical coordinate
    ! gradients at ONLY the two declared nodes. No map/P2 completeness claim.
    call close(dot_product(coordinates(:, 1),matmul(matrix,coordinates(:, 1))), &
        want_energy, 2.0e-15_dp, 'physical radial coordinate finite-node energy')
    call close(dot_product(coordinates(:, 2),matmul(matrix,coordinates(:, 2))), &
        want_energy, 2.0e-15_dp, 'physical vertical coordinate finite-node energy')
    call close(dot_product(coordinates(:, 1),matmul(matrix,coordinates(:, 2))), &
        0.0_dp, 2.0e-15_dp, 'physical coordinate cross energy')
    ! Behavioral mutants use different explicit barycentric/cofactor oracles.
    wrong_matrix = 0; wrong_load = 0; polygon_load = 0; no_det_load = 0; sum_load = 0
    do q = 1, 2
        call independent_shapes(nodes(:, q), phi, reference)
        jg = want_jacobians(:, :, q)
        det = jg(1,1)*jg(2,2)-jg(1,2)*jg(2,1)
        inverse = reshape([jg(2,2),-jg(2,1),-jg(1,2),jg(1,1)], [2,2])/det
        gradient = matmul(transpose(inverse),reference)
        wrong_gradient = matmul(inverse,reference)
        native_point = vertices(:,1)+(vertices(:,2)-vertices(:,1))*nodes(1,q) &
            +(vertices(:,3)-vertices(:,1))*nodes(2,q)
        source = mu0*pp*physical_points(1,q)+ff/physical_points(1,q)
        old_source = mu0*pp*native_point(1)+ff/native_point(1)
        area = (vertices(1,2)-vertices(1,1))*(vertices(2,3)-vertices(2,1)) &
            -(vertices(1,3)-vertices(1,1))*(vertices(2,2)-vertices(2,1))
        sum_load = sum_load+weights(q)*det*source
        do j = 1, 6
            wrong_load(j) = wrong_load(j)+weights(q)*det*old_source*phi(j)
            polygon_load(j) = polygon_load(j)+weights(q)*area*old_source*phi(j)
            no_det_load(j) = no_det_load(j)+weights(q)*source*phi(j)
            do i = 1, 6
                wrong_matrix(i,j) = wrong_matrix(i,j)+weights(q)*det &
                    /physical_points(1,q)*dot_product(wrong_gradient(:,i), &
                    wrong_gradient(:,j))
            end do
        end do
    end do
    call close(sum(load), sum_load, 3.0e-14_dp, &
        'unit partition load and physical measure')
    call require(maxval(abs(matrix-wrong_matrix)) > 0.01_dp, &
        'inverse transpose mutant')
    call require(maxval(abs(load-wrong_load)) > 1.0e-6_dp, &
        'wrong physical R source mutant')
    call require(maxval(abs(load-polygon_load)) > 1.0e-4_dp, &
        'retained polygon load mutant')
    call require(maxval(abs(load-no_det_load)) > 0.1_dp, &
        'omitted detJ load mutant')
    weights = 2*weights
    call evaluate(map)
    call close(maxval(abs(matrix-2*want_matrix)), 0.0_dp, 4.0e-10_dp, &
        'caller measure preserved without automatic normalization')
    weights = weights/2
    ! Existing affine P2 engine, same native seven-point nodes/Delta weights.
    affine_vertices = reshape([2.98_dp,-0.02_dp,3.02_dp,-0.02_dp, &
        2.98_dp,0.02_dp], [2,3])
    call initialize_gs(affine_state,affine_vertices,1.0_dp,status,element_order=2)
    call require(status == 0, 'native affine P2 cell initialized without solve')
    call require(size(affine_state%triangles,2) == 1, 'single affine parity cell')
    affine_vertices = affine_state%points(:,affine_state%triangles(1:3,1))
    call rule%init(5)
    qnodes(1,:) = rule%xi; qnodes(2,:) = rule%eta
    call solovev_cell_p2_forms(map,affine_vertices,qnodes,rule%weights,mu0,pp, &
        0.0_dp,matrix,load,status,certified)
    call require(status == geometry_ok .and. .not. certified, 'identity region forms')
    call close(maxval(abs(matrix-affine_state%stiffness(:,:,1))),0.0_dp, &
        8.0e-13_dp,'identity geometry native affine stiffness parity')
    ar = affine_vertices(1,:)
    do i = 1, 3
        want_load(i) = mu0*pp*(3*ar(i)-sum(ar))/120
    end do
    want_load(4) = mu0*pp*(2*ar(1)+2*ar(2)+ar(3))/30
    want_load(5) = mu0*pp*(ar(1)+2*ar(2)+2*ar(3))/30
    want_load(6) = mu0*pp*(2*ar(1)+ar(2)+2*ar(3))/30
    want_load = want_load*2*affine_state%area(1)
    call close(maxval(abs(load-want_load)),0.0_dp,2.0e-16_dp, &
        'identity polynomial load exact barycentric moments')
    ! Later-node interface must discard the already accumulated first node.
    nodes(:,1) = [0.25_dp,0.125_dp]; nodes(:,2) = [0.5_dp,0.125_dp]
    vertices(:,1) = [3.075_dp,0.05_dp] &
        -0.5_dp*[0.3_dp,0.2_dp]-0.125_dp*[-0.025_dp,0.0375_dp]
    vertices(:,2) = vertices(:,1)+[0.3_dp,0.2_dp]
    vertices(:,3) = vertices(:,1)+[-0.025_dp,0.0375_dp]
    call evaluate(map)
    call rejected(geometry_interface, 'late inner interface')
    vertices = saved_vertices; nodes = saved_nodes
    nodes(1,2) = nan; call evaluate(map)
    call rejected(geometry_invalid, 'late NaN node')
    nodes = saved_nodes; nodes(:,2) = [0.8_dp,0.8_dp]; call evaluate(map)
    call rejected(geometry_invalid, 'outside Delta')
    nodes = saved_nodes
    nodes(:,2) = [nearest(0.5_dp, 1.0_dp), 0.5_dp]
    call evaluate(map)
    call rejected(geometry_invalid, 'late rounded-sum outside Delta')
    nodes = saved_nodes; weights(2) = 0; call evaluate(map)
    call rejected(geometry_invalid, 'zero weight')
    weights(2) = -0.375_dp; call evaluate(map)
    call rejected(geometry_invalid, 'negative weight')
    weights(2) = infinity; call evaluate(map)
    call rejected(geometry_invalid, 'nonfinite weight')
    weights(2) = 0.375_dp; mu0 = 0; call evaluate(map)
    call rejected(geometry_invalid, 'nonpositive mu0')
    mu0 = 2; pp = nan; call evaluate(map)
    call rejected(geometry_invalid, 'nonfinite source input')
    pp = -7.0_dp/5; ff = infinity; call evaluate(map)
    call rejected(geometry_invalid, 'nonfinite FFprime')
    ff = 3.0_dp/7; call evaluate(uninitialized)
    call rejected(geometry_invalid, 'uninitialized map')
    vertices(:,3) = vertices(:,2); call evaluate(map)
    call rejected(geometry_invalid, 'degenerate cell')
    vertices = saved_vertices; pp = huge(1.0_dp); call evaluate(map)
    call rejected(geometry_arithmetic, 'source overflow')
    pp = -7.0_dp/5
    weights(2) = huge(1.0_dp)
    call evaluate(map)
    call rejected(geometry_arithmetic, 'late stiffness accumulation overflow')
    weights(2) = 0.375_dp
    call solovev_cell_p2_forms(map,vertices,nodes(:,:1),weights,mu0,pp,ff, &
        matrix,load,status,certified)
    call rejected(geometry_invalid,'node weight shape mismatch')
    call solovev_cell_p2_forms(map,vertices,nodes(:1,:),weights,mu0,pp,ff, &
        matrix,load,status,certified)
    call rejected(geometry_invalid,'node component dimension')
    call solovev_cell_p2_forms(map,vertices,nodes(:,:0),weights(:0),mu0,pp,ff, &
        matrix,load,status,certified)
    call rejected(geometry_invalid,'empty quadrature')
    print *, 'Solovev mapped P2 cell weak-form checks:',passed
contains
    subroutine evaluate(selected_map)
        type(solovev_geometry_t), intent(in) :: selected_map
        call solovev_cell_p2_forms(selected_map,vertices,nodes,weights,mu0,pp,ff, &
            matrix,load,status,certified)
    end subroutine
    subroutine structure()
        call close(maxval(abs(matrix-transpose(matrix))),0.0_dp,0.0_dp,'symmetric sum')
        call close(maxval(abs(sum(matrix,dim=1))), &
            0.0_dp,8.0e-13_dp,'constant annihilation')
        do k = 1, 6
            do j = 1, 6
                v(j) = real((-1)**(j+k)*(j+2*k),dp)
            end do
            call require(dot_product(v,matmul(matrix,v)) >= -1.0e-12_dp, &
                'positive finite-node energy')
        end do
    end subroutine
    subroutine independent_shapes(x, shape, g)
        real(dp), intent(in) :: x(2)
        real(dp), intent(out) :: shape(6), g(2,6)
        real(dp) :: l, a, b
        a=x(1); b=x(2); l=1-a-b
        shape=[l*(2*l-1),a*(2*a-1),b*(2*b-1),4*l*a,4*a*b,4*l*b]
        g(1,:)=[1-4*l,4*a-1,0.0_dp,4*(l-a),4*b,-4*b]
        g(2,:)=[1-4*l,0.0_dp,4*b-1,-4*a,4*a,4*(l-b)]
    end subroutine
    subroutine rejected(expected,label)
        integer, intent(in) :: expected
        character(*), intent(in) :: label
        call require(status == expected,label//' status')
        call require(all(matrix == 0) .and. all(load == 0),label//' clears forms')
        call require(.not. certified,label//' uncertified')
    end subroutine
    subroutine require(ok,label)
        logical, intent(in) :: ok
        character(*), intent(in) :: label
        if (.not. ok) then
            print *, 'FAIL: ',label
            error stop 1
        end if
        passed=passed+1
    end subroutine
    subroutine close(actual,expected,tolerance,label)
        real(dp), intent(in) :: actual,expected,tolerance
        character(*), intent(in) :: label
        call require(abs(actual-expected) <= tolerance,label)
    end subroutine
end program
