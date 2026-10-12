! Shared leaf-writing support for every KIN6D generator.
! One emit path (primal, JVP/VJP and strict interval leaves) and one FortSym
! revision source: KIN6D_FORTSYM_REVISION in CMakeLists.txt, delivered as the
! configured include file kin6d_fortsym_revision.inc.
module kin6d_gen_support
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym
    use fortsym_products, only: jvp, vjp
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_string, only: str_t, str, chars
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_batch_ir
    use fortsym_rigorous_emit, only: rigorous_kernel_spec_t, rigorous_runtime_t, &
                                     interval_runtime, emit_rigorous_kernel
    implicit none
    private
    include 'kin6d_fortsym_revision.inc'
    public :: fortsym_revision, leaf_writer_t

    !> Writes one generated leaf module. The generator owns the arena; the
    !> writer owns the output unit and the leaf conventions.
    type :: leaf_writer_t
        type(arena_t), pointer :: arena => null()
        type(native_engine_t) :: engine
        type(rigorous_runtime_t) :: runtime
        character(:), allocatable :: generator
        integer :: unit = -1
        !> Simplify primal roots before emission (used by gen_debye_sphere).
        logical :: simplify_roots = .false.
    contains
        procedure :: open => writer_open
        procedure :: line => writer_line
        procedure :: close => writer_close
        procedure :: emit => writer_emit
        procedure :: emit_batch => writer_emit_batch
    end type
contains
    subroutine writer_open(w, arena, generator, path)
        class(leaf_writer_t), intent(inout) :: w
        type(arena_t), intent(inout), target :: arena
        character(*), intent(in) :: generator, path

        w%arena => arena
        w%engine = make_native_engine(arena)
        w%generator = generator
        w%runtime = interval_runtime('fortnum_interval', 'interval_t')
        w%runtime%sin = str('sin')
        w%runtime%cos = str('cos')
        w%runtime%abs = str('abs')
        w%runtime%sincos = str('sincos')
        w%runtime%exp = str('exp')
        open (newunit=w%unit, file=path, status='replace', action='write')
    end subroutine

    subroutine writer_line(w, text)
        class(leaf_writer_t), intent(in) :: w
        character(*), intent(in) :: text
        write (w%unit, '(a)') text
    end subroutine

    subroutine writer_close(w)
        class(leaf_writer_t), intent(inout) :: w
        close (w%unit)
        w%unit = -1
    end subroutine

    !> Emit a floating leaf; optionally its JVP/VJP (ad), a strict interval
    !> twin (interval), only the interval leaf (rigorous), or an interval
    !> twin of the VJP (vjp_interval). quad emits a real(real128) floating
    !> leaf; the leaf module must use iso_fortran_env real128.
    recursive subroutine writer_emit(w, name, args, outputs, roots, ad, interval, &
                                     rigorous, vjp_interval, quad)
        class(leaf_writer_t), intent(inout) :: w
        character(*), intent(in) :: name, args(:), outputs(:)
        type(expr_t), intent(in) :: roots(:)
        logical, intent(in), optional :: ad, interval, rigorous, vjp_interval, quad
        type(expr_t), allocatable :: clean(:)
        type(engine_result_t) :: simplified
        type(kernel_spec_t) :: spec
        type(rigorous_kernel_spec_t) :: rspec
        type(str_t) :: source
        logical :: ok
        character(:), allocatable :: message
        integer :: j

        clean = roots
        if (w%simplify_roots) then
            do j = 1, size(roots)
                simplified = w%engine%simplify(roots(j))
                if (.not. simplified%ok) error stop chars(simplified%message)
                clean(j) = simplified%value
            end do
        end if
        if (flag(rigorous)) then
            rspec%name = str(name)
            rspec%generator = str(w%generator//'; FortSym '//fortsym_revision)
            rspec%runtime = w%runtime
            rspec%elemental_procedure = .true.
            allocate (rspec%args(size(args)), rspec%outputs(size(outputs)))
            do j = 1, size(args)
                rspec%args(j) = str(trim(args(j)))
            end do
            do j = 1, size(outputs)
                rspec%outputs(j) = str(trim(outputs(j)))
            end do
            source = emit_rigorous_kernel(clean, rspec, ok, message)
            if (.not. ok) error stop message
            write (w%unit, '(a)') chars(source)
            return
        end if
        spec%name = str(name)
        spec%generator = str(w%generator)
        spec%generator_revision = str(fortsym_revision)
        spec%regenerate_command = str(w%generator//' from pinned FortSym build')
        spec%pure_procedure = .true.
        spec%elemental_procedure = .true.
        spec%openmp_declare_target = .false.
        spec%openacc_routine_seq = .false.
        if (flag(quad)) spec%scalar_type = str('real(real128)')
        allocate (spec%args(size(args)), spec%outputs(size(outputs)))
        do j = 1, size(args)
            spec%args(j) = str(trim(args(j)))
        end do
        do j = 1, size(outputs)
            spec%outputs(j) = str(trim(outputs(j)))
        end do
        source = emit_kernel(clean, spec, ok, message=message)
        if (.not. ok) error stop message
        write (w%unit, '(a)') chars(source)
        if (flag(ad)) call emit_derivatives(w, name, args, outputs, clean, &
                                            flag(vjp_interval))
        if (flag(interval)) call w%emit(name//'_interval', args, outputs, roots, &
                                        rigorous=.true.)
    end subroutine

    !> Emit `<name>_batch_mod`: a `declare simd` leaf and a batched host driver
    !> over rank-one arrays (FortSym FY1, SIMD mode). `uniform` names the
    !> arguments that are scalars across the batch. The driver opens no thread
    !> team: the consumer calls it from inside an element-parallel region, so
    !> threading stays at the element level.
    subroutine writer_emit_batch(w, name, args, outputs, roots, uniform)
        class(leaf_writer_t), intent(inout) :: w
        character(*), intent(in) :: name, args(:), outputs(:), uniform(:)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_ir_t) :: ir
        type(kernel_emit_spec_t) :: spec
        type(str_t) :: source
        logical :: ok
        character(:), allocatable :: message
        integer :: j

        call lower_kernel_ir(roots, ir, ok, message)
        if (.not. ok) error stop message
        spec%name = str(name)
        spec%generator = str(w%generator)
        spec%generator_revision = str(fortsym_revision)
        spec%regenerate_command = str(w%generator//' from pinned FortSym build')
        spec%pure_procedure = .true.
        spec%batch_threaded = .false.
        spec%batch_deterministic = .false.
        allocate (spec%args(size(args)), spec%outputs(size(outputs)), &
                  spec%uniform(size(uniform)))
        do j = 1, size(args)
            spec%args(j) = str(trim(args(j)))
        end do
        do j = 1, size(outputs)
            spec%outputs(j) = str(trim(outputs(j)))
        end do
        do j = 1, size(uniform)
            spec%uniform(j) = str(trim(uniform(j)))
        end do
        source = emit_fortran_batch_ir(ir, spec, ok, message)
        if (.not. ok) error stop message
        write (w%unit, '(a)') chars(source)
    end subroutine

    recursive subroutine emit_derivatives(w, name, args, outputs, roots, &
                                          vjp_interval)
        class(leaf_writer_t), intent(inout) :: w
        character(*), intent(in) :: name, args(:), outputs(:)
        type(expr_t), intent(in) :: roots(:)
        logical, intent(in) :: vjp_interval
        type(expr_t), allocatable :: variables(:), directions(:), weights(:)
        type(expr_t), allocatable :: derivative_roots(:)
        character(32), allocatable :: tangent_args(:), tangent_outputs(:)
        character(32), allocatable :: adjoint_args(:), adjoint_outputs(:)
        integer :: j, n, m

        n = size(args)
        m = size(outputs)
        allocate (variables(n), directions(n), weights(m))
        allocate (tangent_args(2*n), tangent_outputs(m))
        allocate (adjoint_args(n + m), adjoint_outputs(n))
        do j = 1, n
            variables(j) = sym(w%arena, trim(args(j)))
            directions(j) = sym(w%arena, 'd'//trim(args(j)))
            tangent_args(j) = args(j)
            tangent_args(n + j) = 'd'//trim(args(j))
            adjoint_args(j) = args(j)
            adjoint_outputs(j) = 'g'//trim(args(j))
        end do
        do j = 1, m
            weights(j) = sym(w%arena, 'b'//trim(outputs(j)))
            tangent_outputs(j) = 'd'//trim(outputs(j))
            adjoint_args(n + j) = 'b'//trim(outputs(j))
        end do
        derivative_roots = jvp(roots, variables, directions)
        call simplify_all(w, derivative_roots)
        call w%emit(name//'_jvp', tangent_args, tangent_outputs, derivative_roots)
        derivative_roots = vjp(roots, variables, weights)
        call simplify_all(w, derivative_roots)
        call w%emit(name//'_vjp', adjoint_args, adjoint_outputs, derivative_roots)
        if (vjp_interval) call w%emit(name//'_vjp_interval', adjoint_args, &
                                      adjoint_outputs, derivative_roots, &
                                      rigorous=.true.)
    end subroutine

    subroutine simplify_all(w, roots)
        class(leaf_writer_t), intent(inout) :: w
        type(expr_t), intent(inout) :: roots(:)
        type(engine_result_t) :: simplified
        integer :: j
        do j = 1, size(roots)
            simplified = w%engine%simplify(roots(j))
            if (.not. simplified%ok) error stop chars(simplified%message)
            roots(j) = simplified%value
        end do
    end subroutine

    pure logical function flag(value)
        logical, intent(in), optional :: value
        flag = .false.
        if (present(value)) flag = value
    end function
end module kin6d_gen_support
