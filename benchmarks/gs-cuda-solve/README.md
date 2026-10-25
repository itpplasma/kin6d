# Host-staged FP64 CUDA fixed-source GS prototype

This optional benchmark solves the same host-assembled Grad–Shafranov system
with FortSparse's native CPU and CUDA SPD factors. Geometry, source assembly,
residuals, equilibrated-flux reconstruction and strict audit remain on the CPU.
Factors remain on the GPU, while each solve stages its RHS and answer through
host memory. This is not a resident GPU equilibrium implementation or a
benchmark of the indefinite ideal-map Newton system.

Configure the existing `cpu` preset with
`KIN6D_ENABLE_CUDA_GS_SOLVE=ON`, an explicitly qualified CUDA compiler and
`CMAKE_CUDA_ARCHITECTURES` for the allocated device. Build only
`test_gs_cuda_fixed_source` and the affected native-adapter and layering tests.
The option defaults to OFF. It enables FortSparse's existing native CUDA
backend without changing the default CPU device selection.

The registered tests are `grad_shafranov_cuda_fixed_source_a3` and
`grad_shafranov_cuda_affine_strict_audit`. Run them through CTest with four CPUs
and one admitted GPU. The allocation must set `CUDA_VISIBLE_DEVICES` so that
logical device zero identifies that GPU. No test falls back to a CPU factor.
An explicit GPU request for an indefinite matrix is rejected before dispatch.
Nonlinear profile laws are outside this prototype.

`test_gs_cuda_fixed_source a3 N DEVICE` selects the existing public Solov'ev
test fixture at P2, with `N` boundary nodes from 8 through 64. `affine N DEVICE`
selects a positive-R rectangular constant-profile problem. Every call assembles
one shared matrix and source, factors them once on each backend, and records
five sequential solve measurements. The first solve and GPU setup are retained
separately from subsequent solves. Source and numerical physics leaves are
shared with the CPU path.

The checks use the independent CPU element action rather than the factor's CSR
multiplication. CPU and GPU relative residuals must be at most 1e-9, as must
their relative flux and discrete-energy differences. Equation-level
tangent/adjoint pairing must agree to 1e-10. The CPU-evaluated sampled energy
estimates must agree to 1e-7 relative. The A3 exact energy error is a separate
high-order quadrature diagnostic that includes the represented-boundary defect.
The supported affine case must also pass the existing outward-arithmetic
continuous energy audit. The curved A3 audit must report unavailable, and its
sampled estimate is not promoted to a certificate.

Standard output retains raw solve times, setup/factor wall times, separate CPU
residual, estimator and audit times, residual/error diagnostics, and the active
device and explicit device/host
payload reported by FortSparse. GPU setup includes transfers and runtime
initialization. Solve times include both RHS and answer transfers. Separate
transfer timing is unavailable, and the reported transfer volume is the
FP64 RHS plus answer, 16 bytes per degree of freedom per scalar solve.
The 1 GiB device limit covers explicit backend-owned payload, not driver,
vendor-library or process RSS. The native operator must admit actual host/GPU
headroom separately.

Retain stdout, effective inputs, exact source and dependency revisions,
compiler/flags, hardware/driver, device visibility and resource receipts outside
Git. The `*_first_checked_result_s` fields are aggregate sums of recorded
host assembly, setup/factor, first solve, CPU residual, estimator and audit
attempt phases at the same claim status. They are not stopwatch measurements
of first-result latency: checking occurs after the CPU-first repeated solves.
Exact-oracle and cross-backend validation costs are reported separately.
A GPU solve is admitted only when the reported active device and payload
prove that dispatch occurred. Compare setup-inclusive and repeated-solve cost
at the same mesh, error and estimator status. No speed advantage is asserted
without those measurements.

The explicit GPU route constructs CSC from the shared CSR with FortSparse
`csc_from_triplet`, including explicit row indices, because
the pinned budgeted raw-conversion API admits CPU execution only. Its copy time
is included in GPU setup wall time, its bytes are reported separately, and the
1 GiB explicit backend device-payload limit remains active.
Constructor workspace is included in whole-process RSS measurements and setup
wall time; persistent CSC and explicit triplet-row payload are also reported.

## Public fixed-source solver entry

`initialize_gs` also accepts an explicit `device` and optional positive
`device_memory_limit_bytes`. Device -1 is the unchanged CPU default; a
nonnegative device requests the native SPD accelerator factor, with a default
1 GiB explicit payload cap. `gs_requested_device` reports that policy, not
proof of a successful device factor. The retained factor statistics report the
actual active device and payload.

The public `solve_gs`, `gs_linear_response` and `gs_apply_inverse` entries use
that cached factor. Its device policy is preserved by state copies and by a
factor reset. An explicit device failure or a failed strict linear residual
returns `GS_DEVICE_FAILURE`; it never falls back to CPU CG. CPU default and
the existing CPU iterative selection retain their previous behavior.

`grad_shafranov_cuda_public_fixed_source` tests the public A3 solve, independent
element residual, CPU/GPU flux and energy agreement, implicit inverse pairing,
reuse with a different source, and reset. Its first-solve times are stopwatch
measurements of the public calls, including factor setup and residual checks;
mesh initialization is reported separately. The bound status remains that of
the original fixed-source case, without a new continuous certificate.

This selection does not implement GPU free-boundary or cubic-profile Newton
solves. Those entries explicitly return `GS_FB_UNSUPPORTED_DEVICE` or
`GS_UNSUPPORTED_DEVICE` rather than running their indefinite or updated
operators on CPU under a GPU label. The public FreeGS parity test remains a
CPU regression with its original profiles, current constraints and tolerances.
