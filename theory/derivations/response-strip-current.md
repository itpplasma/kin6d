# Physical current and artificial interface sheets

[V] In the periodic-cylinder limit, let a unit-normal stream-function sheet
have amplitude Lambda at r=rs and harmonic m>=1. Choose magnetic scalar
potentials U_in=-Lambda (r/rs)^m/2 and U_out=Lambda (rs/r)^m/2.
Their potential jump is Lambda and their radial magnetic fields are continuous.
With B=mu0 grad U and the response's peak-flux convention Phi=2 R0 r B_r,

    Phi_in  = -mu0 m R0 Lambda (r/rs)^m,
    Phi_out = -mu0 m R0 Lambda (rs/r)^m.

Both satisfy L_m Phi=0 away from rs, where
L_m=(1/r)d_r(r d_r)-m^2/r^2. Integration through the sheet gives

    integral r L_m(Phi) dr = rs [Phi'] = 2 mu0 m^2 R0 Lambda.

`gen/check_response_strip_current.f90` checks these identities exactly in
FortSym for m=1..6, using the native zero test and refusing unknown results.
There is no fitted normalization or external-code input. The physical current
row optionally returned by `strip_dirichlet_map` applies the existing strip
current integral to each of the two Dirichlet solutions. It is linear in the
edge fluxes; its fixed-parameter JVP is itself and its storage VJP is its
Hermitian transpose. Parameter derivatives and enclosures remain unavailable.

[V] The two artificial sheets used by the outer magnetic Schur construction
replace the strip's interior field with a current-free continuation. Write their
radii as a*rs and b*rs, a<1<b. To reproduce the inner and outer magnetic moments
of a sheet Lambda_s at rs, their amplitudes satisfy

    Lambda_in a^m    + Lambda_out b^m    = Lambda_s,
    Lambda_in a^(-m) + Lambda_out b^(-m) = Lambda_s.

For m=2 their unweighted sum is

    (Lambda_in + Lambda_out)/Lambda_s = (a^2 b^2 + 1)/(a^2+b^2).

For a=1-h, b=1+h its difference from one is
h^2(h^2-4)/(2(1+h^2)), tending to -2 h^2. FortSym checks both moment
identities and these two reductions exactly. Thus the interface-amplitude sum
cannot serve as the physical strip-current observable at finite width.

The former R-B2 gate (a) used that sum at
`test/test_mhd_closure.f90`, after `lam_l4` is computed in `fe_strip`.
It therefore included the magnetic replacement's finite-width moment defect.
The revised test also evaluates the actual resolved-strip current using its
edge fluxes and the analytic cylinder normalization above. It retains the
original 1e-3 target at xd=0.02 and the previous half-width and convergence
checks. The R-B2 strip in this test is still single-helicity and zero pressure;
this derivation supplies no finite-pressure toroidal, Glasser or TC24 claim.
