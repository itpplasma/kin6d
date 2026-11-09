# Graph of plasma models

The models of MILCA and their reductions: an arrow runs from a parent model to a
model derived from it. The status of each node is that of its documentation page;
the derivations behind the arrows are planned in the [theory map](theory-map.md).

```mermaid
flowchart TB
    classDef run fill:#ffffff,stroke:#0072B2,stroke-width:3px
    classDef part fill:#ffffff,stroke:#0072B2,stroke-width:2px,stroke-dasharray:6 3
    classDef plan fill:#ffffff,stroke:#D55E00,stroke-width:2px,stroke-dasharray:3 3

    NB["relativistic Maxwell-Lorentz<br/>many-particle system"]:::part
    CL["classical pairs, Coulomb encounters,<br/>one-component plasma, Debye sphere"]:::run
    BBGKY["BBGKY, correlation hierarchy"]:::part
    QM["many-electron quantum mechanics"]:::part
    QF["finite-basis fermions, jellium,<br/>Hartree-Fock, Fermi gas"]:::run
    QW["Wigner, quantum kinetics"]:::part
    FK["Vlasov-Maxwell-Fokker-Planck"]:::part
    LAB["1X1V Vlasov-Poisson-<br/>Fokker-Planck laboratory"]:::run
    COLL["Landau, Rosenbluth,<br/>linearised collision operators"]:::run
    GK["gyrokinetic"]:::plan
    DK["drift kinetic: thin and<br/>finite orbit width"]:::run
    BA["bounce-averaged, resonant"]:::run
    NEO["neoclassical transport:<br/>NTV, bootstrap, E_r"]:::run
    GF["gyrofluid, gyromoments,<br/>Landau fluid"]:::plan
    MF["multi-fluid, Hall, extended MHD"]:::plan
    MHD["visco-resistive MHD,<br/>3D and helical 2D"]:::part
    LIN["linear MHD response,<br/>ideal energy"]:::run
    EQ2["axisymmetric equilibria:<br/>Grad-Shafranov, free boundary, inverse"]:::run
    EQ3["3D equilibria: frozen-in maps,<br/>non-nested, islands, Taylor"]:::run
    EDGE["edge: sheath fluid, neutrals,<br/>full-f drift kinetics"]:::plan
    RX["reacting, partially ionised,<br/>radiation"]:::plan
    TR["macro transport, multiscale coupling"]:::plan

    NB --> CL
    NB --> BBGKY
    QM --> QF
    QM --> QW
    QW --> FK
    BBGKY --> FK
    BBGKY --> COLL
    FK --> LAB
    FK --> GK
    FK --> GF
    COLL --> DK
    GK --> DK
    DK --> BA
    DK --> NEO
    BA --> NEO
    GF --> MF
    MF --> MHD
    MHD --> LIN
    MHD --> EQ2
    MHD --> EQ3
    MF --> EDGE
    DK --> EDGE
    MF --> RX
    NEO --> TR
    MHD --> TR
```

Solid blue: running with documented evidence; dashed blue: partly running;
dotted orange: designed.

| Model | Status | Page |
|---|---|---|
| relativistic many-particle system | partly running | [microscopic many-body](microscopic-many-body.md) |
| classical pairs, Coulomb encounters | running | [classical pair](classical-pair.md), [Coulomb encounter](coulomb-encounter.md) |
| one-component plasma, Debye sphere | running | [OCP](ocp-microscopic.md), [Debye sphere](debye-sphere.md) |
| BBGKY, correlation hierarchy | partly running | [correlation hierarchy](correlation-hierarchy.md) |
| finite-basis fermions, jellium, Hartree-Fock, Fermi gas | running | [quantum Fock](quantum-fock.md), [jellium](jellium-cell.md), [mean field](mean-field-cell.md), [Fermi gas](fermi-gas.md) |
| Wigner, quantum kinetics | partly running | [screening and Wigner](screening-wigner.md), [quantum plasma](quantum-plasma.md) |
| Vlasov-Maxwell-Fokker-Planck | partly running (1X1V) | [first example](first-example.md), [steady laboratory](steady-laboratory.md) |
| collision operators | running | [collisions](collisions.md) |
| gyrokinetic | designed | [model family](model-family.md) |
| drift kinetic, bounce-averaged, neoclassical | running | [drift kinetics](drift-kinetics.md), [extensions](drift-kinetics-extensions.md) |
| gyrofluid, gyromoments, Landau fluid | designed | [moment models](moment-models.md) |
| multi-fluid, Hall, extended MHD | designed | [model family](model-family.md) |
| visco-resistive MHD | partly running | [nonlinear resistive MHD](nonlinear-resistive-mhd.md), [section](vrmhd-section.md) |
| linear MHD response, ideal energy | running | [linear response](mhd-linear-response.md), [ideal energy](ideal-xi.md) |
| axisymmetric equilibria | running | [Grad-Shafranov](grad-shafranov.md), [free boundary](gs-free-boundary.md), [inverse](gs-inverse.md) |
| 3D equilibria | running | [ideal maps](ideal-map-3d.md), [non-nested](equilibrium-nonnested.md), [Taylor](equilibrium-beltrami.md) |
| edge | designed | [edge physics](edge-physics.md) |
| reacting plasmas, radiation | designed | [reacting](reacting-partially-ionized.md), [radiation](radiation.md) |
| macro transport | designed | [transport and multiscale](transport-multiscale.md) |
