# Lessons learned

Each campaign that brings a physics area of MILCA to agreement with a
reference leaves one page here. It lists every mismatch that was found, its
cause, how it was found, and the check that now catches it. Read the page of
the physics area before you start work on it, and add to it when a new cause
appears.

| Page | Physics area | Campaign |
|---|---|---|
| [Drift-kinetic transport](drift-kinetic-transport.md) | NTV torque, axisymmetric transport matrix, finite orbit width | agreement with NEO-2-QL, NEO-RT and POTATO on ITER TC24 and ASDEX Upgrade 30835, November 2026 |
| [Global ideal-MHD stability](ideal-mhd-stability.md) | delta W, hybrid radial spaces, spectra | P1e, from November 2026 |

## Page structure

A page opens with the physics area, the references and the outcome in a few
sentences. The causes follow, grouped by kind:

1. wrong formula (a term derived or transcribed wrongly, including signs);
2. missing term;
3. different model or convention in the reference;
4. resolution of MILCA or of the reference;
5. wrong limit or reference in a check;
6. inconsistent input data;
7. implementation and toolchain defects.

Each entry says what happened and its size, how it was found, and the check
or rule that prevents it now, with the pull request that made the change.
The page closes with the rules that follow from the campaign.
