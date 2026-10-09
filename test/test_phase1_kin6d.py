"""Behavioral test of the Phase 1 producer against the exact E0 physical LCFS."""
from pathlib import Path
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from equilibrium.phase1.producers.kin6d import run

BINARY, WORK = Path(sys.argv[1]).resolve(), Path(sys.argv[2]).resolve()
del sys.argv[1:]


class ProducerTest(unittest.TestCase):
    def test_exact_physical_boundary(self):
        os.environ["KIN6D_BINARY"] = str(BINARY)
        r0, a, f = 6.2, 0.62, 32.86
        c = f / (12 * r0**3)
        mu0 = 4e-7 * np.pi
        theta = np.arange(1024) * 2 * np.pi / 1024
        radius = np.sqrt(r0*r0 + 2*r0*a*np.cos(theta))
        case = {
            "id": "E0_solovev_physical_lcfs", "R0": r0,
            "boundary": {"R": radius, "Z": r0*a*np.sin(theta)/radius},
            "source": {"kind": "solovev", "pprime": -16*c/mu0, "ffprime": 0},
            "psi_edge": 0, "F_edge": f,
        }
        result = run(case, {"boundary_nodes": 96}, WORK)
        self.assertTrue(result.native_converged)
        self.assertGreater(result.dof, 0)
        self.assertGreater(result.wall_s, 0)
        self.assertLess(result.wall_s, 300)
        self.assertLess(np.linalg.norm(np.array(result.axis)-[r0, 0]), 1e-4)
        self.assertLess(result.meta["flux_consistency"], 1e-5)
        # Off-grid signed fields from the exported native FE state.
        rr, zz = np.meshgrid(np.linspace(r0-.4, r0+.4, 9), np.linspace(-.35, .35, 7))
        value = result.evaluate(rr, zz)
        exact_psi = c*((rr*rr-r0*r0)**2 + 4*rr*rr*zz*zz - 4*r0*r0*a*a)
        exact_br = -8*c*rr*zz
        exact_bz = 4*c*(rr*rr-r0*r0 + 2*zz*zz)
        self.assertLess(np.max(np.abs(value["psi"]-exact_psi)), 1e-5)
        self.assertLess(np.max(np.abs(value["BR"]-exact_br)), 3e-4)
        self.assertLess(np.max(np.abs(value["BZ"]-exact_bz)), 3e-4)
        np.testing.assert_allclose(value["Bphi"], f/rr, atol=1e-13, rtol=0)
        outside = result.evaluate(np.array([r0+2, r0]), np.array([0, 2]))
        self.assertTrue(all(np.isnan(v).all() for v in outside.values()))
        # Independent exact toroidal-flux and q map, parameterized by s_pol.
        spol = np.linspace(0, 1, 1025)
        r = np.sqrt(r0*r0 + 2*r0*a*np.sqrt(spol[:, None])*np.cos(theta))
        phi = 2*np.pi*f*r0*r0*a*a*spol*np.mean(np.sin(theta)**2/r**3, axis=1)
        exact_q = 1.5*r0**3*np.mean(r**-3, axis=1)
        targets = np.linspace(.05, .98, 29)
        oracle = np.interp(targets, phi/phi[-1], exact_q)
        self.assertLess(np.max(np.abs(result.q(targets)/oracle-1)), 2e-4)
        self.assertLess(abs(result.phi_edge/phi[-1]-1), 1e-6)
        volume = 2*np.pi**2*np.mean(radius*r0*a*np.cos(theta)
                                  + r0*r0*a*a*np.sin(theta)**2/radius)
        self.assertLess(abs(result.volume/volume-1), 1e-6)
        print(f"producer: DOF={result.dof}, wall_s={result.wall_s:.6f}, "
              f"q max relative error={np.max(np.abs(result.q(targets)/oracle-1)):.3e}")


if __name__ == "__main__":
    unittest.main()
