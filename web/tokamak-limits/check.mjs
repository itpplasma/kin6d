// Runs the WebAssembly build on the ITER-like case of test_tokamak_discharge
// and prints the same quantities, so the LFortran/emscripten build can be
// compared with the native gfortran build: node check.mjs OUTDIR/milca_tokamak.mjs
import {pathToFileURL} from "node:url";

const mod = await (await import(pathToFileURL(process.argv[2]).href)).default();
const ints = mod._malloc(6 * 4);
mod._milca_td_sizes(ints, ints + 4, ints + 8, ints + 12, ints + 16, ints + 20);
const [nx, npar, nout, nser, nprof, nramp] = Array.from(mod.HEAP32.subarray(ints / 4, ints / 4 + 6));
const f64 = (n) => mod._malloc(8 * n);
const par = f64(npar), out = f64(nout), prof = f64(nx * nprof), st = mod._malloc(4);
const P = [6.2, 2.0, 1.75, 0.33, 5.3, 15.0, 1.0, 1.0, 0.3, 20.0, 1.0, 1.7, 0.005, 1.3, 0.3, 1.0, 2.5, 0.1];
mod.HEAPF64.fill(0, par / 8, par / 8 + npar);
mod.HEAPF64.set(P, par / 8);
let t = performance.now();
mod._milca_td_limits(par, out, prof, st);
const o = mod.HEAPF64.subarray(out / 8, out / 8 + nout);
console.log(`limits ${(performance.now() - t).toFixed(0)} ms status ${mod.HEAP32[st / 4]}`);
console.log(`q0 q95 ${o[0].toFixed(3)} ${o[1].toFixed(3)} li ${o[2].toFixed(3)} betaN ${o[4].toFixed(3)} fGW ${o[7].toFixed(3)}`);
console.log(`2/1 r_s Delta' ${o[29].toFixed(4)} 3/2 ${o[31].toFixed(4)} ballooning ratio ${o[32].toFixed(3)} kink dW m=4 ${o[19].toExponential(3)}`);
const maxrec = 6000, ramp = f64(nramp), nrec = mod._malloc(4);
const ser = f64(nser * maxrec), pro = f64(nx * nprof * maxrec);
const scenarios = {"density ramp": [4, 2.4, -1, 15, 8, 1], "power ramp": [3, 1, 400, 15, 6, 1],
  "current ramp": [3, 1, -1, 24, 6, 1], "safe ramp": [2, 1.05, -1, 15, 4, 0.5]};
for (const [name, r] of Object.entries(scenarios)) {
  mod.HEAPF64.fill(0, ramp / 8, ramp / 8 + nramp);
  mod.HEAPF64.set(r, ramp / 8);
  t = performance.now();
  mod._milca_td_run(par, ramp, maxrec, nrec, ser, pro, st);
  const n = mod.HEAP32[nrec / 4], s = (k, i) => mod.HEAPF64[ser / 8 + i * nser + k - 1];
  let kq = -1;
  for (let i = 0; i < n; i++) if (s(21, i) >= 3 && s(21, i) <= 5) { kq = i; break; }
  console.log(`${name}: records ${n} t_end ${s(1, n - 1).toFixed(5)} phase ${s(21, n - 1)} cause ${s(22, n - 1)} (${(performance.now() - t).toFixed(0)} ms)`);
  if (kq >= 0) console.log(`  TQ at t= ${s(1, kq).toFixed(5)} fGW ${s(6, kq).toFixed(3)} betaN ${s(7, kq).toFixed(3)} q95 ${s(9, kq).toFixed(2)}`);
}
