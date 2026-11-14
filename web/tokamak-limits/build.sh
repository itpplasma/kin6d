#!/bin/sh
# Builds the teaching discharge model (milca_tokamak_discharge with the
# cylindrical stability and s-alpha ballooning leaves) for the browser:
#   build.sh OUTDIR
# Writes OUTDIR/milca_tokamak.mjs, an ES module with the WebAssembly embedded
# (emscripten SINGLE_FILE), exporting the factory of the emscripten module with
# _milca_td_sizes, _milca_td_limits and _milca_td_run.
# Needs LFortran built with WITH_TARGET_WASM (LFORTRAN, default lfortran) and an
# emscripten SDK in EMSDK_PATH.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
top=$(cd "$here/../.." && pwd)
out=$(mkdir -p "${1:-$top/build/web-tokamak}" && cd "${1:-$top/build/web-tokamak}" && pwd)
emcc=${EMCC:-${EMSDK_PATH:?set EMSDK_PATH to the emscripten SDK}/upstream/emscripten/emcc}
lf=${LFORTRAN:-lfortran}
runtime=$(dirname "$(command -v "$lf")")/../lib/lfortran_runtime_wasm_emcc.o
flags="--target=wasm32-unknown-emscripten --separate-compilation --realloc-lhs-arrays -O3"
rm -rf "$out/obj"
mkdir -p "$out/obj"
cd "$out/obj"
for f in mhd/ballooning_s_alpha mhd/cylinder_stability continuum/tokamak_discharge; do
    "$lf" -c $flags "$top/src/$f.f90" -o "$(basename "$f").o"
done
"$emcc" --target=wasm32-unknown-emscripten -O2 -sALLOW_MEMORY_GROWTH \
    -sSTACK_SIZE=8mb -sMODULARIZE -sEXPORT_ES6 -sSINGLE_FILE -sENVIRONMENT=web,node \
    -sEXPORT_NAME=createMilcaTokamak \
    -sEXPORTED_FUNCTIONS=_malloc,_free,_milca_td_sizes,_milca_td_limits,_milca_td_run \
    -sEXPORTED_RUNTIME_METHODS=HEAPF64,HEAP32 \
    -o "$out/milca_tokamak.mjs" ./*.o "$runtime"
rm -rf "$out/obj"
echo "$out/milca_tokamak.mjs"
