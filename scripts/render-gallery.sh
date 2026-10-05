#!/usr/bin/env bash
# Render already-built Fortran demos, preserving numerical evidence alongside media.
# This script never configures, builds, tests, uploads, or changes a Fo session.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="$repo_root/build"
gallery_dir="$repo_root/output/gallery-$(date -u +%Y%m%dT%H%M%SZ)"
snapshot=false
available=false
selection=all
verification_receipt=

usage() {
    cat <<'HELP'
Usage: scripts/render-gallery.sh [options]
  --build-dir DIR             Existing validated executable directory (default build)
  --output-dir DIR            New gallery directory (default output/gallery-UTC)
  --only TARGET[,TARGET...]   Select registered demos (default all)
  --snapshot                 PNG/data only; do not encode movies or GIFs
  --available                Record and skip executables that are not built yet
  --verification-receipt FILE Copy a controller's numerical verification receipt
  --list                     List executable targets and primary movie names
The script runs demos sequentially from the repository directory. GIFs use
12 fps, width 960 and a shared palette. MP4s and data remain unchanged.
HELP
}

# target | primary movie | model | bound scope | snapshot argument supported
entries=(
    'interaction_demo|collective-response|Neutral elastic contact and unscreened Coulomb comparison|Numerical/physical certification scope is recorded by the application evidence; no self-consistent collective plasma claim|yes'
    'streaming_moments_demo|particles-to-moments|Zero-field free streaming on v[-1,1], exact finite Fourier/P1 initial data|Continuous normalized density/current and bulk velocity bounds when the enclosed density is strictly positive; no collisions, fields or sensitivity certificate|yes'
    'debye_sphere_demo|debye-potential-reduction|Linearized spherical Debye-Huckel finite permeable source equilibrium|Strict equilibrium reference intervals, regular-center series remainder and radial chord margins; parameter sweep is not a transient|yes'
    'shielding_dynamics_demo|debye-shielding|Periodic planar linear drift-diffusion-Poisson, overdamped electrons|Continuous finite-Fourier residual/contractivity bounds for electron excess, potential and electric field; no microscopic trajectory or Vlasov claim|yes'
    'slab_oscillation_demo|plasma-oscillation|Numerically integrated cold electron slab against fixed positive ions|Independent residual certificate for the FortNum Cash-Karp oscillator trajectory; analytic cos/sin overlay is an oracle; no thermal, collisional or kinetic damping claim|yes'
    'microscopic_cloud_demo|microscopic-cloud|Isolated three-dimensional 32-electron Coulomb cloud in a fixed positive sphere, with and without a central impurity|Exploratory candidate trajectories and empirical profiles; no trajectory, thermodynamic Debye shielding or finite-N fluctuation certificate|yes'
    'streaming_comparison_demo|particles-vlasov|Periodic prescribed-field-free 1X1V streaming, weighted characteristics and Eulerian reconstruction|Continuous finite Fourier/P1 distribution and normalized density enclosures; deterministic particle Fourier quadrature bound; no weak sensitivity certificate|yes'
    'kin6d_demo|free-streaming|Prescribed-field 1X1V reference and free streaming|Continuous pointwise reference enclosure; renderer does not rerun a formal proof|yes'
    'quantum_reference_demo|quantum-correlations|Interacting opposite-spin electron pair in a periodic 3D neutralized Coulomb box, finite plane-wave basis|Continuous-time finite-basis Hilbert residual and observable/interpolation bounds; continuum projection error unavailable|yes'
    'classical_pair_demo|classical-coulomb-pair|Classical ensemble of true 3D two-electron trajectories in the same periodic neutralized Coulomb box|Simultaneous sampling confidence intervals for numerical bin probabilities; force/time-step errors are convergence diagnostics, not deterministic trajectory bounds|yes'
    'steady_laboratory_demo||Bounded OU/FP and manufactured coupled VPFP reference laboratory|Reviewed continuous reference-profile interval bounds; no general nonlinear uniqueness or coupled formal theorem|no'
    'transport_demo||Arbitrary finite Fourier/P1 damped transport|Continuous residual/stability state and directional tangent bounds, connected plot-cell bands and limiting-case comparison|no'
)

while (($#)); do
    case "$1" in
        --build-dir|--output-dir|--only|--verification-receipt)
            (($# >= 2)) || { usage >&2; exit 2; }
            case "$1" in
                --build-dir) build_dir="$2" ;;
                --output-dir) gallery_dir="$2" ;;
                --only) selection="$2" ;;
                --verification-receipt) verification_receipt="$2" ;;
            esac
            shift 2 ;;
        --snapshot) snapshot=true; shift ;;
        --available) available=true; shift ;;
        --list)
            for entry in "${entries[@]}"; do
                IFS='|' read -r target movie _ <<<"$entry"
                printf '%s\t%s\n' "$target" "${movie:-static plots}"
            done
            exit 0 ;;
        -h|--help) usage; exit 0 ;;
        *) usage >&2; exit 2 ;;
    esac
done

cd "$repo_root"
for command in jq sha256sum git; do
    command -v "$command" >/dev/null || { echo "Missing command: $command" >&2; exit 2; }
done
if ! "$snapshot"; then
    command -v ffmpeg >/dev/null || { echo 'Missing command: ffmpeg' >&2; exit 2; }
    command -v ffprobe >/dev/null || { echo 'Missing command: ffprobe' >&2; exit 2; }
fi
build_dir="$(realpath "$build_dir")"
# GNU time reports resources with -v; BSD/macOS time uses -l.
time_flag=-v
[[ "$(uname -s)" != Darwin ]] || time_flag=-l
mkdir -p "$gallery_dir"
gallery_dir="$(realpath "$gallery_dir")"
[[ "$gallery_dir" != "$repo_root/output" ]] || { echo 'Choose an output subdirectory' >&2; exit 2; }
if [[ -e "$gallery_dir" ]]; then
    [[ -d "$gallery_dir" ]] || { echo 'Gallery path is not a directory' >&2; exit 2; }
    [[ -z "$(find "$gallery_dir" -mindepth 1 -maxdepth 1 -print -quit)" ]] || {
        echo 'Choose a new or empty gallery directory' >&2; exit 2;
    }
fi
if [[ "$selection" != all ]]; then
    IFS=',' read -ra requested <<<"$selection"
    for item in "${requested[@]}"; do
        found=false
        for entry in "${entries[@]}"; do
            [[ "${entry%%|*}" == "$item" ]] && found=true
        done
        "$found" || { echo "Unknown gallery target: $item" >&2; exit 2; }
    done
fi
mkdir -p "$gallery_dir/evidence"
records="$gallery_dir/evidence/artifacts.jsonl"
runs="$gallery_dir/evidence/runs.jsonl"
: > "$records"
: > "$runs"
git rev-parse HEAD > "$gallery_dir/evidence/source-base-commit.txt"
git diff --binary HEAD | sha256sum > "$gallery_dir/evidence/source-worktree-patch.sha256"
git ls-files -co --exclude-standard -z -- '*.f90' 'CMakeLists.txt' \
    'CMakePresets.json' 'cmake/*' 'test/*.cmake' 'docs/*' 'scripts/*' | sort -zu | xargs -0 sha256sum \
    > "$gallery_dir/evidence/source-manifest.sha256"
if [[ -n "$verification_receipt" ]]; then
    cp -- "$verification_receipt" "$gallery_dir/evidence/numerical-verification-receipt.txt"
fi
if ! "$snapshot"; then
    ffmpeg -version > "$gallery_dir/evidence/ffmpeg-version.txt"
    ffprobe -version > "$gallery_dir/evidence/ffprobe-version.txt"
fi
started="$(date -u +%FT%TZ)"
finalize() {
    result=$?
    jq -n --arg started "$started" --arg completed "$(date -u +%FT%TZ)" \
        --arg source "$(cat "$gallery_dir/evidence/source-base-commit.txt")" \
        --arg build "$build_dir" --argjson exitcode "$result" \
        --slurpfile artifacts "$records" --slurpfile runs "$runs" \
        '{schema:"kin6d-render-gallery-v1",started_utc:$started,completed_utc:$completed,
          source_snapshot_at_render:$source,build_directory:$build,exitcode:$exitcode,
          numerical_verification:"No tests run by renderer; see supplied receipt if present",
          display_scope:"Encoded pixels are a visualization; NetCDF-4 results (kind netcdf_result, with embedded configuration/provenance) and CSV/interval data retain the numerical bounds",
          runs:$runs,artifacts:$artifacts}' > "$gallery_dir/manifest.json"
    echo "Gallery: $gallery_dir (exit $result)"
}
trap finalize EXIT

record_artifact() {
    local path="$1" kind="$2" elapsed="$3"
    jq -nc --arg path "${path#"$gallery_dir/"}" --arg target "$target" \
        --arg model "$model" --arg bounds "$bounds" --arg kind "$kind" \
        --arg sha "$(sha256sum "$path" | cut -d' ' -f1)" --arg duration "$elapsed" \
        '{path:$path,target:$target,model:$model,bounds_scope:$bounds,
          kind:$kind,sha256:$sha,render_or_conversion_wall_seconds:$duration}' >> "$records"
}

wall_seconds() {
    awk -v start="$1" -v finish="$2" 'BEGIN {printf "%.6f",finish-start}'
}

encode_movie() {
    local movie_path="$1" prefix metadata gif_start gif_end elapsed frames middle last
    prefix="${movie_path%.mp4}"
    metadata="$prefix.ffprobe.json"
    ffprobe -v error -select_streams v:0 -show_entries \
        stream=codec_name,width,height,nb_frames,r_frame_rate:format=duration \
        -of json "$movie_path" > "$metadata"
    jq -e '.streams[0].width > 0 and .streams[0].height > 0 and
        (.format.duration|tonumber) > 0' "$metadata" >/dev/null
    frames="$(jq -r '.streams[0].nb_frames // "N/A"' "$metadata")"
    if [[ ! "$frames" =~ ^[0-9]+$ ]]; then
        frames="$(ffprobe -v error -count_frames -select_streams v:0 \
            -show_entries stream=nb_read_frames -of default=nw=1:nk=1 "$movie_path")"
    fi
    [[ "$frames" =~ ^[0-9]+$ ]] && ((frames > 0))
    middle=$((frames/2)); last=$((frames-1))
    gif_start="$(date +%s.%N)"
    ffmpeg -hide_banner -loglevel error -threads 1 -filter_complex_threads 1 \
        -y -i "$movie_path" -filter_complex \
        'fps=12,scale=960:-1:flags=lanczos,split[a][b];[a]palettegen=stats_mode=diff[p];[b][p]paletteuse=dither=sierra2_4a' \
        -loop 0 "$prefix.gif"
    gif_end="$(date +%s.%N)"
    elapsed="$(wall_seconds "$gif_start" "$gif_end")"
    ffprobe -v error -show_entries stream=width,height,nb_frames -of json \
        "$prefix.gif" > "$prefix.gif.ffprobe.json"
    jq -e '.streams[0].width == 960 and (.streams[0].nb_frames|tonumber) > 0' \
        "$prefix.gif.ffprobe.json" >/dev/null
    record_artifact "$prefix.gif" animated_gif "$elapsed"
    # Decode exact first/middle/last frame indices, retaining native MP4 dimensions.
    for position in first middle last; do
        case "$position" in first) index=0 ;; middle) index="$middle" ;; last) index="$last" ;; esac
        ffmpeg -hide_banner -loglevel error -threads 1 -y -i "$movie_path" \
            -vf "select=eq(n\,$index)" -frames:v 1 "$prefix-$position.png"
        [[ -s "$prefix-$position.png" ]]
        record_artifact "$prefix-$position.png" movie_frame_snapshot ''
    done
    record_artifact "$metadata" video_metadata ''
    record_artifact "$prefix.gif.ffprobe.json" gif_metadata ''
}

for entry in "${entries[@]}"; do
    IFS='|' read -r target movie model bounds supports_snapshot <<<"$entry"
    [[ "$selection" == all || ",$selection," == *",$target,"* ]] || continue
    binary="$build_dir/$target"
    if [[ ! -x "$binary" ]]; then
        jq -nc --arg target "$target" '{target:$target,status:"missing_executable"}' >> "$runs"
        if "$available"; then echo "Skipped unavailable target: $target"; continue; fi
        echo "Missing executable: $binary" >&2
        exit 1
    fi
    destination="$gallery_dir/$target"
    mkdir -p "$destination"
    marker="$destination/render-start.marker"
    touch "$marker"
    binary_sha="$(sha256sum "$binary" | cut -d' ' -f1)"
    args=()
    if "$snapshot" && [[ "$supports_snapshot" == yes ]]; then args=(--snapshot); fi
    echo "Rendering $target"
    run_start="$(date +%s.%N)"
    if /usr/bin/time "$time_flag" -o "$destination/resources.txt" "$binary" ${args[@]+"${args[@]}"} \
        > "$destination/run.log" 2>&1; then
        run_exit=0
    else
        run_exit=$?
    fi
    run_end="$(date +%s.%N)"
    elapsed="$(wall_seconds "$run_start" "$run_end")"
    jq -nc --arg target "$target" --arg binary "$binary" --arg sha "$binary_sha" \
        --arg start "$run_start" --arg end "$run_end" --arg wall "$elapsed" \
        --argjson exitcode "$run_exit" --argjson args "$(printf '%s\n' ${args[@]+"${args[@]}"} | jq -Rsc 'split("\n")|map(select(length>0))')" \
        '{target:$target,executable:$binary,sha256:$sha,args:$args,
          start_epoch_seconds:$start,end_epoch_seconds:$end,
          application_wall_seconds:$wall,exitcode:$exitcode}' >> "$runs"
    ((run_exit == 0)) || { echo "Failed $target; see $destination/run.log" >&2; exit "$run_exit"; }
    png_count=0
    while IFS= read -r -d '' artifact; do
        name="${artifact##*/}"
        cp -- "$artifact" "$destination/$name"
        case "$name" in
            *.png) kind=snapshot; png_count=$((png_count+1)) ;;
            *.mp4) kind=source_movie ;;
            *.nc) kind=netcdf_result ;;
            *) kind=numerical_evidence ;;
        esac
        record_artifact "$destination/$name" "$kind" "$elapsed"
    done < <(find output -maxdepth 1 -type f -newer "$marker" -print0 | sort -z)
    ((png_count > 0)) || { echo "No fresh snapshot from $target" >&2; exit 1; }
    # F5 contract: every application writes output/<target>.nc with its
    # literal/resolved configuration and provenance; sha256 is recorded above.
    [[ -s "$destination/$target.nc" ]] || { echo "No fresh $target.nc result" >&2; exit 1; }
    if ! "$snapshot" && [[ -n "$movie" ]]; then
        [[ -s "$destination/$movie.mp4" ]] || { echo "No fresh primary movie from $target" >&2; exit 1; }
        for clip in "$destination"/*.mp4; do encode_movie "$clip"; done
    fi
    record_artifact "$destination/run.log" application_log ''
    record_artifact "$destination/resources.txt" process_resources ''
done
