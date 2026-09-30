// isentropic_interp.cpp
//
// Linear interpolation (in theta) of model-level fields onto a fixed set of
// target isentropic levels. Native backend for enm.interpolation.isentropic.
//
// Bracketing rule (matches the IDL routine):
//   For each target theta_t, take the LAST model level k (counting from the
//   surface upwards) with theta[k] < theta_t. Every level above k then has
//   theta >= theta_t. This is the same as IDL's  where(theta lt thm)  applied
//   in top->ground order, taking the first hit, i.e. the TOPMOST crossing.
//   The result is bracketed by (k, k+1), with theta[k] < theta_t <= theta[k+1],
//   so the denominator is never zero.
//   Columns need NOT be monotonic (statically unstable layers are fine).
//   No bracket (theta_t below every level, or above the top level) -> fill_value.
//   Extrapolation above the top level, which the IDL does, is NOT done.
//   Lower-boundary masking (theta_t < theta_lb) is applied in Python.
//
// Preconditions (enforced by the Python wrapper):
//   - th is (nlev, ncols), row-major, level axis ordered surface -> top
//     (index 0 = lowest model level).
//   - vars is (nvars, nlev, ncols) with the same level ordering.
//   - thlev is monotonically increasing, length ntarget.
//
// Algorithm: the last index with th[l] < t equals the largest l such that
// sufmin[l] < t, where sufmin[l] = min(th[l..nlev-1]). sufmin is
// non-decreasing in l and the targets are increasing, so a forward-only
// merge pointer finds every bracket in O(nlev + ntarget) per column.

#include <cstddef>
#include <vector>

extern "C" {

void enm_interp_theta_levels(
    const float* th,      // (nlev, ncols)
    const float* vars,    // (nvars, nlev, ncols)
    long nlev,
    long ncols,
    long nvars,
    const float* thlev,   // (ntarget,)
    long ntarget,
    float* out,           // (nvars, ntarget, ncols)
    float fill_value
) {
    if (nlev < 2) {
        for (long v = 0; v < nvars; ++v)
            for (long k = 0; k < ntarget; ++k)
                for (long c = 0; c < ncols; ++c)
                    out[(v * ntarget + k) * ncols + c] = fill_value;
        return;
    }

    #pragma omp parallel
    {
        std::vector<float> sufmin(static_cast<size_t>(nlev));

        #pragma omp for schedule(static)
        for (long c = 0; c < ncols; ++c) {
            const float* th_col = th + c;   // th_col[l * ncols] == th(l, c)

            sufmin[nlev - 1] = th_col[(nlev - 1) * ncols];
            for (long l = nlev - 2; l >= 0; --l) {
                const float x = th_col[l * ncols];
                sufmin[l] = (x < sufmin[l + 1]) ? x : sufmin[l + 1];
            }

            long lev = -1;   // last level with th < target; only moves forward

            for (long k = 0; k < ntarget; ++k) {
                const float target = thlev[k];

                while (lev + 1 < nlev && sufmin[lev + 1] < target) {
                    ++lev;
                }

                // lev == -1: target at/below every level (below ground)
                // lev == nlev-1: target above the top level
                const bool valid = (lev >= 0) && (lev < nlev - 1);

                for (long v = 0; v < nvars; ++v) {
                    float* out_elem = out + ((v * ntarget + k) * ncols) + c;
                    if (!valid) {
                        *out_elem = fill_value;
                        continue;
                    }
                    const float* var_col = vars + (v * nlev * ncols) + c;
                    const float th0 = th_col[lev * ncols];
                    const float th1 = th_col[(lev + 1) * ncols];
                    const float f0 = var_col[lev * ncols];
                    const float f1 = var_col[(lev + 1) * ncols];
                    const float w = (th1 > th0) ? (target - th0) / (th1 - th0) : 0.0f;
                    *out_elem = f0 + w * (f1 - f0);
                }
            }
        }
    }
}

} // extern "C"


// // isentropic_interp.cpp
// //
// // Linear interpolation of model-level fields onto a fixed set of target
// // isentropic (theta) levels.
// //
// // This is the native backend for enm.interpolation.isentropic; the pure
// // NumPy version there is the fallback used when this extension hasn't
// // been built. Called via ctypes -- plain extern "C" function, raw
// // pointers, no STL types crossing the boundary.
// //
// // Preconditions (enforced by the Python wrapper, not here):
// //   - `th` is laid out as (nlev, ncols), row-major, and for every column
// //     is monotonically INCREASING along the level axis. The wrapper
// //     reorders the native (possibly decreasing) level axis before calling
// //     in, since the direction is a global modeling convention, not
// //     something that varies column to column.
// //   - `vars` is laid out as (nvars, nlev, ncols), using that same level
// //     ordering as `th`.
// //   - `thlev` is monotonically increasing, length `ntarget`.
// //
// // Algorithm: for each column, `th` (length nlev) and `thlev` (length
// // ntarget) are both sorted, so the bracketing model level for every
// // target can be found with a single linear merge pass over the two
// // sequences (the position pointer into `th` only ever moves forward as
// // we sweep through increasing target levels) -- O(nlev + ntarget) per
// // column, rather than re-searching from scratch for every target level.
// // Columns are independent, so the outer loop is parallelized with OpenMP.

// #include <cstddef>

// extern "C" {

// void enm_interp_theta_levels(
//     const float* th,      // (nlev, ncols)
//     const float* vars,    // (nvars, nlev, ncols)
//     long nlev,
//     long ncols,
//     long nvars,
//     const float* thlev,   // (ntarget,)
//     long ntarget,
//     float* out,            // (nvars, ntarget, ncols)
//     float fill_value
// ) {
//     #pragma omp parallel for schedule(static)
//     for (long c = 0; c < ncols; ++c) {
//         const float* th_col = th + c;   // th_col[l * ncols] == th(l, c)
//         long lev = 0;                   // merge pointer, only moves forward

//         for (long k = 0; k < ntarget; ++k) {
//             const float target = thlev[k];

//             // Advance while the level above `lev` is still <= target.
//             while (lev < nlev - 1 && th_col[(lev + 1) * ncols] <= target) {
//                 ++lev;
//             }

//             const bool valid =
//                 (nlev >= 2) &&
//                 (lev < nlev - 1) &&
//                 (th_col[lev * ncols] <= target) &&
//                 (target <= th_col[(lev + 1) * ncols]);

//             for (long v = 0; v < nvars; ++v) {
//                 float* out_elem = out + ((v * ntarget + k) * ncols) + c;
//                 if (!valid) {
//                     *out_elem = fill_value;
//                     continue;
//                 }
//                 const float* var_col = vars + (v * nlev * ncols) + c;
//                 const float th0 = th_col[lev * ncols];
//                 const float th1 = th_col[(lev + 1) * ncols];
//                 const float f0 = var_col[lev * ncols];
//                 const float f1 = var_col[(lev + 1) * ncols];
//                 const float w = (th1 > th0) ? (target - th0) / (th1 - th0) : 0.0f;
//                 *out_elem = f0 + w * (f1 - f0);
//             }
//         }
//     }
// }

// } // extern "C"