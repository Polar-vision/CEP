# Synthetic Parallax Short-Baseline Test

This test checks whether the CEP parallax point parameterization can recover reliable 3D structure when image geometry has very weak baseline.

The test is intentionally small and controlled:

- Camera poses are fixed to the synthetic ground truth.
- Only landmark parameters are optimized.
- Three parameterizations are compared on exactly the same observations and initial points:
  - `xyz`
  - `inverse_depth`
  - `parallax`
- Image observations include Gaussian pixel noise.

This isolates the point parameterization behavior from full bundle-adjustment gauge freedom.

## Build

Build the test from the existing `example_v2` project:

```bat
cmake --build C:\zuo\Projects\repos\CEP\example_v2\build --target synthetic_parallax_test -j
```

If `example_v2` has not been configured yet, configure it first using the commands in the main README.

## Run

```bat
C:\zuo\Projects\repos\CEP\example_v2\build\synthetic_parallax_test.exe
```

The executable writes synthetic CEP-format input files and a CSV summary to:

```text
C:\zuo\Projects\repos\CEP\example_v2\build\synthetic_data
```

The generated scenes are:

- `normal_baseline`: 5 cameras, 1.0 m step, about 4.37 deg endpoint parallax
- `short_baseline`: 5 cameras, 0.02 m step, about 0.088 deg endpoint parallax
- `ultra_short_baseline`: 5 cameras, 0.005 m step, about 0.022 deg endpoint parallax

Each scene contains:

- `cal.txt`
- `Cam.txt`
- `Feature.txt`
- `XYZ.txt`
- `XYZ_truth.txt`

## Metrics

The test prints and writes:

- `initial_rms_px`: initial reprojection RMS in pixels
- `final_rms_px`: final reprojection RMS in pixels
- `point_rmse_m`: 3D point RMSE against synthetic truth
- `median_depth_relative_error`: median relative depth error
- `iterations`
- `termination`

## Current Result

With 0.3 px image noise, the observed behavior is:

```text
Scene: normal_baseline
  xyz / inverse_depth / parallax all converge to about 0.254 px final RMS
  point RMSE is about 0.35 m
  median depth relative error is about 0.4%

Scene: short_baseline
  all three methods still reach about 0.254 px final RMS
  point RMSE grows to about 38 m
  median depth relative error grows to about 20%

Scene: ultra_short_baseline
  all three methods still reach about 0.254 px final RMS
  inverse_depth and parallax point RMSE are about 470 m
  xyz becomes far less stable in 3D
  median depth relative error is about 55%
```

## Interpretation

Parallax parameterization improves numerical behavior compared with raw XYZ in extremely weak geometry, especially by avoiding the huge XYZ blow-up seen in the ultra-short-baseline case. However, it does not truly solve the short-baseline depth observability problem.

When endpoint parallax is around 0.02 to 0.09 degrees and image noise is 0.3 px, many very different depths explain the image observations almost equally well. The final reprojection RMS can be excellent while the recovered 3D depth is still poor. This is an information/geometry limit, not just a parameterization issue.

For real BA experiments, parallax is useful as a more stable parameterization for low-parallax tracks, but reliable depth still needs stronger geometry, lower noise, extra priors, known scale/depth constraints, or additional observations with wider baseline.
