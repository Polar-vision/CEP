#!/usr/bin/env python3
import argparse
import csv
import math
import random
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path


METHOD_ORDER = [
    "world_cartesian_xyz",
    "world_xy_inverse_z",
    "world_bearing_range",
    "world_bearing_inverse_range",
    "anchored_cartesian_xyz",
    "anchored_xy_inverse_z",
    "anchored_bearing_range",
    "anchored_bearing_inverse_range",
    "dual_anchor_bearing_parallax",
]

NOISE_LEVELS = [
    ("clean", 0.0, 0.0, 0.0, [0]),
    ("tiny", 0.02, 0.01, 0.02, [0, 1, 2]),
    ("small", 0.05, 0.03, 0.05, [0, 1, 2]),
    ("medium", 0.20, 0.10, 0.20, [0, 1, 2]),
    ("large", 0.50, 0.30, 0.50, [0, 1, 2]),
    ("severe", 1.00, 0.70, 1.00, [0, 1, 2]),
]


def read_numeric_rows(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        rows.append([float(x) for x in line.split()])
    return rows


def write_cam(path, rows):
    with Path(path).open("w", newline="\n") as f:
        for row in rows:
            values = [f"{v:.9f}" for v in row[:6]]
            if len(row) > 6:
                values.append(str(int(round(row[6]))))
            f.write(" ".join(values) + "\n")


def write_xyz(path, rows):
    with Path(path).open("w", newline="\n") as f:
        for row in rows:
            f.write(" ".join(f"{v:.9f}" for v in row[:3]) + "\n")


def euler_to_rotation(e):
    ey, ex, ez = e[:3]
    c1, c2, c3 = math.cos(ey), math.cos(ex), math.cos(ez)
    s1, s2, s3 = math.sin(ey), math.sin(ex), math.sin(ez)
    return [
        [c1 * c3 - s1 * s2 * s3, c2 * s3, s1 * c3 + c1 * s2 * s3],
        [-c1 * s3 - s1 * s2 * c3, c2 * c3, -s1 * s3 + c1 * s2 * c3],
        [-s1 * c2, -s2, c1 * c2],
    ]


def rotation_error_deg(a, b):
    ra = euler_to_rotation(a)
    rb = euler_to_rotation(b)
    trace = 0.0
    for i in range(3):
        for j in range(3):
            trace += ra[i][j] * rb[i][j]
    cos_angle = max(-1.0, min(1.0, (trace - 1.0) / 2.0))
    return math.degrees(math.acos(cos_angle))


def rmse_vectors(a, b):
    n = min(len(a), len(b))
    if n == 0:
        return math.nan
    total = 0.0
    for i in range(n):
        total += sum((a[i][j] - b[i][j]) ** 2 for j in range(3))
    return math.sqrt(total / n)


def read_ply_xyz(path):
    path = Path(path)
    if not path.exists():
        return []
    rows = []
    in_body = False
    for line in path.read_text(errors="replace").splitlines():
        if in_body:
            parts = line.split()
            if len(parts) >= 3:
                rows.append([float(parts[0]), float(parts[1]), float(parts[2])])
        elif line.strip() == "end_header":
            in_body = True
    return rows


def read_output_cam(path):
    path = Path(path)
    if not path.exists():
        return []
    return [row[:6] for row in read_numeric_rows(path)]


def generate_variant(gt_dir, out_dir, label, rot_noise_deg, cam_noise_m, point_noise_m, seed):
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)

    shutil.copy2(gt_dir / "cal.txt", out_dir / "cal.txt")
    shutil.copy2(gt_dir / "Feature.txt", out_dir / "Feature.txt")

    rot_sigma = math.radians(rot_noise_deg)
    cam_gt = read_numeric_rows(gt_dir / "Cam.txt")
    xyz_gt = read_numeric_rows(gt_dir / "XYZ.txt")

    cam_noisy = []
    for row in cam_gt:
        noisy = row[:]
        for i in range(3):
            noisy[i] += rng.gauss(0.0, rot_sigma)
        for i in range(3, 6):
            noisy[i] += rng.gauss(0.0, cam_noise_m)
        cam_noisy.append(noisy)

    xyz_noisy = []
    for row in xyz_gt:
        xyz_noisy.append([row[i] + rng.gauss(0.0, point_noise_m) for i in range(3)])

    write_cam(out_dir / "Cam.txt", cam_noisy)
    write_xyz(out_dir / "XYZ.txt", xyz_noisy)
    return cam_gt, xyz_gt, cam_noisy, xyz_noisy


def parse_comparison_csv(path):
    with Path(path).open(newline="") as f:
        return list(csv.DictReader(f))


def method_output_stem(code_name):
    return f"{code_name}_euler_angle_uv"


def add_accuracy_metrics(row, variant_dir, gt_cam, gt_xyz, initial_cam, initial_xyz):
    code_name = row["code_name"]
    stem = method_output_stem(code_name)
    final_cam = read_output_cam(variant_dir / f"Cam_{stem}.txt")
    final_xyz = read_ply_xyz(variant_dir / f"XYZ_{stem}.ply")

    row["initial_point_rmse_m"] = rmse_vectors(initial_xyz, gt_xyz)
    row["final_point_rmse_m"] = rmse_vectors(final_xyz, gt_xyz)
    row["initial_camera_center_rmse_m"] = rmse_vectors(
        [r[3:6] for r in initial_cam], [r[3:6] for r in gt_cam]
    )
    row["final_camera_center_rmse_m"] = rmse_vectors(
        [r[3:6] for r in final_cam], [r[3:6] for r in gt_cam]
    )

    n_rot_init = min(len(initial_cam), len(gt_cam))
    n_rot_final = min(len(final_cam), len(gt_cam))
    row["initial_rotation_rmse_deg"] = math.sqrt(
        sum(rotation_error_deg(initial_cam[i], gt_cam[i]) ** 2 for i in range(n_rot_init))
        / max(1, n_rot_init)
    )
    row["final_rotation_rmse_deg"] = math.sqrt(
        sum(rotation_error_deg(final_cam[i], gt_cam[i]) ** 2 for i in range(n_rot_final))
        / max(1, n_rot_final)
    )


def mean(values):
    values = [float(v) for v in values if v != "" and not math.isnan(float(v))]
    return sum(values) / len(values) if values else math.nan


def aggregate(rows):
    grouped = defaultdict(list)
    for row in rows:
        grouped[(row["noise_label"], row["method"])].append(row)

    summaries = []
    for (noise_label, method), items in grouped.items():
        summaries.append(
            {
                "noise_label": noise_label,
                "method": method,
                "anchor_mode": items[0]["anchor_mode"],
                "code_name": items[0]["code_name"],
                "rot_noise_deg": items[0]["rot_noise_deg"],
                "camera_noise_m": items[0]["camera_noise_m"],
                "point_noise_m": items[0]["point_noise_m"],
                "runs": len(items),
                "convergence_rate": sum(i["termination"] == "CONVERGENCE" for i in items)
                / len(items),
                "mean_iterations": mean(i["iterations"] for i in items),
                "mean_time_sec": mean(i["total_time_sec"] for i in items),
                "mean_initial_rms_px": mean(i["initial_rms_px"] for i in items),
                "mean_final_rms_px": mean(i["final_rms_px"] for i in items),
                "mean_final_point_rmse_m": mean(i["final_point_rmse_m"] for i in items),
                "mean_final_camera_center_rmse_m": mean(
                    i["final_camera_center_rmse_m"] for i in items
                ),
                "mean_final_rotation_rmse_deg": mean(i["final_rotation_rmse_deg"] for i in items),
            }
        )
    summaries.sort(key=lambda r: (noise_index(r["noise_label"]), METHOD_ORDER.index(r["method"])))
    return summaries


def noise_index(label):
    for i, level in enumerate(NOISE_LEVELS):
        if level[0] == label:
            return i
    return 999


def write_csv(path, rows, fieldnames):
    with Path(path).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def fmt(value, digits=4):
    value = float(value)
    if math.isnan(value):
        return "nan"
    return f"{value:.{digits}g}"


def best_methods_for_noise(summaries, label):
    rows = [r for r in summaries if r["noise_label"] == label]
    return sorted(rows, key=lambda r: (float(r["mean_final_rms_px"]), float(r["mean_time_sec"])))


def write_report(path, summaries, all_rows, levels, dataset_label):
    path = Path(path)
    with path.open("w", newline="\n") as f:
        f.write("# BA Noise Parameterization Experiment\n\n")
        f.write(f"Dataset: `{dataset_label}`\n\n")
        f.write(
            "Noise is added to the initial `Cam.txt` Euler angles/camera centers and "
            "initial `XYZ.txt` point coordinates. `Feature.txt` and `cal.txt` are kept fixed.\n\n"
        )
        f.write("## Noise Levels\n\n")
        f.write("| Level | Rotation sigma(deg) | Camera sigma(m) | Point sigma(m) | Seeds |\n")
        f.write("| --- | ---: | ---: | ---: | ---: |\n")
        for label, rot, cam, point, seeds in levels:
            f.write(f"| {label} | {rot:g} | {cam:g} | {point:g} | {len(seeds)} |\n")

        f.write("\n## Winners By Final Reprojection RMS\n\n")
        f.write("| Noise | Best method | Anchor | Final RMS(px) | Iter | Time(s) | Point RMSE(m) | Camera RMSE(m) | Rot RMSE(deg) |\n")
        f.write("| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for label, *_ in levels:
            best = best_methods_for_noise(summaries, label)[0]
            f.write(
                f"| {label} | {best['method']} | {best['anchor_mode']} | "
                f"{fmt(best['mean_final_rms_px'])} | {fmt(best['mean_iterations'])} | "
                f"{fmt(best['mean_time_sec'])} | {fmt(best['mean_final_point_rmse_m'])} | "
                f"{fmt(best['mean_final_camera_center_rmse_m'])} | "
                f"{fmt(best['mean_final_rotation_rmse_deg'])} |\n"
            )

        f.write("\n## Full Method Summary\n\n")
        f.write("| Noise | Anchor | Method | Code | Conv. | Iter | Time(s) | Init RMS(px) | Final RMS(px) | Point RMSE(m) | Cam RMSE(m) | Rot RMSE(deg) |\n")
        f.write("| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in summaries:
            f.write(
                f"| {row['noise_label']} | {row['anchor_mode']} | {row['method']} | {row['code_name']} | "
                f"{fmt(row['convergence_rate'])} | {fmt(row['mean_iterations'])} | "
                f"{fmt(row['mean_time_sec'])} | {fmt(row['mean_initial_rms_px'])} | "
                f"{fmt(row['mean_final_rms_px'])} | {fmt(row['mean_final_point_rmse_m'])} | "
                f"{fmt(row['mean_final_camera_center_rmse_m'])} | "
                f"{fmt(row['mean_final_rotation_rmse_deg'])} |\n"
            )

        f.write("\n## Interpretation Notes\n\n")
        f.write(
            "- `world_xy_inverse_z` and `anchored_xy_inverse_z` use inverse of the world/anchor Z component, "
            "not bearing plus inverse range. They are coordinate-axis dependent.\n"
        )
        f.write(
            "- `world_bearing_inverse_range` is the CEP historical `inverse_depth` code path, "
            "but geometrically it is bearing plus inverse Euclidean range in the world frame.\n"
        )
        f.write(
            "- `anchored_bearing_inverse_range` and `dual_anchor_bearing_parallax` normally describe nearly the same "
            "depth weak direction. Parallax adds an associate anchor and expresses that weak direction as an angle.\n"
        )

        by_noise = defaultdict(list)
        for row in summaries:
            by_noise[row["noise_label"]].append(row)
        f.write("\n## Compact Conclusions\n\n")
        for label, *_ in levels:
            rows = sorted(by_noise[label], key=lambda r: float(r["mean_final_rms_px"]))
            best = rows[0]
            runners = ", ".join(r["method"] for r in rows[:3])
            f.write(
                f"- `{label}`: best final RMS is `{best['method']}` "
                f"({fmt(best['mean_final_rms_px'])} px). Top three: {runners}.\n"
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gt-dir",
        default=r"C:\zuo\Projects\BA Datasets\Close-Range\CR1-problem-11-9611\Ground Truth",
    )
    parser.add_argument(
        "--out-root",
        default=r"C:\zuo\Projects\BA Datasets\Close-Range\CR1-problem-11-9611\Noise Experiments",
    )
    parser.add_argument(
        "--example-exe",
        default=r"C:\zuo\Projects\repos\CEP\example_v2\build\example.exe",
    )
    parser.add_argument("--force", action="store_true")
    parser.add_argument(
        "--levels",
        default=",".join(level[0] for level in NOISE_LEVELS),
        help="Comma-separated noise labels to run, for example: clean,large,severe",
    )
    parser.add_argument(
        "--timeout-sec",
        type=int,
        default=300,
        help="Timeout for each dataset variant. One variant runs all 9 BA methods.",
    )
    args = parser.parse_args()

    gt_dir = Path(args.gt_dir)
    out_root = Path(args.out_root)
    example_exe = Path(args.example_exe)
    selected_names = {name.strip() for name in args.levels.split(",") if name.strip()}
    selected_levels = [level for level in NOISE_LEVELS if level[0] in selected_names]
    unknown_levels = selected_names - {level[0] for level in NOISE_LEVELS}
    if unknown_levels:
        raise ValueError(f"Unknown levels: {', '.join(sorted(unknown_levels))}")
    if not selected_levels:
        raise ValueError("No noise levels selected")

    if args.force and out_root.exists():
        shutil.rmtree(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    all_rows = []
    for label, rot_noise_deg, cam_noise_m, point_noise_m, seeds in selected_levels:
        for seed in seeds:
            variant_dir = out_root / f"{label}_seed{seed}"
            gt_cam, gt_xyz, initial_cam, initial_xyz = generate_variant(
                gt_dir, variant_dir, label, rot_noise_deg, cam_noise_m, point_noise_m, seed
            )

            log_path = variant_dir / "run.log"
            print(f"Running {label} seed {seed}: {variant_dir}")
            with log_path.open("w", encoding="utf-8", errors="replace") as log:
                completed = subprocess.run(
                    [str(example_exe), str(variant_dir / "cal.txt")],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    text=True,
                    cwd=str(example_exe.parent),
                    timeout=args.timeout_sec,
                )
            if completed.returncode != 0:
                raise RuntimeError(f"example.exe failed for {variant_dir}, see {log_path}")

            comparison = parse_comparison_csv(variant_dir / "BA-comparison.csv")
            for row in comparison:
                row["noise_label"] = label
                row["seed"] = seed
                row["rot_noise_deg"] = rot_noise_deg
                row["camera_noise_m"] = cam_noise_m
                row["point_noise_m"] = point_noise_m
                row["dataset_dir"] = str(variant_dir)
                add_accuracy_metrics(row, variant_dir, gt_cam, gt_xyz, initial_cam, initial_xyz)
                all_rows.append(row)

    run_fields = [
        "noise_label",
        "seed",
        "rot_noise_deg",
        "camera_noise_m",
        "point_noise_m",
        "anchor_mode",
        "method",
        "code_name",
        "termination",
        "iterations",
        "total_time_sec",
        "initial_rms_px",
        "final_rms_px",
        "initial_point_rmse_m",
        "final_point_rmse_m",
        "initial_camera_center_rmse_m",
        "final_camera_center_rmse_m",
        "initial_rotation_rmse_deg",
        "final_rotation_rmse_deg",
        "dataset_dir",
    ]
    write_csv(out_root / "noise-experiment-runs.csv", all_rows, run_fields)

    summaries = aggregate(all_rows)
    summary_fields = [
        "noise_label",
        "rot_noise_deg",
        "camera_noise_m",
        "point_noise_m",
        "anchor_mode",
        "method",
        "code_name",
        "runs",
        "convergence_rate",
        "mean_iterations",
        "mean_time_sec",
        "mean_initial_rms_px",
        "mean_final_rms_px",
        "mean_final_point_rmse_m",
        "mean_final_camera_center_rmse_m",
        "mean_final_rotation_rmse_deg",
    ]
    write_csv(out_root / "noise-experiment-summary.csv", summaries, summary_fields)
    write_report(out_root / "noise-experiment-report.md", summaries, all_rows, selected_levels, str(gt_dir))
    print(f"Wrote {out_root / 'noise-experiment-report.md'}")


if __name__ == "__main__":
    main()
