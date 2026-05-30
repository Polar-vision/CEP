#include "PBAImp_v2.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr double kPi = 3.14159265358979323846;

struct CameraPose {
  std::array<double, 3> euler;
  std::array<double, 3> center;
};

struct Point3 {
  double x;
  double y;
  double z;
};

struct Observation {
  int camera;
  double u;
  double v;
};

struct Track {
  Point3 truth;
  Point3 initial;
  std::vector<Observation> observations;
};

struct Scene {
  std::string name;
  double baseline_step;
  double noise_px;
  std::vector<CameraPose> cameras;
  std::vector<Track> tracks;
};

struct Metrics {
  std::string parameterization;
  double initial_rms_px = 0.0;
  double final_rms_px = 0.0;
  double point_rmse_m = 0.0;
  double median_depth_relative_error = 0.0;
  int iterations = 0;
  std::string termination;
};

double Squared(double x) { return x * x; }

void EulerToRotation(const std::array<double, 3>& euler, double R[9]) {
  const double ey = euler[0];
  const double ex = euler[1];
  const double ez = euler[2];
  const double c1 = std::cos(ey);
  const double c2 = std::cos(ex);
  const double c3 = std::cos(ez);
  const double s1 = std::sin(ey);
  const double s2 = std::sin(ex);
  const double s3 = std::sin(ez);
  R[0] = c1 * c3 - s1 * s2 * s3;
  R[1] = c2 * s3;
  R[2] = s1 * c3 + c1 * s2 * s3;
  R[3] = -c1 * s3 - s1 * s2 * c3;
  R[4] = c2 * c3;
  R[5] = -s1 * s3 + c1 * s2 * c3;
  R[6] = -s1 * c2;
  R[7] = -s2;
  R[8] = c1 * c2;
}

std::array<double, 2> Project(const CameraPose& camera,
                              const Point3& point,
                              double fx,
                              double fy,
                              double cx,
                              double cy) {
  double R[9];
  EulerToRotation(camera.euler, R);
  const double Xc[3] = {
      point.x - camera.center[0],
      point.y - camera.center[1],
      point.z - camera.center[2],
  };
  const double p[3] = {
      R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2],
      R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2],
      R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2],
  };
  return {fx * p[0] / p[2] + cx, fy * p[1] / p[2] + cy};
}

Point3 DirectionDepthToPoint(const double direction_depth[3]) {
  const double da = direction_depth[0];
  const double ha = direction_depth[1];
  const double depth = direction_depth[2];
  return {std::sin(da) * std::cos(ha) * depth,
          std::sin(ha) * depth,
          std::cos(da) * std::cos(ha) * depth};
}

Point3 DirectionInverseDepthToPoint(const double direction_inverse_depth[3]) {
  const double da = direction_inverse_depth[0];
  const double ha = direction_inverse_depth[1];
  const double inv_depth = direction_inverse_depth[2];
  return {std::sin(da) * std::cos(ha) / inv_depth,
          std::sin(ha) / inv_depth,
          std::cos(da) * std::cos(ha) / inv_depth};
}

Point3 ParallaxToPoint(const double direction[2],
                       double parallax,
                       const std::array<double, 3>& main_center,
                       const std::array<double, 3>& associate_center) {
  const std::array<double, 3> ray = {
      std::sin(direction[0]) * std::cos(direction[1]),
      std::sin(direction[1]),
      std::cos(direction[0]) * std::cos(direction[1]),
  };
  const std::array<double, 3> baseline = {
      associate_center[0] - main_center[0],
      associate_center[1] - main_center[1],
      associate_center[2] - main_center[2],
  };
  const double baseline_norm =
      std::sqrt(Squared(baseline[0]) + Squared(baseline[1]) + Squared(baseline[2]));
  const double dot = ray[0] * baseline[0] + ray[1] * baseline[1] + ray[2] * baseline[2];
  const double w2 = std::acos(std::clamp(dot / baseline_norm, -1.0, 1.0));
  const double scale = baseline_norm * std::sin(w2 + parallax) / std::sin(parallax);
  return {main_center[0] + scale * ray[0],
          main_center[1] + scale * ray[1],
          main_center[2] + scale * ray[2]};
}

double ComputeReprojectionRms(const Scene& scene,
                              const std::vector<Point3>& points,
                              double fx,
                              double fy,
                              double cx,
                              double cy) {
  double squared_sum = 0.0;
  int residual_count = 0;
  for (size_t i = 0; i < scene.tracks.size(); ++i) {
    for (const auto& obs : scene.tracks[i].observations) {
      const auto projected = Project(scene.cameras[obs.camera], points[i], fx, fy, cx, cy);
      squared_sum += Squared(projected[0] - obs.u) + Squared(projected[1] - obs.v);
      residual_count += 2;
    }
  }
  return std::sqrt(squared_sum / residual_count);
}

double ComputePointRmse(const Scene& scene, const std::vector<Point3>& points) {
  double squared_sum = 0.0;
  for (size_t i = 0; i < scene.tracks.size(); ++i) {
    squared_sum += Squared(points[i].x - scene.tracks[i].truth.x);
    squared_sum += Squared(points[i].y - scene.tracks[i].truth.y);
    squared_sum += Squared(points[i].z - scene.tracks[i].truth.z);
  }
  return std::sqrt(squared_sum / scene.tracks.size());
}

double MedianDepthRelativeError(const Scene& scene, const std::vector<Point3>& points) {
  std::vector<double> errors;
  errors.reserve(points.size());
  for (size_t i = 0; i < points.size(); ++i) {
    errors.push_back(std::abs(points[i].z - scene.tracks[i].truth.z) / scene.tracks[i].truth.z);
  }
  std::sort(errors.begin(), errors.end());
  return errors[errors.size() / 2];
}

Scene GenerateScene(const std::string& name,
                    double baseline_step,
                    double noise_px,
                    unsigned seed,
                    int camera_count,
                    int point_count,
                    double fx,
                    double fy,
                    double cx,
                    double cy) {
  Scene scene;
  scene.name = name;
  scene.baseline_step = baseline_step;
  scene.noise_px = noise_px;

  const double center_offset = 0.5 * (camera_count - 1);
  for (int i = 0; i < camera_count; ++i) {
    CameraPose camera;
    camera.euler = {0.0, 0.0, 0.0};
    camera.center = {(i - center_offset) * baseline_step, 0.0, 0.0};
    scene.cameras.push_back(camera);
  }

  std::mt19937 rng(seed);
  std::uniform_real_distribution<double> xy_distribution(-2.0, 2.0);
  std::uniform_real_distribution<double> z_distribution(35.0, 65.0);
  std::normal_distribution<double> image_noise(0.0, noise_px);
  std::normal_distribution<double> lateral_noise(0.0, 0.25);
  std::normal_distribution<double> depth_noise(0.0, 8.0);

  for (int i = 0; i < point_count; ++i) {
    Track track;
    track.truth = {xy_distribution(rng), xy_distribution(rng), z_distribution(rng)};
    track.initial = {track.truth.x + lateral_noise(rng),
                     track.truth.y + lateral_noise(rng),
                     std::max(5.0, track.truth.z + depth_noise(rng))};

    for (int camera_idx = 0; camera_idx < camera_count; ++camera_idx) {
      const auto projected = Project(scene.cameras[camera_idx], track.truth, fx, fy, cx, cy);
      track.observations.push_back(
          {camera_idx, projected[0] + image_noise(rng), projected[1] + image_noise(rng)});
    }
    scene.tracks.push_back(track);
  }

  return scene;
}

void WriteSceneFiles(const Scene& scene,
                     const fs::path& directory,
                     double fx,
                     double fy,
                     double cx,
                     double cy) {
  fs::create_directories(directory);

  {
    std::ofstream cal(directory / "cal.txt");
    cal << std::fixed << std::setprecision(9);
    cal << fx << " 0 " << cx << "\n";
    cal << "0 " << fy << " " << cy << "\n";
    cal << "0 0 1\n";
  }

  {
    std::ofstream cam(directory / "Cam.txt");
    cam << std::fixed << std::setprecision(9);
    for (const auto& camera : scene.cameras) {
      cam << camera.euler[0] << " " << camera.euler[1] << " " << camera.euler[2] << " "
          << camera.center[0] << " " << camera.center[1] << " " << camera.center[2] << " 1\n";
    }
  }

  {
    std::ofstream xyz(directory / "XYZ.txt");
    xyz << std::fixed << std::setprecision(9);
    for (const auto& track : scene.tracks) {
      xyz << track.initial.x << " " << track.initial.y << " " << track.initial.z << "\n";
    }
  }

  {
    std::ofstream truth(directory / "XYZ_truth.txt");
    truth << std::fixed << std::setprecision(9);
    for (const auto& track : scene.tracks) {
      truth << track.truth.x << " " << track.truth.y << " " << track.truth.z << "\n";
    }
  }

  {
    std::ofstream feature(directory / "Feature.txt");
    feature << std::fixed << std::setprecision(9);
    for (const auto& track : scene.tracks) {
      feature << track.observations.size();
      for (const auto& obs : track.observations) {
        feature << " " << obs.camera << " " << obs.u << " " << obs.v;
      }
      feature << "\n";
    }
  }
}

Metrics SolveWithXyz(const Scene& scene, double fx, double fy, double cx, double cy) {
  Metrics metrics;
  metrics.parameterization = "xyz";
  std::vector<Point3> points;
  points.reserve(scene.tracks.size());
  for (const auto& track : scene.tracks) {
    points.push_back(track.initial);
  }
  metrics.initial_rms_px = ComputeReprojectionRms(scene, points, fx, fy, cx, cy);

  std::vector<std::array<double, 3>> eulers;
  std::vector<std::array<double, 3>> centers;
  eulers.reserve(scene.cameras.size());
  centers.reserve(scene.cameras.size());
  for (const auto& camera : scene.cameras) {
    eulers.push_back(camera.euler);
    centers.push_back(camera.center);
  }

  ceres::Problem problem;
  for (size_t i = 0; i < scene.tracks.size(); ++i) {
    for (const auto& obs : scene.tracks[i].observations) {
      auto* cost = PBA::xyz_euler_angle_uv::Create(obs.u, obs.v, fx, fy, cx, cy);
      problem.AddResidualBlock(
          cost, nullptr, eulers[obs.camera].data(), centers[obs.camera].data(), &points[i].x);
    }
  }
  for (size_t i = 0; i < scene.cameras.size(); ++i) {
    problem.SetParameterBlockConstant(eulers[i].data());
    problem.SetParameterBlockConstant(centers[i].data());
  }

  ceres::Solver::Options options;
  options.logging_type = ceres::SILENT;
  options.linear_solver_type = ceres::DENSE_QR;
  options.max_num_iterations = 100;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);
  metrics.iterations = static_cast<int>(summary.iterations.size());
  metrics.termination = summary.termination_type == ceres::CONVERGENCE ? "CONVERGENCE" : "NO_CONVERGENCE";
  metrics.final_rms_px = ComputeReprojectionRms(scene, points, fx, fy, cx, cy);
  metrics.point_rmse_m = ComputePointRmse(scene, points);
  metrics.median_depth_relative_error = MedianDepthRelativeError(scene, points);
  return metrics;
}

Metrics SolveWithInverseDepth(const Scene& scene, double fx, double fy, double cx, double cy) {
  Metrics metrics;
  metrics.parameterization = "inverse_depth";
  std::vector<std::array<double, 3>> params;
  std::vector<Point3> initial_points;
  params.reserve(scene.tracks.size());
  initial_points.reserve(scene.tracks.size());
  for (const auto& track : scene.tracks) {
    initial_points.push_back(track.initial);
    const double x = track.initial.x;
    const double y = track.initial.y;
    const double z = track.initial.z;
    const double depth = std::sqrt(Squared(x) + Squared(y) + Squared(z));
    params.push_back({std::atan2(x, z), std::atan2(y, std::sqrt(Squared(x) + Squared(z))), 1.0 / depth});
  }
  metrics.initial_rms_px = ComputeReprojectionRms(scene, initial_points, fx, fy, cx, cy);

  std::vector<std::array<double, 3>> eulers;
  std::vector<std::array<double, 3>> centers;
  eulers.reserve(scene.cameras.size());
  centers.reserve(scene.cameras.size());
  for (const auto& camera : scene.cameras) {
    eulers.push_back(camera.euler);
    centers.push_back(camera.center);
  }

  ceres::Problem problem;
  for (size_t i = 0; i < scene.tracks.size(); ++i) {
    for (const auto& obs : scene.tracks[i].observations) {
      auto* cost = PBA::inverse_depth_euler_angle_uv::Create(obs.u, obs.v, fx, fy, cx, cy);
      problem.AddResidualBlock(
          cost, nullptr, eulers[obs.camera].data(), centers[obs.camera].data(), params[i].data());
    }
  }
  for (size_t i = 0; i < scene.cameras.size(); ++i) {
    problem.SetParameterBlockConstant(eulers[i].data());
    problem.SetParameterBlockConstant(centers[i].data());
  }

  ceres::Solver::Options options;
  options.logging_type = ceres::SILENT;
  options.linear_solver_type = ceres::DENSE_QR;
  options.max_num_iterations = 100;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  std::vector<Point3> points;
  points.reserve(params.size());
  for (const auto& param : params) {
    points.push_back(DirectionInverseDepthToPoint(param.data()));
  }

  metrics.iterations = static_cast<int>(summary.iterations.size());
  metrics.termination = summary.termination_type == ceres::CONVERGENCE ? "CONVERGENCE" : "NO_CONVERGENCE";
  metrics.final_rms_px = ComputeReprojectionRms(scene, points, fx, fy, cx, cy);
  metrics.point_rmse_m = ComputePointRmse(scene, points);
  metrics.median_depth_relative_error = MedianDepthRelativeError(scene, points);
  return metrics;
}

Metrics SolveWithParallax(const Scene& scene, double fx, double fy, double cx, double cy) {
  Metrics metrics;
  metrics.parameterization = "parallax";
  const int main_anchor = 0;
  const int associate_anchor = static_cast<int>(scene.cameras.size()) - 1;

  std::vector<std::array<double, 2>> directions;
  std::vector<std::array<double, 1>> parallaxes;
  std::vector<Point3> initial_points;
  directions.reserve(scene.tracks.size());
  parallaxes.reserve(scene.tracks.size());
  initial_points.reserve(scene.tracks.size());
  for (const auto& track : scene.tracks) {
    initial_points.push_back(track.initial);
    const auto& main_center = scene.cameras[main_anchor].center;
    const auto& associate_center = scene.cameras[associate_anchor].center;
    const std::array<double, 3> main_ray = {
        track.initial.x - main_center[0],
        track.initial.y - main_center[1],
        track.initial.z - main_center[2],
    };
    const std::array<double, 3> associate_ray = {
        track.initial.x - associate_center[0],
        track.initial.y - associate_center[1],
        track.initial.z - associate_center[2],
    };
    const double main_norm =
        std::sqrt(Squared(main_ray[0]) + Squared(main_ray[1]) + Squared(main_ray[2]));
    const double associate_norm =
        std::sqrt(Squared(associate_ray[0]) + Squared(associate_ray[1]) + Squared(associate_ray[2]));
    const double dot = main_ray[0] * associate_ray[0] + main_ray[1] * associate_ray[1] +
                       main_ray[2] * associate_ray[2];
    directions.push_back({std::atan2(main_ray[0], main_ray[2]),
                          std::atan2(main_ray[1], std::sqrt(Squared(main_ray[0]) + Squared(main_ray[2])))});
    parallaxes.push_back({std::acos(std::clamp(dot / (main_norm * associate_norm), -1.0, 1.0))});
  }
  metrics.initial_rms_px = ComputeReprojectionRms(scene, initial_points, fx, fy, cx, cy);

  std::vector<std::array<double, 3>> eulers;
  std::vector<std::array<double, 3>> centers;
  eulers.reserve(scene.cameras.size());
  centers.reserve(scene.cameras.size());
  for (const auto& camera : scene.cameras) {
    eulers.push_back(camera.euler);
    centers.push_back(camera.center);
  }

  ceres::Problem problem;
  for (size_t i = 0; i < scene.tracks.size(); ++i) {
    for (const auto& obs : scene.tracks[i].observations) {
      if (obs.camera == main_anchor) {
        auto* cost = PBA::parallax_euler_angle_uv_nM::Create(obs.u, obs.v, fx, fy, cx, cy);
        problem.AddResidualBlock(cost, nullptr, eulers[obs.camera].data(), directions[i].data());
      } else if (obs.camera == associate_anchor) {
        auto* cost = PBA::parallax_euler_angle_uv_nA::Create(obs.u, obs.v, fx, fy, cx, cy);
        problem.AddResidualBlock(cost,
                                 nullptr,
                                 eulers[obs.camera].data(),
                                 centers[obs.camera].data(),
                                 centers[main_anchor].data(),
                                 directions[i].data(),
                                 parallaxes[i].data());
      } else {
        auto* cost = PBA::parallax_euler_angle_uv_nP::Create(obs.u, obs.v, fx, fy, cx, cy);
        problem.AddResidualBlock(cost,
                                 nullptr,
                                 eulers[obs.camera].data(),
                                 centers[obs.camera].data(),
                                 centers[main_anchor].data(),
                                 centers[associate_anchor].data(),
                                 directions[i].data(),
                                 parallaxes[i].data());
      }
    }
  }
  for (size_t i = 0; i < scene.cameras.size(); ++i) {
    problem.SetParameterBlockConstant(eulers[i].data());
    problem.SetParameterBlockConstant(centers[i].data());
  }

  ceres::Solver::Options options;
  options.logging_type = ceres::SILENT;
  options.linear_solver_type = ceres::DENSE_QR;
  options.max_num_iterations = 100;
  ceres::Solver::Summary summary;
  ceres::Solve(options, &problem, &summary);

  std::vector<Point3> points;
  points.reserve(directions.size());
  for (size_t i = 0; i < directions.size(); ++i) {
    points.push_back(ParallaxToPoint(
        directions[i].data(), parallaxes[i][0], centers[main_anchor], centers[associate_anchor]));
  }

  metrics.iterations = static_cast<int>(summary.iterations.size());
  metrics.termination = summary.termination_type == ceres::CONVERGENCE ? "CONVERGENCE" : "NO_CONVERGENCE";
  metrics.final_rms_px = ComputeReprojectionRms(scene, points, fx, fy, cx, cy);
  metrics.point_rmse_m = ComputePointRmse(scene, points);
  metrics.median_depth_relative_error = MedianDepthRelativeError(scene, points);
  return metrics;
}

void PrintMetrics(const Scene& scene, const std::vector<Metrics>& metrics) {
  const double total_span = scene.baseline_step * (scene.cameras.size() - 1);
  const double mean_depth = std::accumulate(
      scene.tracks.begin(), scene.tracks.end(), 0.0, [](double s, const Track& t) {
        return s + t.truth.z;
      }) / scene.tracks.size();
  const double approx_endpoint_parallax_deg = std::atan2(total_span, mean_depth) * 180.0 / kPi;

  std::cout << "\nScene: " << scene.name << "\n";
  std::cout << "  baseline step: " << scene.baseline_step << " m, total span: " << total_span
            << " m, approx endpoint parallax: " << approx_endpoint_parallax_deg
            << " deg, noise: " << scene.noise_px << " px\n";
  std::cout << "  method              init_px   final_px   point_rmse_m   median_depth_err   iters  status\n";
  for (const auto& m : metrics) {
    std::cout << "  " << std::left << std::setw(18) << m.parameterization << std::right
              << std::setw(8) << std::setprecision(4) << std::fixed << m.initial_rms_px
              << std::setw(11) << m.final_rms_px << std::setw(15) << m.point_rmse_m
              << std::setw(19) << m.median_depth_relative_error << std::setw(8)
              << m.iterations << "  " << m.termination << "\n";
  }
}

void WriteMetricsCsvRows(std::ofstream& csv,
                         const Scene& scene,
                         const std::vector<Metrics>& metrics) {
  const double total_span = scene.baseline_step * (scene.cameras.size() - 1);
  const double mean_depth = std::accumulate(
      scene.tracks.begin(), scene.tracks.end(), 0.0, [](double s, const Track& t) {
        return s + t.truth.z;
      }) / scene.tracks.size();
  const double approx_endpoint_parallax_deg = std::atan2(total_span, mean_depth) * 180.0 / kPi;

  for (const auto& m : metrics) {
    csv << scene.name << "," << scene.baseline_step << "," << total_span << ","
        << approx_endpoint_parallax_deg << "," << scene.noise_px << ","
        << m.parameterization << "," << m.initial_rms_px << "," << m.final_rms_px << ","
        << m.point_rmse_m << "," << m.median_depth_relative_error << ","
        << m.iterations << "," << m.termination << "\n";
  }
}

}  // namespace

int main(int argc, char* argv[]) {
  (void)argc;

  const double fx = 800.0;
  const double fy = 800.0;
  const double cx = 640.0;
  const double cy = 480.0;
  const double noise_px = 0.3;
  const int camera_count = 5;
  const int point_count = 80;

  const fs::path exe_dir = fs::absolute(fs::path(argv[0])).parent_path();
  const fs::path output_dir = exe_dir / "synthetic_data";
  fs::create_directories(output_dir);

  std::ofstream csv(output_dir / "summary.csv");
  csv << "scene,baseline_step_m,total_span_m,approx_endpoint_parallax_deg,noise_px,"
         "parameterization,initial_rms_px,final_rms_px,point_rmse_m,"
         "median_depth_relative_error,iterations,termination\n";

  const std::vector<Scene> scenes = {
      GenerateScene("normal_baseline", 1.0, noise_px, 7, camera_count, point_count, fx, fy, cx, cy),
      GenerateScene("short_baseline", 0.02, noise_px, 7, camera_count, point_count, fx, fy, cx, cy),
      GenerateScene("ultra_short_baseline", 0.005, noise_px, 7, camera_count, point_count, fx, fy, cx, cy),
  };

  for (const auto& scene : scenes) {
    WriteSceneFiles(scene, output_dir / scene.name, fx, fy, cx, cy);
    std::vector<Metrics> metrics;
    metrics.push_back(SolveWithXyz(scene, fx, fy, cx, cy));
    metrics.push_back(SolveWithInverseDepth(scene, fx, fy, cx, cy));
    metrics.push_back(SolveWithParallax(scene, fx, fy, cx, cy));
    PrintMetrics(scene, metrics);
    WriteMetricsCsvRows(csv, scene, metrics);
  }

  std::cout << "\nSynthetic input files were written to: " << output_dir.string() << "\n";
  std::cout << "Summary CSV: " << (output_dir / "summary.csv").string() << "\n";
  return 0;
}
