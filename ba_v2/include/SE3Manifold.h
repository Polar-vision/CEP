// #include <ceres/manifold.h>
// #include <sophus/se3.hpp>
// #include <Eigen/Core>

// class SE3Manifold : public ceres::Manifold {
// public:
//     SE3Manifold() {}
//     virtual ~SE3Manifold() {}

//     bool Plus(const double* x, const double* delta, double* x_plus_delta) const override {
//         Eigen::Map<const Eigen::Matrix<double, 6, 1>> lie(x);
//         Eigen::Map<const Eigen::Matrix<double, 6, 1>> delta_lie(delta);

//         Sophus::SE3d T = Sophus::SE3d::exp(lie);
//         Sophus::SE3d delta_T = Sophus::SE3d::exp(delta_lie);

//         Eigen::Matrix<double, 6, 1> x_plus_delta_lie = (delta_T * T).log();//李代数平移+轴角
//         Eigen::Map<Eigen::Matrix<double, 6, 1>> x_plus_delta_map(x_plus_delta);
//         x_plus_delta_map = x_plus_delta_lie;

//         return true;
//     }

//     bool PlusJacobian(const double* x, double* jacobian) const override {
//         // Eigen::Map<const Eigen::Matrix<double, 6, 1>> xi(x);
//         // Eigen::Matrix<double, 6, 6> J = Sophus::SE3d::leftJacobian(xi);
//         // Eigen::Map<Eigen::Matrix<double, 6, 6, Eigen::RowMajor>> J_out(jacobian);
//         // J_out = J;
//         // return true;

//         Eigen::Map<const Eigen::Matrix<double, 6, 1>> lie(x);
//         Sophus::SE3d T = Sophus::SE3d::exp(lie);

//         // Adjoint(T^{-1}) = Adj(T.inverse())
//         Sophus::Matrix6d Ad_T_inv = T.inverse().Adj();

//         // Ceres 默认使用行优先存储
//         Eigen::Map<Eigen::Matrix<double, 6, 6, Eigen::RowMajor>> jacobian_(jacobian);
//         jacobian_ = Ad_T_inv;
//         return true;
//     }

//     bool Minus(const double* y, const double* x, double* y_minus_x) const override {
//         // Eigen::Map<const Eigen::Matrix<double, 6, 1>> y_lie(y);
//         // Eigen::Map<const Eigen::Matrix<double, 6, 1>> x_lie(x);
//         // Sophus::SE3d Y = Sophus::SE3d::exp(y_lie);
//         // Sophus::SE3d X = Sophus::SE3d::exp(x_lie);
//         // Eigen::Matrix<double, 6, 1> delta = (Y * X.inverse()).log();
//         // Eigen::Map<Eigen::Matrix<double, 6, 1>> y_minus_x_map(y_minus_x);
//         // y_minus_x_map = delta;
//         return true;
//     }

//     bool MinusJacobian(const double* x, double* jacobian) const override {
//         // Eigen::Map<const Eigen::Matrix<double, 6, 1>> xi(x);
//         // Eigen::Matrix<double, 6, 6> J = Sophus::SE3d::leftJacobianInv(xi);
//         // Eigen::Map<Eigen::Matrix<double, 6, 6, Eigen::RowMajor>> J_out(jacobian);
//         // J_out = J;
//         return true;
//     }

//     int AmbientSize() const override { return 6; }
//     int TangentSize() const override { return 6; }
// };
