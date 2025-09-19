#pragma once
#include "BAExporter_v2.h"
#include "IBA_v2.h"
#include <map>
#include <vector>
#include <set>
#include <algorithm>
// #include "SE3Manifold.h"
#include "sophus/ceres_manifold.hpp"

using namespace std;

class PBA : public IBA
{
public:
	PBA(void);
	~PBA(void);
	virtual bool ba_run(char* szCam,
		                char* szFea, 
						char* szXYZ, 
						char* szCalib, 
						char* szReport,
						char* szPose, 
						char* sz3D,
						objectpointtype optype,
    					rotation3dtype r3dtype,
						imagepointtype iptype,
						parametertype paramtype,
    					manifoldtype manitype);


	virtual bool ba_initialize( char* szCamera, char* szFeature, char* szCalib = NULL, char* szXYZ = NULL );

	void	pba_readAndInitialize( char *camsfname, char *ptsfname, char *calibfname, int *ncams, int *n3Dpts, int *n2Dprojs,double **motstruct, 
				double **imgpts, int **archor, char **vmask, char **umask,int **nphoto, int** nfeature, int** archorSort );
	void	pba_readProjectionAndInitilizeFeature(	FILE *fp, double *params, double *projs, char *vmask, int ncams, 
				int *archor,char* umask,int* nphoto, int* nfeature, int* archorSort );
	//initialize feature. You can provide xyz or system also provide them by itself;
	bool    pba_initializeMainArchor( double* imgpts, double* camera,double* K,double* feature, int nP, int FID, double* KR );
	bool    pba_initializeAssoArchor( double* imgpts, int* photo, double* camera,double* K,double* feature,int nMI, int nAI, int FID, bool bLast );
	bool	pba_initializeOtheArchors( double* imgpts, int* photo, double* camera,double* K,double* feature,int* archorSort,int nfeacout, int nOI, int FID );
	double lc_rep();
	void narrow_test();
	void parallax2xyz();
	void xy_inverse_z2xyz();
	void depth2xyz();
	void inverse_depth2xyz();
	void archored_inverse_depth2xyz();
	void archored_depth2xyz();
	void archored_xy_inverse_z2xyz();
	void archored_xyz2xyz();
	void IFeature();
	void get_yaw_from_polar();
	double costFuc();
	void randCam(ParameterType paramtype);
	//Parameterization of image point
	//rotation+translation
	struct lc_euler_angle_rc {
		lc_euler_angle_rc(double lco, double x, double y, double z)
			: lco(lco), x(x), y(y), z(z){
		}

	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T c1 = cos(ey);
		T c2 = cos(ex);
		T s1 = sin(ey);
		T s2 = sin(ex);
		R[0]=c1;       R[1]=T(0);        R[2]=s1;
		R[3]=-s1*s2;   R[4]=c2;          R[5]=c1*s2;
		R[6]=-s1*c2;   R[7]=-s2;         R[8]=c1*c2;

		// R_[0]=cos(ey);          R_[1]=T(0);       R_[2]=sin(ey);
		// R_[3]=-sin(ey)*sin(ex); R_[4]=cos(ex);    R_[5]=cos(ey)*sin(ex);
		// R_[6]=-sin(ey)*cos(ex); R_[7]=-sin(ex);   R_[8]=cos(ey)*cos(ex);

		T Xc[3];
		Xc[0] = T(x) - camera_center[0];
		Xc[1] = T(y) - camera_center[1];
		Xc[2] = T(z) - camera_center[2];

		T p[3];
		p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
		p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
		p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z

		T lcp = atan(sqrt(p[0] * p[0] + p[1] * p[1]) / p[2]);

		// Residual
		residuals[0] = lcp - T(lco);

		return true;
	}

	static ceres::CostFunction *Create(double lco, double x, double y, double z)
	{
		return new ceres::AutoDiffCostFunction<lc_euler_angle_rc, 1, 2, 3>(
			new lc_euler_angle_rc(lco, x, y, z));
	}

	double lco, x, y, z;
	};
	struct uv_euler_angle_rc {
	uv_euler_angle_rc(double observed_u, double observed_v,
								   double fx, double fy, double cx, double cy,
				   				   double x, double y, double z)
		: u(observed_u), v(observed_v), 
		fx(fx), fy(fy), cx(cx), cy(cy), 
		x(x), y(y), z(z) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T Xc[3];
		Xc[0] = T(x) - camera_center[0];
		Xc[1] = T(y) - camera_center[1];
		Xc[2] = T(z) - camera_center[2];

		// P = R * (X - C)
		T p[3];
		p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
		p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
		p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									  double fx, double fy, double cx, double cy,
								      double x, double y, double z) {		
		return new ceres::AutoDiffCostFunction<uv_euler_angle_rc, 2, 3, 3>(
			new uv_euler_angle_rc(u, v, fx, fy, cx, cy, x, y, z));
	}

	double u, v;
	double fx, fy, cx, cy;
	double x, y, z;
	};
	
	//rotation
	struct lc_euler_angle_r {
		lc_euler_angle_r(double lco, double x, double y, double z, double xc, double yc, double zc)
			: lco(lco), x(x), y(y), z(z), xc(xc), yc(yc), zc(zc){
		}

	template <typename T>
	bool operator()(const T* const euler_angles,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T c1 = cos(ey);
		T c2 = cos(ex);
		T s1 = sin(ey);
		T s2 = sin(ex);
		R[0]=c1;       R[1]=T(0);        R[2]=s1;
		R[3]=-s1*s2;   R[4]=c2;          R[5]=c1*s2;
		R[6]=-s1*c2;   R[7]=-s2;         R[8]=c1*c2;

		T Xc[3];
		Xc[0] = T(x) - T(xc);
		Xc[1] = T(y) - T(yc);
		Xc[2] = T(z) - T(zc);

		T p[3];
		p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
		p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
		p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z

		T lcp = atan(sqrt(p[0] * p[0] + p[1] * p[1]) / p[2]);

		// Residual
		residuals[0] = lcp - T(lco);

		return true;
	}

	static ceres::CostFunction *Create(double lco, double x, double y, double z, double xc, double yc, double zc)
	{
		return new ceres::AutoDiffCostFunction<lc_euler_angle_r, 1, 2>(
			new lc_euler_angle_r(lco, x, y, z, xc, yc, zc));
	}

	double lco, x, y, z, xc, yc, zc;
	};
	struct uv_euler_angle_r {
	uv_euler_angle_r(double observed_u, double observed_v,
								   double fx, double fy, double cx, double cy,
				   				   double x, double y, double z,
								   double xc, double yc, double zc)
		: u(observed_u), v(observed_v), 
		fx(fx), fy(fy), cx(cx), cy(cy), 
		x(x), y(y), z(z),
		xc(xc), yc(yc), zc(zc) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const euler_angles,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T Xc[3];
		Xc[0] = T(x) - T(xc);
		Xc[1] = T(y) - T(yc);
		Xc[2] = T(z) - T(zc);

		T p[3];
		p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];
		p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];
		p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];

		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									  double fx, double fy, double cx, double cy,
								      double x, double y, double z,
									  double xc, double yc, double zc) {		
		return new ceres::AutoDiffCostFunction<uv_euler_angle_r, 2, 3>(
			new uv_euler_angle_r(u, v, fx, fy, cx, cy, x, y, z, xc, yc, zc));
	}

	double u, v;
	double fx, fy, cx, cy;
	double x, y, z;
	double xc, yc, zc;
	};


	//Parameterization of object point
	//zero archor
	struct xyz_euler_angle_uv {
	xyz_euler_angle_uv(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center,
					const T* const point3D,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		// T c0, c1, c2, s0, s1, s2;
		// c0 = cos(euler_angles[0]);
		// c1 = cos(euler_angles[1]);
		// c2 = cos(euler_angles[2]);
		// s0 = sin(euler_angles[0]);
		// s1 = sin(euler_angles[1]);
		// s2 = sin(euler_angles[2]);
		// R[0] = c1 * c0;
		// R[1] = c1 * s0;
		// R[2] = -s1;
		// R[3] = s2 * s1 * c0 - c2 * s0;
		// R[4] = s2 * s1 * s0 + c2 * c0;
		// R[5] = s2 * c1;
		// R[6] = c2 * s1 * c0 + s2 * s0;
		// R[7] = c2 * s1 * s0 - s2 * c0;
		// R[8] = c2 * c1;

		T Xc[3];
		Xc[0] = point3D[0] - camera_center[0];
		Xc[1] = point3D[1] - camera_center[1];
		Xc[2] = point3D[2] - camera_center[2];

		// P = R * (X - C)
		T p[3];
		p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
		p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
		p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<xyz_euler_angle_uv, 2, 3, 3, 3>(
			new xyz_euler_angle_uv(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct xy_inverse_z_euler_angle_uv {
	xy_inverse_z_euler_angle_uv(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center,
					const T* const point3D_inverse_z,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T Xc[3];
		Xc[0] = point3D_inverse_z[0] - camera_center[0];
		Xc[1] = point3D_inverse_z[1] - camera_center[1];
		Xc[2] = T(1)/point3D_inverse_z[2] - camera_center[2];

		// P = R * (X - C)
		T p[3];
		p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
		p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
		p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<xy_inverse_z_euler_angle_uv, 2, 3, 3, 3>(
			new xy_inverse_z_euler_angle_uv(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct depth_euler_angle_uv {
	depth_euler_angle_uv(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center,
					const T* const point3D,//direction+depth
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T ptXj[3];//take world system as reference
		ptXj[0] = sin(point3D[0]) * cos(point3D[1]);
		ptXj[1] = sin(point3D[1]);
		ptXj[2] = cos(point3D[0]) * cos(point3D[1]);

		T Xc[3];
		Xc[0] = ptXj[0]*point3D[2] - camera_center[0];
		Xc[1] = ptXj[1]*point3D[2] - camera_center[1];
		Xc[2] = ptXj[2]*point3D[2] - camera_center[2];

		// P = R * (X - C)
		T p[3];
		p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
		p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
		p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<depth_euler_angle_uv, 2, 3, 3, 3>(
			new depth_euler_angle_uv(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct inverse_depth_euler_angle_uv {
	inverse_depth_euler_angle_uv(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center,
					const T* const point3D,//direction+inverse depth
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T ptXj[3];//take world system as reference
		ptXj[0] = sin(point3D[0]) * cos(point3D[1]);
		ptXj[1] = sin(point3D[1]);
		ptXj[2] = cos(point3D[0]) * cos(point3D[1]);

		T Xc[3];
		Xc[0] = ptXj[0]/point3D[2] - camera_center[0];
		Xc[1] = ptXj[1]/point3D[2] - camera_center[1];
		Xc[2] = ptXj[2]/point3D[2] - camera_center[2];

		// P = R * (X - C)
		T p[3];
		p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
		p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
		p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<inverse_depth_euler_angle_uv, 2, 3, 3, 3>(
			new inverse_depth_euler_angle_uv(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	
	//one archor
	struct archored_xyz_euler_angle_uv_nM {
	archored_xyz_euler_angle_uv_nM(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const archored_point3D,//相较于主锚点的笛卡尔坐标
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T p[3];
		p[0] = R[0] * archored_point3D[0] + R[1] * archored_point3D[1] + R[2] * archored_point3D[2];//X
		p[1] = R[3] * archored_point3D[0] + R[4] * archored_point3D[1] + R[5] * archored_point3D[2];//Y
		p[2] = R[6] * archored_point3D[0] + R[7] * archored_point3D[1] + R[8] * archored_point3D[2];//Z

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<archored_xyz_euler_angle_uv_nM, 2, 3, 3>(
			new archored_xyz_euler_angle_uv_nM(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct archored_xyz_euler_angle_uv_nP {
	archored_xyz_euler_angle_uv_nP(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const archored_point3D,//相较于主锚点的笛卡尔坐标
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T nM2nP[3];
		nM2nP[0] = camera_center_nP[0] - camera_center_nM[0];
		nM2nP[1] = camera_center_nP[1] - camera_center_nM[1];
		nM2nP[2] = camera_center_nP[2] - camera_center_nM[2];
		T ray[3];
		ray[0] = archored_point3D[0] - nM2nP[0];
		ray[1] = archored_point3D[1] - nM2nP[1];
		ray[2] = archored_point3D[2] - nM2nP[2];

		T p[3];
		p[0] = R[0] * ray[0] + R[1] * ray[1] + R[2] * ray[2];
		p[1] = R[3] * ray[0] + R[4] * ray[1] + R[5] * ray[2];
		p[2] = R[6] * ray[0] + R[7] * ray[1] + R[8] * ray[2];

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<archored_xyz_euler_angle_uv_nP, 2, 3, 3, 3, 3>(
			new archored_xyz_euler_angle_uv_nP(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct archored_xy_inverse_z_euler_angle_uv_nM {
	archored_xy_inverse_z_euler_angle_uv_nM(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const archored_point3D_idt,//相较于主锚点的笛卡尔坐标
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T p[3];
		p[0] = R[0] * archored_point3D_idt[0] + R[1] * archored_point3D_idt[1] + R[2] / archored_point3D_idt[2];//X
		p[1] = R[3] * archored_point3D_idt[0] + R[4] * archored_point3D_idt[1] + R[5] / archored_point3D_idt[2];//Y
		p[2] = R[6] * archored_point3D_idt[0] + R[7] * archored_point3D_idt[1] + R[8] / archored_point3D_idt[2];//Z

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<archored_xy_inverse_z_euler_angle_uv_nM, 2, 3, 3>(
			new archored_xy_inverse_z_euler_angle_uv_nM(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct archored_xy_inverse_z_euler_angle_uv_nP {
	archored_xy_inverse_z_euler_angle_uv_nP(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const archored_point3D_idt,//相较于主锚点的笛卡尔坐标
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T nM2nP[3];
		nM2nP[0] = camera_center_nP[0] - camera_center_nM[0];
		nM2nP[1] = camera_center_nP[1] - camera_center_nM[1];
		nM2nP[2] = camera_center_nP[2] - camera_center_nM[2];
		T ray[3];
		ray[0] = archored_point3D_idt[0] - nM2nP[0];
		ray[1] = archored_point3D_idt[1] - nM2nP[1];
		ray[2] = T(1)/archored_point3D_idt[2] - nM2nP[2];

		T p[3];
		p[0] = R[0] * ray[0] + R[1] * ray[1] + R[2] * ray[2];
		p[1] = R[3] * ray[0] + R[4] * ray[1] + R[5] * ray[2];
		p[2] = R[6] * ray[0] + R[7] * ray[1] + R[8] * ray[2];

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<archored_xy_inverse_z_euler_angle_uv_nP, 2, 3, 3, 3, 3>(
			new archored_xy_inverse_z_euler_angle_uv_nP(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct archored_depth_euler_angle_uv_nM {
	archored_depth_euler_angle_uv_nM(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const point3D_direction,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T ptXj[3];//主锚点到物点的向量
		ptXj[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXj[1] = sin(point3D_direction[1]);
		ptXj[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		T p[3];
		p[0] = R[0] * ptXj[0] + R[1] * ptXj[1] + R[2] * ptXj[2];//X
		p[1] = R[3] * ptXj[0] + R[4] * ptXj[1] + R[5] * ptXj[2];//Y
		p[2] = R[6] * ptXj[0] + R[7] * ptXj[1] + R[8] * ptXj[2];//Z

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<archored_depth_euler_angle_uv_nM, 2, 3, 2>(
			new archored_depth_euler_angle_uv_nM(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct archored_depth_euler_angle_uv_nP {
	archored_depth_euler_angle_uv_nP(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const point3D_direction,
					const T* const point3D_depth,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T ptXj[3];//主锚点到物点的向量
		ptXj[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]) * point3D_depth[0];
		ptXj[1] = sin(point3D_direction[1]) * point3D_depth[0];
		ptXj[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]) * point3D_depth[0];
		T ptXm[3];//主锚点到当前锚点的向量
		ptXm[0] = camera_center_nP[0] - camera_center_nM[0];
		ptXm[1] = camera_center_nP[1] - camera_center_nM[1];
		ptXm[2] = camera_center_nP[2] - camera_center_nM[2];

		T Xc[3];//当前锚点到物点的向量
		Xc[0] = ptXj[0] - ptXm[0];
		Xc[1] = ptXj[1] - ptXm[1];
		Xc[2] = ptXj[2] - ptXm[2];

		// P = R * (X - C)
		T p[3];
		p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
		p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
		p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<archored_depth_euler_angle_uv_nP, 2, 3, 3, 3, 2, 1>(
			new archored_depth_euler_angle_uv_nP(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct archored_inverse_depth_euler_angle_uv_nM {
	archored_inverse_depth_euler_angle_uv_nM(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const point3D_direction,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T ptXj[3];//主锚点到物点的向量
		ptXj[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXj[1] = sin(point3D_direction[1]);
		ptXj[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		T p[3];
		p[0] = R[0] * ptXj[0] + R[1] * ptXj[1] + R[2] * ptXj[2];//X
		p[1] = R[3] * ptXj[0] + R[4] * ptXj[1] + R[5] * ptXj[2];//Y
		p[2] = R[6] * ptXj[0] + R[7] * ptXj[1] + R[8] * ptXj[2];//Z

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<archored_inverse_depth_euler_angle_uv_nM, 2, 3, 2>(
			new archored_inverse_depth_euler_angle_uv_nM(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct archored_inverse_depth_euler_angle_uv_nP {
	archored_inverse_depth_euler_angle_uv_nP(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const point3D_direction,
					const T* const point3D_inverse_depth,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T ptXj[3];//主锚点到物点的向量
		ptXj[0] = sin(point3D_direction[0]) * cos(point3D_direction[1])/point3D_inverse_depth[0];
		ptXj[1] = sin(point3D_direction[1])/point3D_inverse_depth[0];
		ptXj[2] = cos(point3D_direction[0]) * cos(point3D_direction[1])/point3D_inverse_depth[0];
		T ptXm[3];//主锚点到当前锚点的向量
		ptXm[0] = camera_center_nP[0] - camera_center_nM[0];
		ptXm[1] = camera_center_nP[1] - camera_center_nM[1];
		ptXm[2] = camera_center_nP[2] - camera_center_nM[2];

		T Xc[3];//当前锚点到物点的向量
		Xc[0] = ptXj[0] - ptXm[0];
		Xc[1] = ptXj[1] - ptXm[1];
		Xc[2] = ptXj[2] - ptXm[2];

		// P = R * (X - C)
		T p[3];
		p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
		p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
		p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<archored_inverse_depth_euler_angle_uv_nP, 2, 3, 3, 3, 2, 1>(
			new archored_inverse_depth_euler_angle_uv_nP(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};

	//two archors
	struct parallax_euler_angle_uv_nM {
	parallax_euler_angle_uv_nM(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const point3D_direction,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T p[3];
		T ptXj[3];
		ptXj[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXj[1] = sin(point3D_direction[1]);
		ptXj[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);
		p[0] = R[0] * ptXj[0] + R[1] * ptXj[1] + R[2] * ptXj[2];
		p[1] = R[3] * ptXj[0] + R[4] * ptXj[1] + R[5] * ptXj[2];
		p[2] = R[6] * ptXj[0] + R[7] * ptXj[1] + R[8] * ptXj[2];

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<parallax_euler_angle_uv_nM, 2, 3, 2>(
			new parallax_euler_angle_uv_nM(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct parallax_euler_angle_uv_nA {
	parallax_euler_angle_uv_nA(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const point3D_direction,
					const T* const point3D_parallax,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

		T p[3], pti2k[3], ptXUnit[3], ptXk[3];
		//主锚点到副锚点的平移向量
		pti2k[0] = camera_center_nP[0] - camera_center_nM[0];	
		pti2k[1] = camera_center_nP[1] - camera_center_nM[1];	
		pti2k[2] = camera_center_nP[2] - camera_center_nM[2];	

		//主锚点到特征点的单位向量
		ptXUnit[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXUnit[1] = sin(point3D_direction[1]);
		ptXUnit[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		//compute angle w2
		T dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1] + ptXUnit[2]*pti2k[2];
		T dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
		T dW2;
		if (dDot/dDisi2k > T(1))
			dW2 = T(0);
		if (dDot/dDisi2k < T(-1))
			dW2 = T(PI);
		else
			dW2  = acos(dDot/dDisi2k);

		//compute Xk vector according sin theory
		ptXk[0] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[0] - sin(point3D_parallax[0]) * pti2k[0];
		ptXk[1] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[1] - sin(point3D_parallax[0]) * pti2k[1];
		ptXk[2] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[2] - sin(point3D_parallax[0]) * pti2k[2];

		p[0] = R[0] * ptXk[0] + R[1] * ptXk[1] + R[2] * ptXk[2];
		p[1] = R[3] * ptXk[0] + R[4] * ptXk[1] + R[5] * ptXk[2];
		p[2] = R[6] * ptXk[0] + R[7] * ptXk[1] + R[8] * ptXk[2];

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);
		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {	
		// printf("%s\n","ReprojectionError entered");	
		return new ceres::AutoDiffCostFunction<parallax_euler_angle_uv_nA, 2, 3, 3, 3, 2, 1>(
			new parallax_euler_angle_uv_nA(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct parallax_euler_angle_uv_nP {
	parallax_euler_angle_uv_nP(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const euler_angles,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const camera_center_nA,
					const T* const point3D_direction,
					const T* const point3D_parallax,
					T* residuals) const {
		T R[9];
		T ey = euler_angles[0];
		T ex = euler_angles[1];
		T ez = euler_angles[2];
		T c1 = cos(ey);   T c2 = cos(ex);   T c3 = cos(ez);
		T s1 = sin(ey);   T s2 = sin(ex);   T s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;
	
		T p[3], pti2k[3], pti2l[3], ptXUnit[3], ptXk[3];

		//主锚点到副锚点的平移向量
		pti2k[0] = camera_center_nA[0] - camera_center_nM[0];		
		pti2k[1] = camera_center_nA[1] - camera_center_nM[1];		
		pti2k[2] = camera_center_nA[2] - camera_center_nM[2];
		//主锚点到当且锚点的平移向量
		pti2l[0] = camera_center_nP[0] - camera_center_nM[0];	
		pti2l[1] = camera_center_nP[1] - camera_center_nM[1];		
		pti2l[2] = camera_center_nP[2] - camera_center_nM[2];
		
		//XUnit 主锚点到特征点的单位向量
		ptXUnit[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXUnit[1] = sin(point3D_direction[1]);
		ptXUnit[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		//compute angle w2
		T dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1]+ ptXUnit[2]*pti2k[2];
		T dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
		T dW2;
		//dW2  = acos( dDot/dDisi2k );
		if (dDot/dDisi2k > T(1))
			dW2 = T(0);
		if ( dDot/dDisi2k < T(-1))
			dW2 = T(PI);
		else
			dW2  = acos(dDot/dDisi2k);

		//compute Xl vector according sin theory
		ptXk[0] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[0] - sin(point3D_parallax[0]) * pti2l[0];
		ptXk[1] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[1] - sin(point3D_parallax[0]) * pti2l[1];
		ptXk[2] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[2] - sin(point3D_parallax[0]) * pti2l[2];
		
		p[0] = R[0] * ptXk[0] + R[1] * ptXk[1] + R[2] * ptXk[2];
		p[1] = R[3] * ptXk[0] + R[4] * ptXk[1] + R[5] * ptXk[2];
		p[2] = R[6] * ptXk[0] + R[7] * ptXk[1] + R[8] * ptXk[2];

		// Normalize
		T xp = p[0] / p[2];
		T yp = p[1] / p[2];

		// Project
		T predicted_u = fx * xp + cx;
		T predicted_v = fy * yp + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);
		
		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {	
		// printf("%s\n","ReprojectionError entered");	
		return new ceres::AutoDiffCostFunction<parallax_euler_angle_uv_nP, 2, 3, 3, 3, 3, 2, 1>(
			new parallax_euler_angle_uv_nP(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	

	//Parameterization of 3d rotation
	struct axis_angle_xyz_uv {
	axis_angle_xyz_uv(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const axis_angle,
					const T* const camera_center,
					const T* const point3D,
					T* residuals) const {

		T Xc[3];
		Xc[0] = point3D[0] - camera_center[0];
		Xc[1] = point3D[1] - camera_center[1];
		Xc[2] = point3D[2] - camera_center[2];

		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::AngleAxisRotatePoint(axis_angle, Xc, transformed_point.data());
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<axis_angle_xyz_uv, 2, 3, 3, 3>(
			new axis_angle_xyz_uv(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct quaternion_xyz_uv {
	quaternion_xyz_uv(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const quat,
					const T* const camera_center,
					const T* const point3D,
					T* residuals) const {

		T Xc[3];
		Xc[0] = point3D[0] - camera_center[0];
		Xc[1] = point3D[1] - camera_center[1];
		Xc[2] = point3D[2] - camera_center[2];

		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::QuaternionRotatePoint(quat, Xc, transformed_point.data());
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<quaternion_xyz_uv, 2, 4, 3, 3>(
			new quaternion_xyz_uv(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	//xyz_euler_angle_uv

	struct axis_angle_parallax_uv_nM {
	axis_angle_parallax_uv_nM(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const axis_angle,
					const T* const point3D_direction,
					T* residuals) const {

		T ptXj[3];
		ptXj[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXj[1] = sin(point3D_direction[1]);
		ptXj[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::AngleAxisRotatePoint(axis_angle, ptXj, transformed_point.data());
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<axis_angle_parallax_uv_nM, 2, 3, 2>(
			new axis_angle_parallax_uv_nM(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct axis_angle_parallax_uv_nA {
	axis_angle_parallax_uv_nA(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const axis_angle,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const point3D_direction,
					const T* const point3D_parallax,
					T* residuals) const {


		T p[3], pti2k[3], ptXUnit[3], ptXk[3];
		//主锚点到副锚点的平移向量
		pti2k[0] = camera_center_nP[0] - camera_center_nM[0];	
		pti2k[1] = camera_center_nP[1] - camera_center_nM[1];	
		pti2k[2] = camera_center_nP[2] - camera_center_nM[2];	

		//主锚点到特征点的单位向量
		ptXUnit[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXUnit[1] = sin(point3D_direction[1]);
		ptXUnit[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		//compute angle w2
		T dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1] + ptXUnit[2]*pti2k[2];
		T dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
		T dW2;
		if (dDot/dDisi2k > T(1))
			dW2 = T(0);
		if (dDot/dDisi2k < T(-1))
			dW2 = T(PI);
		else
			dW2  = acos(dDot/dDisi2k);

		//compute Xk vector according sin theory
		ptXk[0] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[0] - sin(point3D_parallax[0]) * pti2k[0];
		ptXk[1] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[1] - sin(point3D_parallax[0]) * pti2k[1];
		ptXk[2] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[2] - sin(point3D_parallax[0]) * pti2k[2];

		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::AngleAxisRotatePoint(axis_angle, ptXk, transformed_point.data());
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();
		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;
		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);
		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {	
		// printf("%s\n","ReprojectionError entered");	
		return new ceres::AutoDiffCostFunction<axis_angle_parallax_uv_nA, 2, 3, 3, 3, 2, 1>(
			new axis_angle_parallax_uv_nA(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct axis_angle_parallax_uv_nP {
	axis_angle_parallax_uv_nP(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const axis_angle,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const camera_center_nA,
					const T* const point3D_direction,
					const T* const point3D_parallax,
					T* residuals) const {
	
		T pti2k[3], pti2l[3], ptXUnit[3], ptXk[3];

		//主锚点到副锚点的平移向量
		pti2k[0] = camera_center_nA[0] - camera_center_nM[0];		
		pti2k[1] = camera_center_nA[1] - camera_center_nM[1];		
		pti2k[2] = camera_center_nA[2] - camera_center_nM[2];
		//主锚点到当且锚点的平移向量
		pti2l[0] = camera_center_nP[0] - camera_center_nM[0];	
		pti2l[1] = camera_center_nP[1] - camera_center_nM[1];		
		pti2l[2] = camera_center_nP[2] - camera_center_nM[2];
		
		//XUnit 主锚点到特征点的单位向量
		ptXUnit[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXUnit[1] = sin(point3D_direction[1]);
		ptXUnit[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		//compute angle w2
		T dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1]+ ptXUnit[2]*pti2k[2];
		T dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
		T dW2;
		//dW2  = acos( dDot/dDisi2k );
		if (dDot/dDisi2k > T(1))
			dW2 = T(0);
		if ( dDot/dDisi2k < T(-1))
			dW2 = T(PI);
		else
			dW2  = acos(dDot/dDisi2k);

		//compute Xl vector according sin theory
		ptXk[0] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[0] - sin(point3D_parallax[0]) * pti2l[0];
		ptXk[1] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[1] - sin(point3D_parallax[0]) * pti2l[1];
		ptXk[2] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[2] - sin(point3D_parallax[0]) * pti2l[2];
		
		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::AngleAxisRotatePoint(axis_angle, ptXk, transformed_point.data());
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);
		
		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {	
		// printf("%s\n","ReprojectionError entered");	
		return new ceres::AutoDiffCostFunction<axis_angle_parallax_uv_nP, 2, 3, 3, 3, 3, 2, 1>(
			new axis_angle_parallax_uv_nP(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct quaternion_parallax_uv_nM {
	quaternion_parallax_uv_nM(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const quat,
					const T* const point3D_direction,
					T* residuals) const {

		T ptXj[3];
		ptXj[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXj[1] = sin(point3D_direction[1]);
		ptXj[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::QuaternionRotatePoint(quat, ptXj, transformed_point.data());
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<quaternion_parallax_uv_nM, 2, 4, 2>(
			new quaternion_parallax_uv_nM(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct quaternion_parallax_uv_nA {
	quaternion_parallax_uv_nA(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const quat,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const point3D_direction,
					const T* const point3D_parallax,
					T* residuals) const {


		T p[3], pti2k[3], ptXUnit[3], ptXk[3];
		//主锚点到副锚点的平移向量
		pti2k[0] = camera_center_nP[0] - camera_center_nM[0];	
		pti2k[1] = camera_center_nP[1] - camera_center_nM[1];	
		pti2k[2] = camera_center_nP[2] - camera_center_nM[2];	

		//主锚点到特征点的单位向量
		ptXUnit[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXUnit[1] = sin(point3D_direction[1]);
		ptXUnit[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		//compute angle w2
		T dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1] + ptXUnit[2]*pti2k[2];
		T dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
		T dW2;
		if (dDot/dDisi2k > T(1))
			dW2 = T(0);
		if (dDot/dDisi2k < T(-1))
			dW2 = T(PI);
		else
			dW2  = acos(dDot/dDisi2k);

		//compute Xk vector according sin theory
		ptXk[0] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[0] - sin(point3D_parallax[0]) * pti2k[0];
		ptXk[1] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[1] - sin(point3D_parallax[0]) * pti2k[1];
		ptXk[2] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[2] - sin(point3D_parallax[0]) * pti2k[2];

		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::QuaternionRotatePoint(quat, ptXk, transformed_point.data());
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;
		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);
		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {	
		// printf("%s\n","ReprojectionError entered");	
		return new ceres::AutoDiffCostFunction<quaternion_parallax_uv_nA, 2, 4, 3, 3, 2, 1>(
			new quaternion_parallax_uv_nA(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct quaternion_parallax_uv_nP {
	quaternion_parallax_uv_nP(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const quat,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const camera_center_nA,
					const T* const point3D_direction,
					const T* const point3D_parallax,
					T* residuals) const {
	
		T pti2k[3], pti2l[3], ptXUnit[3], ptXk[3];

		//主锚点到副锚点的平移向量
		pti2k[0] = camera_center_nA[0] - camera_center_nM[0];		
		pti2k[1] = camera_center_nA[1] - camera_center_nM[1];		
		pti2k[2] = camera_center_nA[2] - camera_center_nM[2];
		//主锚点到当且锚点的平移向量
		pti2l[0] = camera_center_nP[0] - camera_center_nM[0];	
		pti2l[1] = camera_center_nP[1] - camera_center_nM[1];		
		pti2l[2] = camera_center_nP[2] - camera_center_nM[2];
		
		//XUnit 主锚点到特征点的单位向量
		ptXUnit[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXUnit[1] = sin(point3D_direction[1]);
		ptXUnit[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		//compute angle w2
		T dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1]+ ptXUnit[2]*pti2k[2];
		T dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
		T dW2;
		//dW2  = acos( dDot/dDisi2k );
		if (dDot/dDisi2k > T(1))
			dW2 = T(0);
		if ( dDot/dDisi2k < T(-1))
			dW2 = T(PI);
		else
			dW2  = acos(dDot/dDisi2k);

		//compute Xl vector according sin theory
		ptXk[0] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[0] - sin(point3D_parallax[0]) * pti2l[0];
		ptXk[1] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[1] - sin(point3D_parallax[0]) * pti2l[1];
		ptXk[2] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[2] - sin(point3D_parallax[0]) * pti2l[2];
		
		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::QuaternionRotatePoint(quat, ptXk, transformed_point.data());
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);
		
		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {	
		// printf("%s\n","ReprojectionError entered");	
		return new ceres::AutoDiffCostFunction<quaternion_parallax_uv_nP, 2, 4, 3, 3, 3, 2, 1>(
			new quaternion_parallax_uv_nP(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};


	//Manifold
	struct xyz_uv_lie_manifold {
	xyz_uv_lie_manifold(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const se3,
		            // const T* const quat,
					// const T* const translation,
					const T* const point3D,
					T* residuals) const {
		
		T quat[4];
		quat[0]=se3[0];
		quat[1]=se3[1];
		quat[2]=se3[2];
		quat[3]=se3[3];
		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::QuaternionRotatePoint(quat, point3D, transformed_point.data());
		Eigen::Map<const Eigen::Matrix<T,3,1>> t(se3+4);
		// Eigen::Map<const Eigen::Matrix<T,3,1>> t(translation);
		transformed_point += t;
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);
		//将李代数se(3)转换为李群SE(3)
		// Eigen::Map<const Eigen::Matrix<T, 7, 1>> se3_vec(se3);
		// Sophus::SE3<T> Tcw = Sophus::SE3<T>::exp(se3_vec);
		// //变换三维点到相机坐标系
		// Eigen::Map<const Eigen::Matrix<T,3,1>> Pw(point3D);
        // Eigen::Matrix<T,3,1> Pc = Tcw.rotationMatrix() * Pw + se3_vec.head<3>();
		// // Eigen::Matrix<T,3,1> Pc = Tcw * Pw;

		// T x = Pc[0] / Pc[2];
		// T y = Pc[1] / Pc[2];
		// // Project
		// T predicted_u = fx * x + cx;
		// T predicted_v = fy * y + cy;

		// // Residual
		// residuals[0] = predicted_u - T(u);
		// residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<xyz_uv_lie_manifold, 2, 7, 3>(
			new xyz_uv_lie_manifold(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct parallax_uv_quaternion_manifold_nM {
	parallax_uv_quaternion_manifold_nM(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const quat,
					const T* const point3D_direction,
					T* residuals) const {

		T ptXj[3];
		ptXj[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXj[1] = sin(point3D_direction[1]);
		ptXj[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::QuaternionRotatePoint(quat, ptXj, transformed_point.data());
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);

		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {		
		return new ceres::AutoDiffCostFunction<parallax_uv_quaternion_manifold_nM, 2, 4, 2>(
			new parallax_uv_quaternion_manifold_nM(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct parallax_uv_quaternion_manifold_nA {
	parallax_uv_quaternion_manifold_nA(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	template <typename T>
	bool operator()(const T* const quat,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const point3D_direction,
					const T* const point3D_parallax,
					T* residuals) const {


		T p[3], pti2k[3], ptXUnit[3], ptXk[3];
		//主锚点到副锚点的平移向量
		pti2k[0] = camera_center_nP[0] - camera_center_nM[0];	
		pti2k[1] = camera_center_nP[1] - camera_center_nM[1];	
		pti2k[2] = camera_center_nP[2] - camera_center_nM[2];	

		//主锚点到特征点的单位向量
		ptXUnit[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXUnit[1] = sin(point3D_direction[1]);
		ptXUnit[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		//compute angle w2
		T dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1] + ptXUnit[2]*pti2k[2];
		T dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
		T dW2;
		if (dDot/dDisi2k > T(1))
			dW2 = T(0);
		if (dDot/dDisi2k < T(-1))
			dW2 = T(PI);
		else
			dW2  = acos(dDot/dDisi2k);

		//compute Xk vector according sin theory
		ptXk[0] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[0] - sin(point3D_parallax[0]) * pti2k[0];
		ptXk[1] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[1] - sin(point3D_parallax[0]) * pti2k[1];
		ptXk[2] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[2] - sin(point3D_parallax[0]) * pti2k[2];

		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::QuaternionRotatePoint(quat, ptXk, transformed_point.data());
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;
		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);
		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {	
		// printf("%s\n","ReprojectionError entered");	
		return new ceres::AutoDiffCostFunction<parallax_uv_quaternion_manifold_nA, 2, 4, 3, 3, 2, 1>(
			new parallax_uv_quaternion_manifold_nA(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
	struct parallax_uv_quaternion_manifold_nP {
	parallax_uv_quaternion_manifold_nP(double observed_u, double observed_v,
					double fx, double fy, double cx, double cy)
		: u(observed_u), v(observed_v),
		fx(fx), fy(fy), cx(cx), cy(cy) {}

	// mutable bool printed_once_ = false;  // 每个 CostFunctor 对象内部标志
	template <typename T>
	bool operator()(const T* const quat,
					const T* const camera_center_nP,
					const T* const camera_center_nM,
					const T* const camera_center_nA,
					const T* const point3D_direction,
					const T* const point3D_parallax,
					T* residuals) const {
	
		T pti2k[3], pti2l[3], ptXUnit[3], ptXk[3];

		//主锚点到副锚点的平移向量
		pti2k[0] = camera_center_nA[0] - camera_center_nM[0];		
		pti2k[1] = camera_center_nA[1] - camera_center_nM[1];		
		pti2k[2] = camera_center_nA[2] - camera_center_nM[2];
		//主锚点到当且锚点的平移向量
		pti2l[0] = camera_center_nP[0] - camera_center_nM[0];	
		pti2l[1] = camera_center_nP[1] - camera_center_nM[1];		
		pti2l[2] = camera_center_nP[2] - camera_center_nM[2];
		
		//XUnit 主锚点到特征点的单位向量
		ptXUnit[0] = sin(point3D_direction[0]) * cos(point3D_direction[1]);
		ptXUnit[1] = sin(point3D_direction[1]);
		ptXUnit[2] = cos(point3D_direction[0]) * cos(point3D_direction[1]);

		//compute angle w2
		T dDot = ptXUnit[0]*pti2k[0] + ptXUnit[1]*pti2k[1]+ ptXUnit[2]*pti2k[2];
		T dDisi2k = sqrt(pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2]);
		T dW2;
		//dW2  = acos( dDot/dDisi2k );
		if (dDot/dDisi2k > T(1))
			dW2 = T(0);
		if ( dDot/dDisi2k < T(-1))
			dW2 = T(PI);
		else
			dW2  = acos(dDot/dDisi2k);

		//compute Xl vector according sin theory
		ptXk[0] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[0] - sin(point3D_parallax[0]) * pti2l[0];
		ptXk[1] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[1] - sin(point3D_parallax[0]) * pti2l[1];
		ptXk[2] = dDisi2k * sin(dW2+point3D_parallax[0]) * ptXUnit[2] - sin(point3D_parallax[0]) * pti2l[2];
		
		Eigen::Matrix<T, 3, 1> transformed_point;
		ceres::QuaternionRotatePoint(quat, ptXk, transformed_point.data());
		const Eigen::Matrix<T, 2, 1> projected_point = transformed_point.hnormalized();

		// Project
		T predicted_u = fx * projected_point.x() + cx;
		T predicted_v = fy * projected_point.y() + cy;

		// Residual
		residuals[0] = predicted_u - T(u);
		residuals[1] = predicted_v - T(v);
		
		return true;
	}

	static ceres::CostFunction* Create(double u, double v,
									double fx, double fy, double cx, double cy) {	
		// printf("%s\n","ReprojectionError entered");	
		return new ceres::AutoDiffCostFunction<parallax_uv_quaternion_manifold_nP, 2, 4, 3, 3, 3, 2, 1>(
			new parallax_uv_quaternion_manifold_nP(u, v, fx, fy, cx, cy));
	}

	double u, v;
	double fx, fy, cx, cy;
	};
};