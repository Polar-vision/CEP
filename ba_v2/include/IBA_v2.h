#pragma once
#include "BAExporter_v2.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <time.h>

#include <Eigen/Core>
// #include <Eigen/Geometry>
// #include <Eigen/LU>
// #include <Eigen/StdVector>
// #include <Eigen/Cholesky>

#include "ceres/ceres.h"
#include "ceres/rotation.h"
#include "ceres/manifold.h"
// #include "ceres/local_parameterization.h"
#include "sophus/se3.hpp"

#define PI  3.1415926535898 
#define MAXARCHOR 0.5
using namespace std;
using namespace Eigen;
#define MAXSTRLEN  2048 /* 2K */
#define SKIP_LINE(f){                                                       \
	char buf[MAXSTRLEN];                                                        \
	while(!feof(f))                                                           \
	if(!fgets(buf, MAXSTRLEN-1, f) || buf[strlen(buf)-1]=='\n') break;      \
}

class IBA
{
public:
    IBA(void);
    ~IBA(void);

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
    					manifoldtype manitype)=0;

	virtual bool ba_initialize(char* szCamera, char* szFeature, char* szCalib = NULL, char* szXYZ = NULL)=0;

	int findNcameras(FILE* fp);
	void ba_readCameraPoseration(char* fname, double* ical);
	void readNpointsAndNprojections(FILE* fp, int* n3Dpts, int pnp, int* nprojs, int mnp);
	void ba_readCameraPose(FILE* fp, double* params, int* m_v);
	void ba_updateKR(double* KR, double* K, double* p);
	int readNInts(FILE* fp, int* vals, int nvals);
	int readNDoubles(FILE* fp, double* vals, int nvals);
	int countNDoubles(FILE* fp);
	int skipNDoubles(FILE* fp, int nvals);
	void ba_readCablibration(FILE* fp, double* K);
	void readNpointsAndNprojectionsFromProj(FILE* fp, int& n3Dpts, int& nprojs);
	void readPointProjections(FILE* fp, double* imgpts, int* photo, int* imgptsSum, int n3Dpts, int n2Dprojs);
	void readImagePts(const char* szProj, double** imgpts, int** photo, int** imgptsSum, int& n3Dpts, int& n2Dprojs);
	int		m_ncams, m_n3Dpts, m_n2Dprojs, m_nS, nc_;  //number of camera, 3D points, 2D projection points, non-zero element of S matrix
	int* m_archor;
	int* m_photo, * m_feature;
	double* m_motstruct, * m_imgpts;			  //6 camera pose and 3 feature parameters/PBA parameter,
	double* m_XYZ;								  //initial XYZ provided 	
	double* m_K;								  //calibration parameters
	//std::vector<double> m_K;
	int* m_V;//camera id
	char* m_vmask, * m_umask, * m_smask;
	int* m_imgptsSum, * m_struct, * m_pnt2main, * m_archorSort;
	double* m_KR, * m_KdA, * m_KdB, * m_KdG, * m_P;
	bool    m_bProvideXYZ, m_bFocal;
	char* m_szCameraInit;
	char* m_szFeatures;
	char* m_szCalibration;
	char* m_szXYZ;
	char* m_szCamePose;
	char* m_sz3Dpts;
	char* m_szReport;

	struct Intrinsic{
		double fx, fy, cx, cy;
	};
	struct Camera{
		double euler_angle[3];
		double axis_angle[3];
		double quat[4];
		double camera_center[3];
		double translation[3];
		double se3[7];
		int camidx;
	};
	struct Point3D{
		double xyz[3];
		double xy_inverse_z[3];
		double archored_xyz[3];//take main archor as reference
		double archored_xy_inverse_z[3];//take main archor as reference
		double inverse_depth;//take main archor as reference
		double depth;//take main archor as reference
		double direction[2];//take main archor as reference

		double world_depth[3];
		double world_inverse_depth[3];

		double parallax;//take main archor and associate archor as references

		int nM;//主锚点
		int nA;//副锚点
	};
	struct Observation{
		int view_idx;
		double u, v;
		double lightcone;
		double polar;
	};
	struct Track{
		int nview;
		std::vector<Observation> obss;
	};
	struct FObservation{
		int feature_idx;
		double u,v;
		double lightcone;
		double polar;
	};
	// struct ImageTrack{
	// 	int idx;
	// 	std::vector<FObservation> fobs;
	// };
	struct Ray{
		double ray[3];
	};
	std::vector<Intrinsic> intrs;
	std::vector<Camera> cams;
	std::vector<Point3D> points;
	std::vector<Track> tracks;
	std::vector<Ray> rays;
	// std::unordered_map<int, std::vector<FObservation>> image_tracks;
	std::map<int, std::vector<FObservation>> image_tracks;
};
