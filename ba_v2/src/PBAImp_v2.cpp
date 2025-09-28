#include "PBAImp_v2.h"

PBA::PBA(void)
{
	m_bProvideXYZ = false;
	m_bFocal = false;
	m_szCameraInit = m_szFeatures = m_szCalibration = m_szXYZ = m_sz3Dpts = m_szCamePose = m_szReport = NULL;
}


PBA::~PBA(void)
{
	if ( m_szCameraInit!= NULL)	free(m_szCameraInit);
	if ( m_szFeatures!= NULL)	free(m_szFeatures);
	if ( m_szCalibration!= NULL)	free(m_szCalibration);
	if ( m_szXYZ!= NULL)		free(m_szXYZ);
	if ( m_szCamePose!= NULL)	free(m_szCamePose);
	if ( m_sz3Dpts!= NULL)		free(m_sz3Dpts);
	if ( m_szReport!= NULL)		free(m_szReport);
}

double light_cone_obs(double u, double v, 
                                  double fx, double fy, 
                                  double cx, double cy) {
    Eigen::Vector3d ray(
        (u - cx) / fx,
        (v - cy) / fy,
        1.0
    );
    
    ray.normalize();
    
    Eigen::Vector3d optical_axis(0.0, 0.0, 1.0);
    
    double cos_theta = optical_axis.dot(ray);
    cos_theta = std::clamp(cos_theta, -1.0, 1.0);
    
    return std::acos(cos_theta);
}


double PBA::lc_rep(){
	double s = 0.0;
	int nobs = 0;
	for (int i = 0; i < tracks.size(); i++)
	{ 
		nobs += tracks[i].nview;
		int nM = points[i].nM;
		int nA = points[i].nA;

		for (int j = 0; j < tracks[i].nview; j++)
		{
			int view_idx = tracks[i].obss[j].view_idx;
			double u = tracks[i].obss[j].u;
			double v = tracks[i].obss[j].v;
			double fx = 0, fy = 0, cx = 0, cy = 0;
			int cam_idx = cams[view_idx].camidx;
			fx = intrs[cam_idx-1].fx;
			fy = intrs[cam_idx-1].fy;
			cx = intrs[cam_idx-1].cx;
			cy = intrs[cam_idx-1].cy;
	
			int nP = view_idx;
			double lco = atan(sqrt(((u - cx) / fx) * ((u - cx) / fx) + ((v - cy) / fy) * ((v - cy) / fy)));

			double R[9];
			double ey = cams[nP].euler_angle[0];
			double ex = cams[nP].euler_angle[1];
			double ez = cams[nP].euler_angle[2];
			double c1 = cos(ey), c2 = cos(ex), c3 = cos(ez);
			double s1 = sin(ey), s2 = sin(ex), s3 = sin(ez);
			R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
			R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
			R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;

			double Xc[3];
			Xc[0] = points[i].xyz[0] - cams[nP].camera_center[0];
			Xc[1] = points[i].xyz[1] - cams[nP].camera_center[1];
			Xc[2] = points[i].xyz[2] - cams[nP].camera_center[2];

			double p[3];
			p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];
			p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];
			p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];

			double lcp = atan(sqrt(p[0] * p[0] + p[1] * p[1]) / p[2]);

			// Residual
			double res2 = (lcp - lco)*(lcp - lco);
			s += res2;
		}
	}
	s /= (2*nobs);
	return sqrt(s);
}

struct CostRecorderCallback_simple : public ceres::IterationCallback {
    std::vector<double>& costs;

    CostRecorderCallback_simple(std::vector<double>& costs_) : costs(costs_) {}

    virtual ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override {
        costs.push_back(summary.cost);  
        return ceres::SOLVER_CONTINUE;
    }
};
struct CostRecorderCallback : public ceres::IterationCallback {
    std::vector<double>& costs;
    std::vector<std::vector<double>>& euler_history; // 每次迭代存全部相机欧拉角
    std::vector<PBA::Camera>& cameras; // 指向相机数组
    int num_cams;    // 相机数量

    CostRecorderCallback(std::vector<double>& costs_,
                         std::vector<std::vector<double>>& euler_history_,
                         std::vector<PBA::Camera>& cams, int num_cams_)
        : costs(costs_), euler_history(euler_history_),
          cameras(cams), num_cams(num_cams_) {}

    virtual ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override {
        costs.push_back(summary.cost);
        // 保存本轮迭代所有相机的欧拉角
        std::vector<double> this_iter;
        this_iter.reserve(num_cams * 6);
        for (int i = 0; i < num_cams; i++) {
            this_iter.push_back(cameras[i].euler_angle[0]);
            this_iter.push_back(cameras[i].euler_angle[1]);
            this_iter.push_back(cameras[i].euler_angle[2]);
			this_iter.push_back(cameras[i].camera_center[0]);
			this_iter.push_back(cameras[i].camera_center[1]);
			this_iter.push_back(cameras[i].camera_center[2]);
        }
        euler_history.push_back(this_iter);

        return ceres::SOLVER_CONTINUE;
    }
};

void PBA::narrow_test(){
	ceres::Problem problem;
	ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);
	//1.6 degree test
	Camera cam;
	cam.camera_center[0]=0;
	cam.camera_center[1]=0;
	cam.camera_center[2]=0;
	cam.euler_angle[0]=0;
	cam.euler_angle[1]=0;
	cam.euler_angle[2]=0;
	cams.push_back(cam);
	cam.camera_center[0]=10;
	cam.camera_center[1]=0;
	cam.camera_center[2]=0;
	cam.euler_angle[0]=0;
	cam.euler_angle[1]=0;
	cam.euler_angle[2]=0;
	cams.push_back(cam);
	Point3D point;
	point.xyz[0]=cams[1].camera_center[0]/2;
	point.xyz[1]=0;
	double ray1[3],ray2[3];
	ray1[0] = point.xyz[0]-cams[0].camera_center[0];
	ray1[1] = point.xyz[1]-cams[0].camera_center[1];
	ray2[0] = point.xyz[0]-cams[1].camera_center[0];
	ray2[1] = point.xyz[1]-cams[1].camera_center[1];
	// for(double i=10000;i>100;i/=1.1){
		point.xyz[2]=355;
		point.nM=0;
		point.nA=1;
		ray1[2] = point.xyz[2]-cams[0].camera_center[2];
		ray2[2] = point.xyz[2]-cams[1].camera_center[2];

		double norm_ray1 = sqrt(ray1[0]*ray1[0]+ray1[1]*ray1[1]+ray1[2]*ray1[2]);
		double norm_ray2 = sqrt(ray2[0]*ray2[0]+ray2[1]*ray2[1]+ray2[2]*ray2[2]);

		double dot_ray = ray1[0]*ray2[0]+ray1[1]*ray2[1]+ray1[2]*ray2[2];
		point.parallax = acos(dot_ray/(norm_ray1*norm_ray2));
		point.direction[0] = atan2( ray1[0], ray1[2] );				
		point.direction[1] = atan2( ray1[1], sqrt(ray1[0]*ray1[0]+ ray1[2]*ray1[2]));
		points.push_back(point);
		// double intersect_angle = acos(dot_ray/(norm_ray1*norm_ray2))*180/PI;
		// printf("%s%f%s%f%s\n","depth:",i,"  intersection angle:",intersect_angle,"degree");
		printf("%s%f%s\n","intersection angle:",point.parallax*180/PI,"degree");
	// }
	// Normalize
	double xp1 = ray1[0] / ray1[2], xp2 = ray2[0] / ray2[2];
	double yp1 = ray1[1] / ray1[2], yp2 = ray2[1] / ray2[2];

	// Project
	double fx = 4209.003931;
	double fy = 4220.634348;
	double cx = 2736;
	double cy = 1824;
	Track track;
	track.nview = 2;

	for(int i=0;i<track.nview;i++){
		Observation obs;
		obs.view_idx = i;
		obs.u = fx * xp1 + cx;
		obs.v = fy * yp1 + cy;
		track.obss.push_back(obs);
		printf("%s%f%s%f%s\n","(u,v)=(",obs.u,",",obs.v,")");
	}
	tracks.push_back(track);
	for(int i=0;i<cams.size();i++){
		printf("%s%f %f %f %f %f %f%s\n","(phi,omega,kappa,x,y,z)=(",cams[i].euler_angle[0],cams[i].euler_angle[1],cams[i].euler_angle[2],
			cams[i].camera_center[0],cams[i].camera_center[1],cams[i].camera_center[2],")");
	}
	for(int i=0;i<points.size();i++){
		printf("%s%f %f %f %f %f %f%s\n","(X,Y,Z,alpha,beta,gamma)=(",points[i].xyz[0],points[i].xyz[1],points[i].xyz[2],
			points[i].direction[0]*180/PI,points[i].direction[1]*180/PI,points[i].parallax*180/PI,")");
	}
	FILE* Gfp1w = NULL,* Gfp2w = NULL;
	Gfp1w = fopen("E:/zuo/projects/backup/data/narrow_test/G-XYZ.txt", "w");
	for (const auto& it : points) {
		fprintf(Gfp1w, "%f %f %f\n", it.xyz[0], it.xyz[1], it.xyz[2]);
	}
	Gfp2w = fopen("E:/zuo/projects/backup/data/narrow_test/G-Cam.txt", "w");
	for(const auto& it:cams){
		fprintf(Gfp2w, "%f %f %f %f %f %f\n", it.euler_angle[0], it.euler_angle[1], it.euler_angle[2],
		                                     it.camera_center[0],it.camera_center[1],it.camera_center[2]);
	}
	fclose(Gfp1w);
	fclose(Gfp2w);

	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> dist_uv(0, 2);
	normal_distribution<> dist_euler(0, 4);
	normal_distribution<> dist_xyz(0, 10);
	normal_distribution<> dist_XYZ(0, 20);

	for(int i=0;i<tracks.size();i++){
		int _nview=tracks[i].nview;
		int nM=points[i].nM;
		int nA=points[i].nA;
		points[i].xyz[0]+=dist_XYZ(gen);
		points[i].xyz[1]+=dist_XYZ(gen);
		points[i].xyz[2]+=dist_XYZ(gen);

		for(int j=0;j<_nview;j++){
			int nP=tracks[i].obss[j].view_idx;
			tracks[i].obss[j].u+=dist_uv(gen);
			tracks[i].obss[j].v+=dist_uv(gen);
			double u=tracks[i].obss[j].u;
			double v=tracks[i].obss[j].v;
			cams[nP].camera_center[0]+=dist_xyz(gen);
			cams[nP].camera_center[1]+=dist_xyz(gen);
			cams[nP].camera_center[2]+=dist_xyz(gen);
			cams[nP].euler_angle[0]+=dist_euler(gen)*PI/180;
			cams[nP].euler_angle[1]+=dist_euler(gen)*PI/180;
			cams[nP].euler_angle[2]+=dist_euler(gen)*PI/180;				

		}
		ray1[0] = points[i].xyz[0]-cams[nM].camera_center[0];
		ray1[1] = points[i].xyz[1]-cams[nM].camera_center[1];
		ray1[2] = points[i].xyz[2]-cams[nM].camera_center[2];
		ray2[0] = points[i].xyz[0]-cams[nA].camera_center[0];
		ray2[1] = points[i].xyz[1]-cams[nA].camera_center[1];
		ray2[2] = points[i].xyz[2]-cams[nA].camera_center[2];

		norm_ray1 = sqrt(ray1[0]*ray1[0]+ray1[1]*ray1[1]+ray1[2]*ray1[2]);
		norm_ray2 = sqrt(ray2[0]*ray2[0]+ray2[1]*ray2[1]+ray2[2]*ray2[2]);

		dot_ray = ray1[0]*ray2[0]+ray1[1]*ray2[1]+ray1[2]*ray2[2];
		points[i].parallax = acos(dot_ray/(norm_ray1*norm_ray2));
		points[i].direction[0] = atan2( ray1[0], ray1[2] );				
		points[i].direction[1] = atan2( ray1[1], sqrt(ray1[0]*ray1[0]+ ray1[2]*ray1[2]));
	}

	for(int i=0;i<tracks.size();i++){
		for(int j=0;j<tracks[i].nview;j++){
			printf("%s%f%s%f%s\n","(u,v)=(",tracks[i].obss[j].u,",",tracks[i].obss[j].v,")");
		}
	}
	for(int i=0;i<cams.size();i++){
		printf("%s%f %f %f %f %f %f%s\n","(phi,omega,kappa,x,y,z)=(",cams[i].euler_angle[0],cams[i].euler_angle[1],cams[i].euler_angle[2],
			cams[i].camera_center[0],cams[i].camera_center[1],cams[i].camera_center[2],")");
	}
	for(int i=0;i<points.size();i++){
		printf("%s%f %f %f %f %f %f%s\n","(X,Y,Z,alpha,beta,gamma)=(",points[i].xyz[0],points[i].xyz[1],points[i].xyz[2],
			points[i].direction[0]*180/PI,points[i].direction[1]*180/PI,points[i].parallax*180/PI,")");
	}
	FILE* fp1w = NULL,* fp2w = NULL,* fp3w = NULL,* fp4w = NULL;
	fp1w = fopen("E:/zuo/projects/backup/data/narrow_test/XYZ.txt", "w");
	for (const auto& it : points) {
		fprintf(fp1w, "%f %f %f\n", it.xyz[0], it.xyz[1], it.xyz[2]);
	}
	fp2w = fopen("E:/zuo/projects/backup/data/narrow_test/Cam.txt", "w");
	for(const auto& it:cams){
		fprintf(fp2w, "%f %f %f %f %f %f\n", it.euler_angle[0], it.euler_angle[1], it.euler_angle[2],
		                                     it.camera_center[0],it.camera_center[1],it.camera_center[2]);
	}
	fp3w = fopen("E:/zuo/projects/backup/data/narrow_test/Feature.txt", "w");
	for(const auto& it:tracks){
		fprintf(fp3w,"%d",it.nview);
		for(int i=0;i<it.nview;i++){
			fprintf(fp3w,"  %d %f %f",it.obss[i].view_idx, it.obss[i].u, it.obss[i].v);
		}
		fprintf(fp3w,"\n");
	}
	fp4w = fopen("E:/zuo/projects/backup/data/narrow_test/cal.txt", "w");
	fprintf(fp4w,"%f %d %f\n",fx,0,cx);
	fprintf(fp4w,"%d %f %f\n",0,fy,cy);
	fprintf(fp4w,"%d %d %d\n",0,0,1);
	fclose(fp1w);
	fclose(fp2w);
	fclose(fp3w);
	fclose(fp4w);
	points.clear();
	cams.clear();
	tracks.clear();
}

void PBA::parallax2xyz(){
	double xj[3], xk[3];
	double Tik[3];
	int nM, nA;
	double Dik;
	double w, w2;
	for (int i = 0; i < points.size(); i++)
	{
		w = points[i].parallax;

		xj[0] = sin(points[i].direction[0]) * cos(points[i].direction[1]);
		xj[1] = sin(points[i].direction[1]);
		xj[2] = cos(points[i].direction[0]) * cos(points[i].direction[1]);

		nM = points[i].nM;
		nA = points[i].nA;

		Tik[0] = cams[nA].camera_center[0]-cams[nM].camera_center[0];
		Tik[1] = cams[nA].camera_center[1]-cams[nM].camera_center[1];
		Tik[2] = cams[nA].camera_center[2]-cams[nM].camera_center[2];
		
		Dik = sqrt(Tik[0] * Tik[0] + Tik[1] * Tik[1] + Tik[2] * Tik[2]);
		
		w2 = acos((xj[0] * Tik[0] + xj[1] * Tik[1] + xj[2] * Tik[2]) / Dik);

		xk[0] = (Dik * sin(w2 + w) * xj[0]) / sin(w);
		xk[1] = (Dik * sin(w2 + w) * xj[1]) / sin(w);
		xk[2] = (Dik * sin(w2 + w) * xj[2]) / sin(w);

		points[i].xyz[0] = cams[nM].camera_center[0] + xk[0];
		points[i].xyz[1] = cams[nM].camera_center[1] + xk[1];
		points[i].xyz[2] = cams[nM].camera_center[2] + xk[2];
	}
}
void PBA::xy_inverse_z2xyz(){
	for(int i=0;i<points.size();i++){
		points[i].xyz[0]=points[i].xy_inverse_z[0];
		points[i].xyz[1]=points[i].xy_inverse_z[1];
		points[i].xyz[2]=1/points[i].xy_inverse_z[2];
	}
}
void PBA::depth2xyz(){
	double xj[3];
	for(int i=0;i<points.size();i++){
		xj[0] = sin(points[i].world_depth[0]) * cos(points[i].world_depth[1]);
		xj[1] = sin(points[i].world_depth[1]);
		xj[2] = cos(points[i].world_depth[0]) * cos(points[i].world_depth[1]);

		points[i].xyz[0] = xj[0] * points[i].world_depth[2];
		points[i].xyz[1] = xj[1] * points[i].world_depth[2];
		points[i].xyz[2] = xj[2] * points[i].world_depth[2];
	}
}
void PBA::inverse_depth2xyz(){
	double xj[3];
	for(int i=0;i<points.size();i++){
		xj[0] = sin(points[i].world_inverse_depth[0]) * cos(points[i].world_inverse_depth[1]);
		xj[1] = sin(points[i].world_inverse_depth[1]);
		xj[2] = cos(points[i].world_inverse_depth[0]) * cos(points[i].world_inverse_depth[1]);

		points[i].xyz[0] = xj[0] / points[i].world_inverse_depth[2];
		points[i].xyz[1] = xj[1] / points[i].world_inverse_depth[2];
		points[i].xyz[2] = xj[2] / points[i].world_inverse_depth[2];
	}
}
void PBA::archored_inverse_depth2xyz(){
	double xj[3],xk[3];
	int nM;
	for(int i=0;i<points.size();i++){
		xj[0] = sin(points[i].direction[0]) * cos(points[i].direction[1]);
		xj[1] = sin(points[i].direction[1]);
		xj[2] = cos(points[i].direction[0]) * cos(points[i].direction[1]);

		nM = points[i].nM;
		xk[0] = xj[0] / points[i].inverse_depth;
		xk[1] = xj[1] / points[i].inverse_depth;
		xk[2] = xj[2] / points[i].inverse_depth;

		points[i].xyz[0] = cams[nM].camera_center[0] + xk[0];
		points[i].xyz[1] = cams[nM].camera_center[1] + xk[1];
		points[i].xyz[2] = cams[nM].camera_center[2] + xk[2];
	}
}
void PBA::archored_depth2xyz(){
	double xj[3],xk[3];
	int nM;
	for(int i=0;i<points.size();i++){
		xj[0] = sin(points[i].direction[0]) * cos(points[i].direction[1]);
		xj[1] = sin(points[i].direction[1]);
		xj[2] = cos(points[i].direction[0]) * cos(points[i].direction[1]);

		nM = points[i].nM;
		xk[0] = xj[0] * points[i].depth;
		xk[1] = xj[1] * points[i].depth;
		xk[2] = xj[2] * points[i].depth;

		points[i].xyz[0] = cams[nM].camera_center[0] + xk[0];
		points[i].xyz[1] = cams[nM].camera_center[1] + xk[1];
		points[i].xyz[2] = cams[nM].camera_center[2] + xk[2];
	}
}
void PBA::archored_xy_inverse_z2xyz(){
	int nM;
	for(int i=0;i<points.size();i++){
		nM = points[i].nM;
		points[i].xyz[0] = cams[nM].camera_center[0] + points[i].archored_xy_inverse_z[0];
		points[i].xyz[1] = cams[nM].camera_center[1] + points[i].archored_xy_inverse_z[1];
		points[i].xyz[2] = cams[nM].camera_center[2] + 1/points[i].archored_xy_inverse_z[2];
	}
}
void PBA::archored_xyz2xyz(){
	int nM;
	for(int i=0;i<points.size();i++){
		nM = points[i].nM;
		points[i].xyz[0] = cams[nM].camera_center[0] + points[i].archored_xyz[0];
		points[i].xyz[1] = cams[nM].camera_center[1] + points[i].archored_xyz[1];
		points[i].xyz[2] = cams[nM].camera_center[2] + points[i].archored_xyz[2];
	}
}
void EulerAnglesToRotationMatrix(double eulerAngles[3],double R[9]){
	double ey = eulerAngles[0];
	double ex = eulerAngles[1];
	double ez = eulerAngles[2];
	double c1 = cos(ey);   double c2 = cos(ex);   double c3 = cos(ez);
	double s1 = sin(ey);   double s2 = sin(ex);   double s3 = sin(ez);
	R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
	R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
	R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;
}
void SimEulerAngleToRotationMatrix(double eulerAngles[2],double R[9]){
	double ey = eulerAngles[0];
	double ex = eulerAngles[1];
	double c1 = cos(ey);
	double c2 = cos(ex);
	double s1 = sin(ey);
	double s2 = sin(ex);
	R[0]=c1;       R[1]=0;           R[2]=s1;
	R[3]=-s1*s2;   R[4]=c2;          R[5]=c1*s2;
	R[6]=-s1*c2;   R[7]=-s2;         R[8]=c1*c2;
}
void SimRotationMatrixToEulerAngle(double R[9],double eulerAngles[2]){
	eulerAngles[0] = atan2(R[2],R[0]);
	eulerAngles[1] = atan2(-R[7],R[4]);
}
void rotationMatrixToEulerAngles_phi_omega_kappa(double R[9], double eulerAngles[3])
{
	//assert(isRotationMatrix(R));
	double sy = sqrt(R[1] * R[1] + R[4] * R[4]);

	bool singular = sy < 1e-6;

	double phi, omega, kappa;
	if (!singular)
	{
		phi = atan2(-R[6], R[8]);
		omega = atan2(-R[7], sy);
		kappa = atan2(R[1], R[4]);
	}
	else
	{
		phi = atan2(-R[3], R[0]);
		omega = atan2(-R[7], sy);
		kappa = 0;
	}
	eulerAngles[0] = phi;
	eulerAngles[1] = omega;
	eulerAngles[2] = kappa;
}
void q2euc(const double q[4],double eu[3])
{
	double a=q[0],b=q[1],c=q[2],d=q[3];
	double R[9];
	R[0] = 1 - 2 * c * c - 2 * d * d;
	R[1] = 2 * b * c - 2 * a * d;
	R[2] = 2 * b * d + 2 * a * c;
	R[3] = 2 * b * c + 2 * a * d;
	R[4] = 1 - 2 * b * b - 2 * d * d;
	R[5] = 2 * c * d - 2 * a * b;
	R[6] = 2 * b * d - 2 * a * c;
	R[7] = 2 * c * d + 2 * a * b;
	R[8] = 1 - 2 * b * b - 2 * c * c;

	rotationMatrixToEulerAngles_phi_omega_kappa(R, eu);
}
double DotProduct(const double x[3],const double y[3]){
	return (x[0]*y[0]+x[1]*y[1]+x[2]*y[2]);
}
void axis2euc(const double angle_axis[3],double eu[3])
{
	double R[9];
	double theta2=DotProduct(angle_axis,angle_axis);
	if(theta2>numeric_limits<double>::epsilon()){
		// We want to be careful to only evaluate the square root if the
		// norm of the angle_axis vector is greater than zero. Otherwise
		// we get a division by zero.
		double theta=sqrt(theta2);
		double wx=angle_axis[0]/theta;
		double wy=angle_axis[1]/theta;
		double wz=angle_axis[2]/theta;

		double costheta=cos(theta);
		double sintheta=sin(theta);

		R[0] =     costheta   + wx*wx*(1 -    costheta);
		R[3] =  wz*sintheta   + wx*wy*(1 -    costheta);
		R[6] = -wy*sintheta   + wx*wz*(1 -    costheta);
		R[1] =  wx*wy*(1 - costheta)     - wz*sintheta;
		R[4] =     costheta   + wy*wy*(1 -    costheta);
		R[7] =  wx*sintheta   + wy*wz*(1 -    costheta);
		R[2] =  wy*sintheta   + wx*wz*(1 -    costheta);
		R[5] = -wx*sintheta   + wy*wz*(1 -    costheta);
		R[8] =     costheta   + wz*wz*(1 -    costheta);
	}else{
		// Near zero, we switch to using the first order Taylor expansion.
		R[0] =  1;
		R[3] =  angle_axis[2];
		R[6] = -angle_axis[1];
		R[1] = -angle_axis[2];
		R[4] =  1;
		R[7] =  angle_axis[0];
		R[2] =  angle_axis[1];
		R[5] = -angle_axis[0];
		R[8] = 1;
	}

	rotationMatrixToEulerAngles_phi_omega_kappa(R, eu);
}
void PBA::IFeature(){
	//feature first to image first
	for (int i = 0; i < tracks.size(); i++)
	{ 
		for (int j = 0; j < tracks[i].nview; j++)
		{ 
			int view_idx = tracks[i].obss[j].view_idx;
			double u = tracks[i].obss[j].u;
			double v = tracks[i].obss[j].v;
			double fx = 0, fy = 0, cx = 0, cy = 0;
			int cam_idx = cams[view_idx].camidx;
			fx = intrs[cam_idx-1].fx;
			fy = intrs[cam_idx-1].fy;
			cx = intrs[cam_idx-1].cx;
			cy = intrs[cam_idx-1].cy;

			double lco = atan(sqrt(((u - cx) / fx) * ((u - cx) / fx) + ((v - cy) / fy) * ((v - cy) / fy)));
			double pol = atan2(u-cx,v-cy);
			FObservation _fobs;
			_fobs.feature_idx=i;
			_fobs.u=u;
			_fobs.v=v;
			_fobs.lightcone=lco;
			_fobs.polar=pol;
			image_tracks[view_idx].push_back(_fobs);
		}
	}
	string originalPath(m_szXYZ);
	size_t pos = originalPath.find_last_of("/\\");
	string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
	string ifeature_path = parentPath + "/IFeature.txt";
	if(ifeature_path.c_str()!=NULL){
		FILE *fp = nullptr;
		fopen_s(&fp, ifeature_path.c_str(), "w");
		for (auto &kv : image_tracks) {
			int img_id = kv.first;
			auto &obs_list = kv.second;
			fprintf(fp,"%d  ",obs_list.size());
			for (auto &obs : obs_list) {
				fprintf(fp,"%d %lf %lf %lf  ",obs.feature_idx,obs.u,obs.v,obs.lightcone*180/PI);
			}
			fprintf(fp,"\n");
    	}
		fclose(fp);
	}
}
bool SolveYawFromProjection(
    double u, double v,
    double fx, double fy, double u0, double v0,
    double Cx, double Cy, double Cz,
    double ey,   
    double ex, 
    double X, double Y, double Z,
    double& yaw_out)
{
    // 1) 归一化像平面坐标
    const double xp = (u - u0) / fx;  // x'
    const double yp = (v - v0) / fy;  // y'

    // 2) 世界向量 v_w = P - C
    const double vx = X - Cx;
    const double vy = Y - Cy;
    const double vz = Z - Cz;

    // 3) 先绕 y 轴
    const double cey  = std::cos(ey);
    const double sey  = std::sin(ey);
    const double v1x = cey * vx + sey * vz;
    const double v1y = vy;
    const double v1z = -sey * vx + cey * vz;

	// 4) 再绕 x 轴
    const double cex = std::cos(ex);
    const double sex = std::sin(ex);
    const double v2x = v1x;
    const double v2y = v1y * cex + v1z * sex;
    const double v2z = -v1y * sex + v1z * cex;

    // 5) 解 yaw：cosψ 与 sinψ
    const double denom = v2x*v2x + v2y*v2y;
    const double eps = 1e-12;
    if (denom < eps) {
        // 退化：w 几乎在 z 轴上，yaw 不可观
        return false;
    }

    const double cez = (v2z * (v2x * xp + v2y * yp)) / denom;
    const double sez = (v2z * (-v2y * xp + v2x * yp)) / denom;

    // 6) atan2 得到 yaw
    yaw_out = std::atan2(sez, cez);
    return std::isfinite(yaw_out);
}
void PBA::get_yaw_from_polar(){
	for (auto &kv : image_tracks) {
		int view_idx = kv.first;
		// double R[9];
		double ey = cams[view_idx].euler_angle[0];
		double ex = cams[view_idx].euler_angle[1];
		double xc = cams[view_idx].camera_center[0];
		double yc = cams[view_idx].camera_center[1];
		double zc = cams[view_idx].camera_center[2];
		// double c1 = cos(ey);
		// double c2 = cos(ex);
		// double s1 = sin(ey);
		// double s2 = sin(ex);
		// R[0]=c1;       R[1]=0;           R[2]=s1;
		// R[3]=-s1*s2;   R[4]=c2;          R[5]=c1*s2;
		// R[6]=-s1*c2;   R[7]=-s2;         R[8]=c1*c2;

		double sum = 0;
		int j=0;
		auto &obs_list = kv.second;
		for (auto &obs : obs_list) {
			// double pa = obs.polar;
			// double lc = obs.lightcone;
			if(j==0){
				int i = obs.feature_idx;
				double 
				x = points[i].xyz[0],
				y = points[i].xyz[1],
				z = points[i].xyz[2];

				double u = obs.u;
				double v = obs.v;
				double fx = 0, fy = 0, cx = 0, cy = 0;
				int cam_idx = cams[view_idx].camidx;
				fx = intrs[cam_idx-1].fx;
				fy = intrs[cam_idx-1].fy;
				cx = intrs[cam_idx-1].cx;
				cy = intrs[cam_idx-1].cy;

				double yaw;
				// printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",u,v,fx,fy,cx,cy,xc,yc,zc,ey,ex,x,y,z);
				SolveYawFromProjection(u,v,fx,fy,cx,cy,xc,yc,zc,ey,ex,x,y,z,yaw);
				sum+=yaw;
				j++;
			}
			else{
				break;
			}
		}
		cams[view_idx].euler_angle[2]=sum/j;
	}
}
double PBA::costFuc(){
	double R[9],u,v,fx,fy,cx,cy,ey,ex,ez,c1,c2,c3,s1,s2,s3,Xc[3],p[3],xp,yp,predicted_u,predicted_v,residuals[2],err=0;
	int nobs=0,cam_idx,view_idx;
	// int k1=0,k2=0;
	// int it_id = 0;
	for (auto &kv : image_tracks) {
		// k1++;
		// k2=0;
		// it_id++;
		// if(it_id!=2){
		// 	continue;
		// }
		
		int view_idx = kv.first;
		auto &obs_list = kv.second;
		nobs+=obs_list.size();
		for (auto &obs : obs_list) {
			// k2++;
			int i=obs.feature_idx;
			u = obs.u;
			v = obs.v;
			cam_idx = cams[view_idx].camidx;
			fx = intrs[cam_idx-1].fx;
			fy = intrs[cam_idx-1].fy;
			cx = intrs[cam_idx-1].cx;
			cy = intrs[cam_idx-1].cy;
			ey=cams[view_idx].euler_angle[0];
			ex=cams[view_idx].euler_angle[1];
			ez=cams[view_idx].euler_angle[2];
			c1 = cos(ey);   c2 = cos(ex);   c3 = cos(ez);
			s1 = sin(ey);   s2 = sin(ex);   s3 = sin(ez);
			R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
			R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
			R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;
			Xc[0] = points[i].xyz[0] - cams[view_idx].camera_center[0];
			Xc[1] = points[i].xyz[1] - cams[view_idx].camera_center[1];
			Xc[2] = points[i].xyz[2] - cams[view_idx].camera_center[2];
			p[0] = R[0] * Xc[0] + R[1] * Xc[1] + R[2] * Xc[2];//X
			p[1] = R[3] * Xc[0] + R[4] * Xc[1] + R[5] * Xc[2];//Y
			p[2] = R[6] * Xc[0] + R[7] * Xc[1] + R[8] * Xc[2];//Z
			// Normalize
			xp = p[0] / p[2];
			yp = p[1] / p[2];

			// Project
			predicted_u = fx * xp + cx;
			predicted_v = fy * yp + cy;

			// if(k1==1 && k2==1){
			// 	printf("%d %f %f %f %f %f %f\n",view_idx, cams[view_idx].euler_angle[0],cams[view_idx].euler_angle[1],cams[view_idx].euler_angle[2],
			// 	cams[view_idx].camera_center[0],cams[view_idx].camera_center[1],cams[view_idx].camera_center[2]);
			// 	printf("%d %f %f %f\n",i, points[i].xyz[0],points[i].xyz[1],points[i].xyz[2]);
			// 	printf("%f %f\n",u,v);
			// 	printf("%f %f\n\n",predicted_u,predicted_v);
			// }
			
			// Residual
			residuals[0] = predicted_u - u;
			residuals[1] = predicted_v - v;
			err+=(residuals[0]*residuals[0]+residuals[1]*residuals[1]);
		}
	}
	return sqrt(err/(2*nobs));
}
void PBA::randCam(ParameterType paramtype){
	std::random_device rd;
	std::mt19937 gen(rd()); 
	std::uniform_real_distribution<double> dist1(-PI/4, PI/4),dist2(-PI/4,PI/4),dist3(-2*PI,2*PI),
	dist4(-6,6),dist5(-2,2),dist6(-2,2);

	string originalPath(m_szXYZ);
	size_t pos = originalPath.find_last_of("/\\");
	string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
	string ifeature_path = parentPath + "/Cam.txt";
	if(ifeature_path.c_str()!=NULL){
		FILE *fp = nullptr;
		fopen_s(&fp, ifeature_path.c_str(), "w");
		for(int i=0;i<cams.size();i++){
			if(paramtype==rotation_translation){
				cams[i].camera_center[0]=dist4(gen);
				cams[i].camera_center[1]=dist5(gen);
				cams[i].camera_center[2]=dist6(gen);
			}
			// cams[i].euler_angle[0]=dist1(gen);
			// cams[i].euler_angle[1]=dist2(gen);
			// cams[i].euler_angle[2]=dist3(gen);

			cams[i].euler_angle[0]=0;
			cams[i].euler_angle[1]=0;
			cams[i].euler_angle[2]=0;
			
			fprintf(fp,"%lf %lf %lf %lf %lf %lf %d\n",cams[i].euler_angle[0],cams[i].euler_angle[1],cams[i].euler_angle[2],
			cams[i].camera_center[0],cams[i].camera_center[1],cams[i].camera_center[2],cams[i].camidx);
		}
		fclose(fp);
	}
}
// 将角度规范到 [-PI/2, PI/2]
double NormalizeHalfPi(double angle) {
    // 对 pitch/roll 来说，超过 ±π/2 需要取补角
    while (angle > PI/2)  angle -= PI;
    while (angle < -PI/2) angle += PI;
    return angle;
}
bool PBA::ba_run(char* szCam,
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
    manifoldtype manitype)
{
	m_szCameraInit = szCam;                                                       
	m_szFeatures   = szFea;                                                       
	m_szCalibration= szCalib;                                                     
	m_szXYZ        = szXYZ;                                                                                                                                                                                                            
	m_szCamePose   = szPose;                                                      
	m_sz3Dpts      = sz3D;                                                        
	m_szReport     = szReport;                                                    

	const char* object_point_type;
	const char* image_point_type;
	const char* rotation_3d_type;
	const char* parameter_type;
	const char* manifold_type;
	switch (optype)
	{
	case xyz:
		object_point_type="xyz";break;
	case xy_inverse_z:	
		object_point_type="xy_inverse_z";break;
	case depth:
		object_point_type="depth";break;
	case inverse_depth:
		object_point_type="inverse_depth";break;
	case archored_xyz:
		object_point_type="archored_xyz";break;
	case archored_xy_inverse_z:
		object_point_type="archored_xy_inverse_z";break;
	case archored_depth:
		object_point_type="archored_depth";break;
	case archored_inverse_depth:
		object_point_type="archored_inverse_depth";break;
	case parallax:
		object_point_type="parallax";break;	
	}
	switch (r3dtype)
	{
	case euler_angle:
		rotation_3d_type="euler_angle";break;
	case axis_angle:
		rotation_3d_type="angle_axis";break;
	case quaternion:
		rotation_3d_type="quaternion";break;
	}
	switch(iptype)
	{
	case uv:
		image_point_type="uv";break;
	case light_cone:
		image_point_type="light_cone";break;
	}
	switch(paramtype)
	{
	case rotation_translation_landmark:
		parameter_type="rotation_translation_landmark";break;
	case rotation_landmark:
		parameter_type="rotation_landmark";break;
	case translation_landmark:
		parameter_type="translation_landmark";break;
	case rotation_translation:
		parameter_type="rotation_translation";break;
	case rotation:
		parameter_type="rotation";break;
	case translation:
		parameter_type="translation";break;
	case landmark:
		parameter_type="landmark";break;
	}
	switch(manitype)
	{
	case lie:
		manifold_type="lie";break;
	case quaternion_manifold:
		manifold_type="quaternion_manifold";break;
	case sphere_manifold:
		manifold_type="sphere_manifold";break;
	case line_manifold:
		manifold_type="line_manifold";break;
	case euclidean_manifold:
		manifold_type="euclidean_manifold";break;
	case none:
		manifold_type="none";break;
	}
	printf("%s\n",m_szXYZ);
	printf("parameterization\n");
	printf("%s %s %s %s %s\n","object point  |  ","3d rotation  |  ","image point  |  ","parameter  |  ","manifold");
	printf("%s %s %s %s %s\n",object_point_type,rotation_3d_type,image_point_type,parameter_type,manifold_type);

	string str = string(object_point_type)+"_"+rotation_3d_type+"_"+image_point_type;
	// narrow_test();
	if(iptype==uv){
		ba_initialize(m_szCameraInit, m_szFeatures, m_szCalibration, m_szXYZ); 
		if(paramtype!=rotation_translation_landmark){
			randCam(paramtype);//change euler angle or camera center and save as Cam.txt
		}
	}else if(iptype==light_cone){
		ba_initialize(m_szCameraInit, m_szFeatures, m_szCalibration, m_szXYZ); 
	}
	
	ceres::Problem problem;
	ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);
	double initial_cost,final_cost,real_final_cost;
	int nobs=0;
	if(paramtype==rotation || paramtype==rotation_translation || paramtype==translation)
	{
		IFeature();
		// int im_id = 0;
		for (auto &kv : image_tracks) {
			// im_id++;
			// if(im_id!=2){
			// 	continue;
			// }
			int view_idx = kv.first;
			auto &obs_list = kv.second;
			nobs+=obs_list.size();
			for (auto &obs : obs_list) {
				int i=obs.feature_idx;
				double u = obs.u;
				double v = obs.v;
				double lco = obs.lightcone;
				double fx = 0, fy = 0, cx = 0, cy = 0;
				int cam_idx = cams[view_idx].camidx;
				fx = intrs[cam_idx-1].fx;
				fy = intrs[cam_idx-1].fy;
				cx = intrs[cam_idx-1].cx;
				cy = intrs[cam_idx-1].cy;
				
				/*Parameterization of image point*/
				if(iptype==light_cone && r3dtype==euler_angle && paramtype==rotation_translation){
					ceres::CostFunction *cost_function = lc_euler_angle_rc::Create(lco,points[i].xyz[0],points[i].xyz[1],points[i].xyz[2]);
					problem.AddResidualBlock(
						cost_function,
						nullptr, //loss_function
						&cams[view_idx].euler_angle[0],
						&cams[view_idx].camera_center[0]);
				}
				else if(iptype==uv && r3dtype==euler_angle && paramtype==rotation_translation){
					ceres::CostFunction *cost_function = uv_euler_angle_rc::Create(u, v, fx, fy, cx, cy, 
																points[i].xyz[0],points[i].xyz[1],points[i].xyz[2]);
					problem.AddResidualBlock(
						cost_function,
						nullptr,
						&cams[view_idx].euler_angle[0],
						&cams[view_idx].camera_center[0]);
				}
				else if(iptype==light_cone && r3dtype==euler_angle && paramtype==rotation){
					ceres::CostFunction *cost_function = lc_euler_angle_r::Create(lco,
						points[i].xyz[0],points[i].xyz[1],points[i].xyz[2],
						cams[view_idx].camera_center[0],cams[view_idx].camera_center[1],cams[view_idx].camera_center[2]);
					problem.AddResidualBlock(
						cost_function,
						// loss_function,
						nullptr, // loss_function
						&cams[view_idx].euler_angle[0]);
				}
				else if(iptype==uv && r3dtype==euler_angle && paramtype==rotation){
					ceres::CostFunction *cost_function = uv_euler_angle_r::Create(u, v, fx, fy, cx, cy, 
												points[i].xyz[0],points[i].xyz[1],points[i].xyz[2],
												cams[view_idx].camera_center[0],cams[view_idx].camera_center[1],cams[view_idx].camera_center[2]);
					problem.AddResidualBlock(
						cost_function,
						// loss_function,
						nullptr,
						&cams[view_idx].euler_angle[0]);
				}
			}
		}
		// if(iptype==light_cone){
		// 	initial_cost=costFuc();
		// }
		initial_cost=costFuc();
		// if(iptype == uv){
		// 	initial_cost = lc_rep();
		// }
	}
	else if(paramtype==rotation_translation_landmark){
		if(manitype == lie){
			auto se3manifold = new Sophus::Manifold<Sophus::SE3>();
			// auto quat_manifold = new ceres::QuaternionManifold();
			for (int i = 0; i < points.size();i++){
				problem.AddParameterBlock(points[i].xyz, 3);
			}
			for (int i = 0; i < cams.size(); i++) {
				problem.AddParameterBlock(cams[i].se3, 7, se3manifold);
				// problem.AddParameterBlock(cams[i].quat, 4, quat_manifold);
				// problem.AddParameterBlock(cams[i].translation, 3);
			}
		}
		else if(manitype == quaternion_manifold){
			auto quat_manifold = new ceres::QuaternionManifold();
			for (int i = 0; i < points.size();i++){
				problem.AddParameterBlock(points[i].xyz, 3);
			}
			for (int i = 0; i < cams.size(); i++) {
				problem.AddParameterBlock(cams[i].quat, 4, quat_manifold);
				problem.AddParameterBlock(cams[i].camera_center, 3);
			}
		}

		for (int i = 0; i < tracks.size(); i++)
		{ 
			nobs += tracks[i].nview;
			int nM = points[i].nM;
			int nA = points[i].nA;

			for (int j = 0; j < tracks[i].nview; j++)
			{ 
				int view_idx = tracks[i].obss[j].view_idx;
				double u = tracks[i].obss[j].u;
				double v = tracks[i].obss[j].v;
				double fx = 0, fy = 0, cx = 0, cy = 0;
				int cam_idx = cams[view_idx].camidx;
				fx = intrs[cam_idx-1].fx;
				fy = intrs[cam_idx-1].fy;
				cx = intrs[cam_idx-1].cx;
				cy = intrs[cam_idx-1].cy;
		
				int nP = view_idx;
				double lco = atan(sqrt(((u - cx) / fx) * ((u - cx) / fx) + ((v - cy) / fy) * ((v - cy) / fy)));
				double pol = atan2(u-cx,v-cy);
				tracks[i].obss[j].lightcone=lco;
				tracks[i].obss[j].polar=pol;
				
				/*Parameterization of object point*/
				//Zero archor
				if(optype==xyz && r3dtype==euler_angle && iptype==uv){
					ceres::CostFunction *cost_function = xyz_euler_angle_uv::Create(u, v, fx, fy, cx, cy);
					problem.AddResidualBlock(
						cost_function,
						nullptr, 
						&cams[view_idx].euler_angle[0],
						&cams[view_idx].camera_center[0],
						&points[i].xyz[0]);
				}
				else if(optype==xy_inverse_z && r3dtype==euler_angle && iptype==uv){
						ceres::CostFunction *cost_function = xy_inverse_z_euler_angle_uv::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].euler_angle[0],
							&cams[nP].camera_center[0],
							&points[i].xy_inverse_z[0]);
				}
				else if(optype==depth && r3dtype==euler_angle && iptype==uv){
					ceres::CostFunction *cost_function = depth_euler_angle_uv::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].euler_angle[0],
							&cams[nP].camera_center[0],
							&points[i].world_depth[0]);
				}
				else if(optype==inverse_depth && r3dtype==euler_angle && iptype==uv){
					ceres::CostFunction *cost_function = inverse_depth_euler_angle_uv::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].euler_angle[0],
							&cams[nP].camera_center[0],
							&points[i].world_inverse_depth[0]);
				}
				
				//One archor
				if(optype==archored_xyz && r3dtype==euler_angle && iptype==uv){
					if(nP==nM){
						ceres::CostFunction *cost_function = archored_xyz_euler_angle_uv_nM::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function,
							&cams[nP].euler_angle[0],
							&points[i].archored_xyz[0]);
					}else{
						ceres::CostFunction *cost_function = archored_xyz_euler_angle_uv_nP::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].euler_angle[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&points[i].archored_xyz[0]);
					}
				}
				else if(optype==archored_xy_inverse_z && r3dtype==euler_angle && iptype==uv){
					if(nP==nM){
						ceres::CostFunction *cost_function = archored_xy_inverse_z_euler_angle_uv_nM::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function,
							&cams[nP].euler_angle[0],
							&points[i].archored_xy_inverse_z[0]);
					}else{
						ceres::CostFunction *cost_function = archored_xy_inverse_z_euler_angle_uv_nP::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].euler_angle[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&points[i].archored_xy_inverse_z[0]);
					}
				}
				else if(optype==archored_depth && r3dtype==euler_angle && iptype==uv){
					if(nP==nM){
						ceres::CostFunction *cost_function = archored_depth_euler_angle_uv_nM::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function,
							&cams[nP].euler_angle[0],
							&points[i].direction[0]);
					}else{
						ceres::CostFunction *cost_function = archored_depth_euler_angle_uv_nP::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].euler_angle[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&points[i].direction[0],
							&points[i].depth);
					}
				}
				else if(optype==archored_inverse_depth && r3dtype==euler_angle && iptype==uv){
					// problem.AddParameterBlock(&points[i].idp[0], 3); 
					if(nP==nM){
						// std::vector<int> constant_parameters = {2};
						// auto* subset_manifold = new ceres::SubsetManifold(3, constant_parameters);
						// problem.SetManifold(&points[i].idp[0], subset_manifold);
						ceres::CostFunction *cost_function = archored_inverse_depth_euler_angle_uv_nM::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function,
							&cams[nP].euler_angle[0],
							&points[i].direction[0]);
					}else{
						ceres::CostFunction *cost_function = archored_inverse_depth_euler_angle_uv_nP::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].euler_angle[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&points[i].direction[0],
							&points[i].inverse_depth);
					}
				}

				//Two archors
				if(optype==parallax && r3dtype==euler_angle && iptype==uv){
					if(nP==nM){
						ceres::CostFunction *cost_function = parallax_euler_angle_uv_nM::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function,
							&cams[nP].euler_angle[0],
							&points[i].direction[0]);
					}else if(nP==nA){
						ceres::CostFunction *cost_function = parallax_euler_angle_uv_nA::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].euler_angle[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&points[i].direction[0],
							&points[i].parallax);
					}else{
						ceres::CostFunction *cost_function = parallax_euler_angle_uv_nP::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].euler_angle[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&cams[nA].camera_center[0],
							&points[i].direction[0],
							&points[i].parallax);
					}
				}

				/*Parameterization of 3d rotation*/
				//axis angle
				if(optype==xyz && r3dtype==axis_angle && iptype==uv){
						ceres::CostFunction *cost_function = axis_angle_xyz_uv::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].axis_angle[0],
							&cams[nP].camera_center[0],
							&points[i].xyz[0]);
				}
				else if(optype==xyz && r3dtype==quaternion && iptype==uv){
						ceres::CostFunction *cost_function = quaternion_xyz_uv::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].quat[0],
							&cams[nP].camera_center[0],
							&points[i].xyz[0]);
				}
				//quaternion
				if(optype==parallax && r3dtype==axis_angle && iptype==uv){
					if(nP==nM){
						ceres::CostFunction *cost_function = axis_angle_parallax_uv_nM::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function,
							&cams[nP].axis_angle[0],
							&points[i].direction[0]);
					}else if(nP==nA){
						ceres::CostFunction *cost_function = axis_angle_parallax_uv_nA::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].axis_angle[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&points[i].direction[0],
							&points[i].parallax);
					}else{
						ceres::CostFunction *cost_function = axis_angle_parallax_uv_nP::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].axis_angle[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&cams[nA].camera_center[0],
							&points[i].direction[0],
							&points[i].parallax);
					}
				}
				else if(optype==parallax && r3dtype==quaternion && iptype==uv){
					if(nP==nM){
						ceres::CostFunction *cost_function = quaternion_parallax_uv_nM::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function,
							&cams[nP].quat[0],
							&points[i].direction[0]);
					}else if(nP==nA){
						ceres::CostFunction *cost_function = quaternion_parallax_uv_nA::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].quat[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&points[i].direction[0],
							&points[i].parallax);
					}else{
						ceres::CostFunction *cost_function = quaternion_parallax_uv_nP::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].quat[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&cams[nA].camera_center[0],
							&points[i].direction[0],
							&points[i].parallax);
					}
				}
				
				/*Manifold*/
				if(optype==xyz && iptype==uv && manitype==lie){
					ceres::CostFunction *cost_function = xyz_uv_lie_manifold::Create(u, v, fx, fy, cx, cy);
					problem.AddResidualBlock(
						cost_function,
						nullptr, // loss_function,
						&cams[nP].se3[0],
						&points[i].xyz[0]);
					// problem.AddResidualBlock(
					// 	cost_function,
					// 	nullptr, // loss_function,
					// 	&cams[nP].quat[0],
					// 	&cams[nP].translation[0],
					// 	&points[i].xyz[0]);
				}
				else if(optype==parallax&& iptype==uv && manitype==quaternion_manifold){
					if(nP==nM){
						ceres::CostFunction *cost_function = parallax_uv_quaternion_manifold_nM::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function,
							&cams[nP].quat[0],
							&points[i].direction[0]);
					}else if(nP==nA){
						ceres::CostFunction *cost_function = parallax_uv_quaternion_manifold_nA::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].quat[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&points[i].direction[0],
							&points[i].parallax);
					}else{
						ceres::CostFunction *cost_function = parallax_uv_quaternion_manifold_nP::Create(u, v, fx, fy, cx, cy);
						problem.AddResidualBlock(
							cost_function,
							nullptr, //loss_function, 
							&cams[nP].quat[0],
							&cams[nP].camera_center[0],
							&cams[nM].camera_center[0],
							&cams[nA].camera_center[0],
							&points[i].direction[0],
							&points[i].parallax);
					}
				}
			}
		}
	}

	// double ini_err = lc_rep();
	
	std::vector<double> cost_per_iteration,cost_per_iteration1;
	int num_cams = m_ncams;
	std::vector<std::vector<double>> euler_history,euler_history1;
	// Camera* cameras = new Camera[num_cams];
	CostRecorderCallback* callback=NULL;
	CostRecorderCallback_simple * callback_simple=NULL;
	if(iptype==light_cone){
		callback = new CostRecorderCallback(cost_per_iteration,euler_history,cams,num_cams);
	}else if(iptype==uv){
		callback_simple=new CostRecorderCallback_simple(cost_per_iteration);
	}
	printf("Ceres Options\n");
	ceres::Solver::Options options;
	options.logging_type = ceres::SILENT;
	options.minimizer_progress_to_stdout = true;
	options.num_threads = 8;
	options.max_num_iterations = 100;  
	// options.trust_region_strategy_type = ceres::DOGLEG;         
	options.linear_solver_type = ceres::SPARSE_SCHUR;            
	options.sparse_linear_algebra_library_type = ceres::CUDA_SPARSE;

	if(iptype==light_cone){
		options.callbacks.push_back(callback);
	}else if(iptype==uv){
		options.callbacks.push_back(callback_simple);
	}
	
	options.update_state_every_iteration = true; 
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.FullReport() << "\n";

	// if(iptype==light_cone){
	// 	printf("get yaw from polar\n");
	// 	get_yaw_from_polar();
	// 	final_cost= costFuc();
	// }
	string originalPath(m_szXYZ);
	size_t pos = originalPath.find_last_of("/\\");
	string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
	final_cost= costFuc();
	if(iptype == uv){
		// final_cost = lc_rep();
	}
	else if(iptype == light_cone){
		//输出结果
		string opt_cam_tmp = parentPath + "/Cam_" + str + "_temp.txt";
		if(opt_cam_tmp.c_str()!=NULL){
			FILE *fp = nullptr;
			fopen_s(&fp, opt_cam_tmp.c_str(), "w");
			if(r3dtype==axis_angle){
				for(int i=0;i<cams.size();i++){
					axis2euc(cams[i].axis_angle,cams[i].euler_angle);
				}
			}
			if(r3dtype==quaternion){
				for(int i=0;i<cams.size();i++){
					q2euc(cams[i].quat,cams[i].euler_angle);
				}
			}
			for (int i = 0; i < cams.size(); i++){
				// double R[9];
				// SimEulerAngleToRotationMatrix(cams[i].euler_angle,R);
				// SimRotationMatrixToEulerAngle(R,cams[i].euler_angle);
				// double ey  = NormalizeHalfPi(cams[i].euler_angle[0]);
    			// double ex = NormalizeHalfPi(cams[i].euler_angle[1]);	
				// cams[i].euler_angle[0] = ey;
				// cams[i].euler_angle[1] = ex;
				// if(i!=1){
				// 	continue;
				// }
				fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", cams[i].euler_angle[0], cams[i].euler_angle[1], cams[i].euler_angle[2],
											cams[i].camera_center[0],cams[i].camera_center[1],cams[i].camera_center[2]);
			}
			fclose(fp);
		}
		printf("get yaw from polar\n");
		get_yaw_from_polar();
		ceres::Problem problem1;
		// int im_id = 0;
		for (auto &kv : image_tracks) {
			int view_idx = kv.first;
			// im_id++;
			// if(im_id!=2){
			// 	continue;
			// }
			auto &obs_list = kv.second;
			for (auto &obs : obs_list) {
				int i=obs.feature_idx;
				double u = obs.u;
				double v = obs.v;
				double fx = 0, fy = 0, cx = 0, cy = 0;
				int cam_idx = cams[view_idx].camidx;
				fx = intrs[cam_idx-1].fx;
				fy = intrs[cam_idx-1].fy;
				cx = intrs[cam_idx-1].cx;
				cy = intrs[cam_idx-1].cy;
				if(paramtype==rotation){
					ceres::CostFunction *cost_function = uv_euler_angle_r::Create(u, v, fx, fy, cx, cy, 
									points[i].xyz[0],points[i].xyz[1],points[i].xyz[2],
									cams[view_idx].camera_center[0],cams[view_idx].camera_center[1],cams[view_idx].camera_center[2]);
					problem1.AddResidualBlock(
						cost_function,
						// loss_function,
						nullptr,
						&cams[view_idx].euler_angle[0]);
				}
				else if(paramtype==rotation_translation){
					ceres::CostFunction *cost_function = uv_euler_angle_rc::Create(u, v, fx, fy, cx, cy, 
											points[i].xyz[0],points[i].xyz[1],points[i].xyz[2]);
					problem1.AddResidualBlock(
						cost_function,
						nullptr,
						&cams[view_idx].euler_angle[0],
						&cams[view_idx].camera_center[0]);
				}
			}
		}
		
		CostRecorderCallback* callback = new CostRecorderCallback(cost_per_iteration1,euler_history1,cams,num_cams);
		printf("total refine\n");
		ceres::Solver::Options options;
		options.logging_type = ceres::SILENT;
		options.minimizer_progress_to_stdout = true;
		options.num_threads = 8;
		options.max_num_iterations = 100;  
		// options.trust_region_strategy_type = ceres::DOGLEG;         
		options.linear_solver_type = ceres::SPARSE_SCHUR;            
		options.sparse_linear_algebra_library_type = ceres::CUDA_SPARSE;

		options.callbacks.push_back(callback);
		options.update_state_every_iteration = true; 
		ceres::Solver::Summary summary;
		ceres::Solve(options, &problem1, &summary);
		std::cout << summary.FullReport() << "\n";
		real_final_cost= costFuc();
	}
	
	// // printf("%f %s %f\n", ini_err, " to ", fin_err);


	string converge_curve = parentPath + "/convergence_" + str + ".txt";
	string opt_point = parentPath + "/XYZ_" + str + ".ply";
	string opt_cam = parentPath + "/Cam_" + str + ".txt";
	if(opt_cam.c_str()!=NULL){
		FILE *fp = nullptr;
		fopen_s(&fp, opt_cam.c_str(), "w");
		if(r3dtype==axis_angle){
			for(int i=0;i<cams.size();i++){
				axis2euc(cams[i].axis_angle,cams[i].euler_angle);
			}
		}
		if(r3dtype==quaternion){
			for(int i=0;i<cams.size();i++){
				q2euc(cams[i].quat,cams[i].euler_angle);
			}
		}
		for (int i = 0; i < cams.size(); i++){
			// if(i!=1){
			// 	continue;
			// }
			double R[9];
			EulerAnglesToRotationMatrix(cams[i].euler_angle,R);
			rotationMatrixToEulerAngles_phi_omega_kappa(R,cams[i].euler_angle);
			// double ey  = NormalizeHalfPi(cams[i].euler_angle[0]);
			// double ex = NormalizeHalfPi(cams[i].euler_angle[1]);	
			// cams[i].euler_angle[0] = ey;
			// cams[i].euler_angle[1] = ex;
			fprintf(fp, "%lf %lf %lf %lf %lf %lf\n", cams[i].euler_angle[0], cams[i].euler_angle[1], cams[i].euler_angle[2],
			                             cams[i].camera_center[0],cams[i].camera_center[1],cams[i].camera_center[2]);
		}
		fclose(fp);
	}
	//save features xyz
	if (opt_point.c_str() != NULL){
		FILE *fp = nullptr;
		fopen_s(&fp, opt_point.c_str(), "w");
		fprintf(fp, "%s\n", "ply");
		fprintf(fp, "%s\n", "format ascii 1.0");
		fprintf(fp, "%s %d\n", "element vertex", m_n3Dpts);
		fprintf(fp, "%s\n", "property float x");
		fprintf(fp, "%s\n", "property float y");
		fprintf(fp, "%s\n", "property float z");
		fprintf(fp, "%s\n", "end_header");
	
		//zero archor
		if(optype==xy_inverse_z){
			xy_inverse_z2xyz();
		}else if(optype==depth){
			depth2xyz();
		}else if(optype==inverse_depth){
			inverse_depth2xyz();
		}

		//one archor
		if(optype==archored_xyz){
			archored_xyz2xyz();
		}else if(optype==archored_xy_inverse_z){
			archored_xy_inverse_z2xyz();
		}else if(optype==archored_depth){
			archored_depth2xyz();
		}else if(optype==archored_inverse_depth){
			archored_inverse_depth2xyz();
		}

		//two archors
		if(optype==parallax){
			parallax2xyz();
		}

		for (int i = 0; i < points.size(); i++)
			fprintf(fp, "%lf %lf %lf\n", points[i].xyz[0], points[i].xyz[1], points[i].xyz[2]);
		fclose(fp);
	}

	// printf("%d\n",num_cams);
	if(converge_curve.c_str()!=NULL){
		FILE *fp=nullptr;
		fopen_s(&fp,converge_curve.c_str(),"w");
		for (int i = 0; i < cost_per_iteration.size(); ++i) {
			if(iptype==light_cone){
				for(int j=0;j<num_cams;j++){
					// if(j!=1){
					// 	continue;
					// }
					cams[j].euler_angle[0] = euler_history[i][j*6 + 0];
					cams[j].euler_angle[1] = euler_history[i][j*6 + 1];
					cams[j].euler_angle[2] = euler_history[i][j*6 + 2];
					if(paramtype==rotation_translation){
						cams[j].camera_center[0] = euler_history[i][j*6 + 3];
						cams[j].camera_center[1] = euler_history[i][j*6 + 4];
						cams[j].camera_center[2] = euler_history[i][j*6 + 5];
					}
				}
				// printf("%d %f\n",i,euler_history[i][0]);
				double cost = costFuc();
				fprintf(fp, "%d %lf\n", i, cost);
			}else{
				fprintf(fp, "%d %lf\n", i, sqrt(cost_per_iteration[i]/nobs));
			}

			// fprintf(fp, "%d %lf %lf\n", i, cost,sqrt(cost_per_iteration[i]/nobs));
			
		}
		if(iptype==light_cone){
			for (int i = 0; i < cost_per_iteration1.size(); ++i) {
				for(int j=0;j<num_cams;j++){
					// if(j!=1){
					// 	continue;
					// }
					cams[j].euler_angle[0] = euler_history1[i][j*6 + 0];
					cams[j].euler_angle[1] = euler_history1[i][j*6 + 1];
					cams[j].euler_angle[2] = euler_history1[i][j*6 + 2];
					if(paramtype==rotation_translation){
						cams[j].camera_center[0] = euler_history1[i][j*6 + 3];
						cams[j].camera_center[1] = euler_history1[i][j*6 + 4];
						cams[j].camera_center[2] = euler_history1[i][j*6 + 5];
					}
				}
				// printf("%d %f\n",i,euler_history[i][0]);
				double cost = costFuc();
				// fprintf(fp, "%d %lf\n", i, sqrt(cost_per_iteration[i]/nobs));
				// fprintf(fp, "%d %lf %lf\n", i+cost_per_iteration.size(), cost,sqrt(cost_per_iteration1[i]/nobs));
				fprintf(fp, "%d %lf\n", i+cost_per_iteration.size(), cost);
			}
		}
	}
	if(iptype == uv){
		// printf("%f %s %f\n", initial_cost, " to ", final_cost);
		initial_cost = summary.initial_cost;
		final_cost = summary.final_cost;
		printf("%f %s %f\n", sqrt(initial_cost/nobs), " to ", sqrt(final_cost/nobs));
	}
	if(iptype == light_cone)
	{
		// initial_cost = summary.initial_cost;
		// final_cost = summary.final_cost;
		// printf("%f %s %f\n", sqrt(initial_cost/nobs), " to ", sqrt(final_cost/nobs));
		printf("%f %s %f %s %f\n", initial_cost, " to ", final_cost," to ",real_final_cost);
	}
	// if(iptype==light_cone){
		
	// }else{
		// initial_cost = summary.initial_cost;
		// final_cost = summary.final_cost;
		// printf("%f %s %f\n", sqrt(initial_cost/nobs), " to ", sqrt(final_cost/nobs));
	// }
	// // printf("%f %s %f\n", initial_cost*2/nobs, " to ", final_cost*2/nobs);
	// printf("%f %s %f\n", sqrt(initial_cost/nobs), " to ", sqrt(final_cost/nobs));

	return true;
}

bool PBA::ba_initialize( char* szCamera, char* szFeature,  char* szCalib, char* szXYZ )
{
	// printf("BA: Bundle Adjustment Version 1.0\n");
	FILE* fp = nullptr;

	//must input initial initial camera pose file and projection image points file
	fopen_s(&fp, szCamera, "r" );
	if ( fp == NULL )
	{
		fprintf( stderr, "BA: Missing initial camera poses file! \n");
		exit(1);
	}
	else
		fclose(fp);

	fopen_s(&fp, szFeature, "r" );
	if ( fp == NULL )
	{	
		fprintf( stderr, "BA: Missing feature projection points file! \n");
		exit(1); 
	}
	else
		fclose(fp);

	if (szCalib != NULL)
	{
		FILE* fpc = nullptr;
		fopen_s(&fpc, szCalib, "r");
		nc_ = findNcameras(fpc) / 3;
		fclose(fpc);
		fpc = nullptr;
		// printf("%d\n", nc_);
		m_bFocal = false;
		m_K = (double*)malloc(9 * nc_ * sizeof(double));
		ba_readCameraPoseration(szCalib, m_K);
	}

	if ( szXYZ != NULL )
		m_bProvideXYZ = true;	
	//read camera pose & features images projs, and initialize features points( three kinds of angle )
	pba_readAndInitialize( szCamera, szFeature,szCalib, &m_ncams, &m_n3Dpts, &m_n2Dprojs,&m_motstruct,//number of camera, 3D points, 2D projection points,6 camera pose and 3 feature parameters
		&m_imgpts, &m_archor, &m_vmask, &m_umask, &m_photo, &m_feature, &m_archorSort);

	return true;
}


void PBA::pba_readProjectionAndInitilizeFeature(FILE *fp,
	double *params, double *projs, char *vmask, int ncams, 
	int *archor,char* umask,int* nphoto, int* nfeature, int* archorSort )
{
	int n;
	int nframes;
	int ptno = 0, cur;

	int nproj2D = 0;
	
	int count = 0;

	int frameno;
	int feastart = 0;

	int nP, nP2;
	
	double* ptr1 = projs;

	int i, j; 
	int  sum, cnp = 6;

	int nFlag;
	
	int *ptr2;
	bool bAdjust;	

	m_smask = (char*)malloc(m_ncams*m_ncams*sizeof(char));
	memset( m_smask, 0, m_ncams*m_ncams*sizeof(char) );

	int* archorEx = new int[m_n3Dpts*2];

	bool bM, bN;
	int max_nframes = -1;
	//read all projection point, initialize three feature angle at the same time
	while(!feof(fp))
	{
		nFlag = 0;
		n = readNInts( fp, &nframes, 1 );  
		if( n!=1 )
			break;

		Track track;
		track.nview = nframes;

		archor[ptno*3] = nframes;
		cur = 0;
		bM = bN = false;
		for( i=0, sum = 0; i<nframes; ++i )
		{
			n = readNInts( fp, &frameno, 1 );
			nphoto[nproj2D] = frameno;
			nfeature[nproj2D] = ptno;
			nproj2D++;

			if(frameno>=ncams)
			{
				fprintf(stderr, "ParallaxBA: the image No. of projection point is out of max image No.\n");
				return;
			}

			n += readNDoubles( fp, ptr1, 2 ); 

			Observation obs;
			obs.view_idx = frameno;
			obs.u = ptr1[0];
			obs.v = ptr1[1];
			track.obss.push_back(obs);

			ptr1+=2;
			if(n!=3)
			{
				fprintf(stderr, "ParallaxBA:reading image projections wrong!\n");
				return;
			}

			if ( bM && bN )
			{
				ptr2 = archorSort+ptno*2;
				//bAdjust = pba_initializeOtheArchors_Mindw( //_Mindw
				bAdjust = pba_initializeOtheArchors( //Maxdw
					projs+feastart*2,
					nphoto+feastart,
					m_motstruct,
					m_K,
					m_motstruct + m_ncams*cnp + ptno * 3, 
					ptr2,
					sum,
					i,
					ptno );
				if ( bAdjust )
				{
					archor[ptno*3+1] = *(nphoto+feastart+ptr2[0]);
					archor[ptno*3+2] = *(nphoto+feastart+ptr2[1]);

					archorEx[ptno*2] = ptr2[0];
					archorEx[ptno*2+1] = ptr2[1];
				}
				sum++;
			}

			if ( bM && !bN )
			{	
				bool bLast = (i == nframes-1);
				bool bT = pba_initializeAssoArchor( 
					projs+feastart*2,	
					nphoto+feastart,	
					m_motstruct,
					m_K,
					m_motstruct+m_ncams*cnp+ptno*3,
					0,
					1,
					ptno,
					bLast );

				if (bT)
				{
					archorSort[ptno*2+1] = i;
					archor[ptno*3+2] = nphoto[count];
					sum++;

					archorEx[ptno*2+1] = i;

					bN = true;
				}
			}

			if ( !bM )
			{
				bool bLast = (i == nframes-2);
				bool bT = pba_initializeMainArchor( 
					projs+feastart*2,	
					m_motstruct,		
					m_K,				
					m_motstruct+m_ncams*cnp+ptno*3,
					nphoto[count],		
					ptno,				
					m_KR );				

				archorSort[ptno*2] = i;
				archor[ptno*3+1] = nphoto[count];
				sum++;

				archorEx[ptno*2] = i;
				bM = true;
			}	
			count++;	
		}

		tracks.push_back(track);
		//set masks for U and S matrix             
		for( i = 0; i < nframes; i++ )
		{
			nP = nphoto[feastart+i];                         
			int nM_ = archor[ptno * 3 + 1];
			int nA_ = archor[ptno * 3 + 2];
			int tmp3 = archor[ptno * 3 + 3];
			
			if (nM_<nP)                         
				umask[nM_*(ncams)+nP] = 1;
			else                                             
				umask[nP*(ncams)+nM_] = 1;

			
			if (nA_<nP)                         
				umask[nA_*(ncams)+nP] = 1;
			else                                             
				umask[nP*(ncams)+nA_] = 1;

			umask[nP*ncams+nP] = 1;

			for ( j = i; j < nframes; j++  )
			{
				nP2 = nphoto[feastart+j];                

				if ( nP == nP2 )                              
					m_smask[nP*m_ncams+nP2] = 1;
				else if ( nP < nP2 )
					m_smask[nP*m_ncams+nP2] = 1;
				else
				{
					m_smask[nP2*m_ncams + nP] = 1;
				}
					
			}
		}					
		feastart += nframes;
		ptno++;
	}


	// count number of non-zero element in S matrix
	m_nS = 0;
	for ( i = 0; i < m_ncams; i++ ) 
	{
		for (j = 0; j < m_ncams; j++)
		{
			if (m_smask[i*m_ncams + j] == 1)
			{
				m_nS++;
			}
		}
	}
}

void PBA::pba_readAndInitialize(char *camsfname, char *ptsfname,char *calibfname, int *ncams,
	int *n3Dpts, int *n2Dprojs,
	double **motstruct, double **imgpts,
	int **archor, char **vmask,
	char **umask, int **nphoto,
	int** nfeature, int** archorSort)
{
	FILE *fpc = nullptr, *fpp = nullptr, *fpXYZ = nullptr;
	int i, tmp1, tmp2;
	double ptMain[3], ptA[3];
	double dW1, dW2;	

	//calculate number of cameras, 3D points and projection points
	fopen_s(&fpc, camsfname, "r" );
	*ncams	=	findNcameras( fpc );
	m_ncams =	*ncams;
	m_V = (int*)malloc(sizeof(int) * m_ncams);

	fopen_s(&fpp, ptsfname, "r" );
	readNpointsAndNprojections( fpp, n3Dpts, 3, n2Dprojs, 2 );

	*motstruct = (double*)malloc((*ncams * 6 + *n3Dpts * 3) * sizeof(double));
	if(	*motstruct==NULL )
	{
		fprintf(stderr, "ParallaxBA error: Memory allocation for 'motstruct' failed \n");
		exit(1);
	}

	*imgpts = (double*)malloc(*n2Dprojs * 2 * sizeof(double));
	if(	*imgpts==NULL )
	{
		fprintf(stderr, "ParallaxBA error: Memory allocation for 'imgpts' failed\n");
		exit(1);
	}

	rewind(fpc);
	rewind(fpp);

	//allocate indicator of U
	*umask = (char*)malloc(*ncams * *ncams );
	memset(*umask, 0, *ncams * *ncams * sizeof(char));

	//allocate main and associate anchors
	*archor = (int*)malloc(*n3Dpts*3*sizeof(int));//
	memset(*archor, -1, *n3Dpts * 3 * sizeof(int)); 

	*nphoto		= (int*)malloc(*n2Dprojs*3*sizeof(int));//
	*nfeature	= (int*)malloc(*n2Dprojs*3*sizeof(int));//
	*archorSort = (int*)malloc(*n3Dpts*3*sizeof(int));

	ba_readCameraPose(fpc, *motstruct, m_V);


	fclose(fpc);
	fpc = NULL;
	
	//Update KR
	m_KR  = (double*)malloc(m_ncams*9*sizeof(double));
	m_KdA = (double*)malloc(m_ncams*9*sizeof(double));
	m_KdB = (double*)malloc(m_ncams*9*sizeof(double));
	m_KdG = (double*)malloc(m_ncams*9*sizeof(double));
	
	ba_updateKR(m_KR, m_K, *motstruct);

	//if XYZ are provided, we can use them as feature initialization.
	// fprintf(stdout, "%s\n", m_bProvideXYZ ? "true" : "false");  
	if (m_bProvideXYZ)
	{
		// printf("%s\n",m_szXYZ);
		fopen_s(&fpXYZ, m_szXYZ, "r");
		m_XYZ = (double*)malloc(m_n3Dpts*3*sizeof(double));

		for( i = 0; i < m_n3Dpts; i++){
			fscanf_s(fpXYZ, "%lf  %lf  %lf", m_XYZ + i * 3, m_XYZ + i * 3 + 1, m_XYZ + i * 3 + 2);

			Point3D p3d;
			p3d.xyz[0] = m_XYZ[i * 3];
			p3d.xyz[1] = m_XYZ[i * 3 + 1];
			p3d.xyz[2] = m_XYZ[i * 3 + 2];
			points.push_back(p3d);

		}
		// for (int i = 0; i < m_n3Dpts;i++){
		// 	printf("%f %f %f\n", points[i].xyz[0], points[i].xyz[1], points[i].xyz[2]);
		// }

		fclose(fpXYZ);

		string originalPath(m_szXYZ);
		size_t pos = originalPath.find_last_of("/\\");
		string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
		string init3D = parentPath + "/" + "XYZ.ply";
		//save features xyz
		if (init3D.c_str() != NULL)
		{
			FILE *fp = nullptr;
			fopen_s(&fp, init3D.c_str(), "w");
			fprintf(fp, "%s\n", "ply");
			fprintf(fp, "%s\n", "format ascii 1.0");
			fprintf(fp, "%s %d\n", "element vertex", m_n3Dpts);
			fprintf(fp, "%s\n", "property float x");
			fprintf(fp, "%s\n", "property float y");
			fprintf(fp, "%s\n", "property float z");
			fprintf(fp, "%s\n", "end_header");
			for (i = 0; i < m_n3Dpts; i++)
				fprintf(fp, "%lf %lf %lf\n", m_XYZ[i * 3], m_XYZ[i * 3 + 1], m_XYZ[i * 3 + 2]);
			fclose(fp);
		}
	}


	pba_readProjectionAndInitilizeFeature(fpp,
		*motstruct + *ncams * 6,
		*imgpts,
		*vmask,
		*ncams,
		*archor,
		*umask,
		*nphoto,
		*nfeature,
		*archorSort);

	fclose(fpp);

	
	int nCount = 0;
	double pti2k[3];
	int cur = 0;

	if (m_bProvideXYZ)
	{
		for (i = 0; i < m_n3Dpts; i++)
		{
			int nM = m_archor[i*3+1];
			int nN = m_archor[i*3+2];
			//printf("%d %d\n", nM, nN);


			ptMain[0] = *(*motstruct + nM*6 + 3);
			ptMain[1] = *(*motstruct + nM*6 + 4);
			ptMain[2] = *(*motstruct + nM*6 + 5);

			ptA[0] = *(*motstruct + nN*6 + 3);
			ptA[1] = *(*motstruct + nN*6 + 4);
			ptA[2] = *(*motstruct + nN*6 + 5);

			pti2k[0] = ptA[0] - ptMain[0];
			pti2k[1] = ptA[1] - ptMain[1];
			pti2k[2] = ptA[2] - ptMain[2];

			double dispti2k;
			dispti2k = sqrt(pti2k[0] * pti2k[0] + pti2k[1] * pti2k[1] + pti2k[2] * pti2k[2]);//����ģ��

			ptMain[0] = m_XYZ[i * 3 + 0] - ptMain[0];
			ptMain[1] = m_XYZ[i * 3 + 1] - ptMain[1];
			ptMain[2] = m_XYZ[i * 3 + 2] - ptMain[2]; 

			ptA[0] = m_XYZ[i * 3 + 0] - ptA[0];
			ptA[1] = m_XYZ[i * 3 + 1] - ptA[1];
			ptA[2] = m_XYZ[i * 3 + 2] - ptA[2];

			dW1 = ptMain[0] * ptMain[0] + ptMain[1] * ptMain[1] + ptMain[2] * ptMain[2];

			dW2 = ptA[0] * ptA[0] + ptA[1] * ptA[1] + ptA[2] * ptA[2];

			double disDot2;
			disDot2 = ptMain[0] * pti2k[0] + ptMain[1] * pti2k[1] + ptMain[2] * pti2k[2]; 
			double dww = disDot2 / (dispti2k * sqrt(dW1));


			double* pKR = m_KR + nM * 9;
			double n[2], n2[2], ptXj[3];

			ptXj[0] = ptMain[0];	ptXj[1] = ptMain[1];	ptXj[2] = ptMain[2];

			n[0] = (pKR[0]*ptXj[0] + pKR[1]*ptXj[1] + pKR[2]*ptXj[2])/
				(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

			n[1] = (pKR[3]*ptXj[0] + pKR[4]*ptXj[1] + pKR[5]*ptXj[2])/
				(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

			pKR = m_KR + nN*9;

			ptXj[0] = ptA[0];	ptXj[1] = ptA[1];	ptXj[2] = ptA[2];
			n2[0] = (pKR[0]*ptXj[0] + pKR[1]*ptXj[1] + pKR[2]*ptXj[2])/
				(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

			n2[1] = (pKR[3]*ptXj[0] + pKR[4]*ptXj[1] + pKR[5]*ptXj[2])/
				(pKR[6] * ptXj[0] + pKR[7] * ptXj[1] + pKR[8] * ptXj[2]);

			//printf("%d %d %d\n", m_archorSort[i * 2], m_archorSort[i * 2 + 1], m_archor[i * 3]);
			
			int id1 = cur + m_archorSort[i * 2];
			int id2 = cur + m_archorSort[i * 2 + 1];
			//printf("%f %f\n", m_imgpts[id1 * 2], m_imgpts[id1 * 2 + 1]);
			//printf("%f %f\n", m_imgpts[id2 * 2], m_imgpts[id2 * 2 + 1]);
			double err1 = (m_imgpts[id1 * 2] - n[0]) * (m_imgpts[id1 * 2] - n[0]) + (m_imgpts[id1 * 2 + 1] - n[1]) * (m_imgpts[id1 * 2 + 1] - n[1]);;
			double err2 = (m_imgpts[id2 * 2] - n2[0]) * (m_imgpts[id2 * 2] - n2[0]) + (m_imgpts[id2 * 2 + 1] - n2[1]) * (m_imgpts[id2 * 2 + 1] - n2[1]);

			//printf("%f %f\n", err1, err2);
			cur += m_archor[i * 3];

			if ((sqrt(dW1) / sqrt(dW2) > 30) || (err1 > err2) && (m_archor[3 * i] == 2))
			{
				nCount++;

				m_archor[i*3+1] = nN;
				m_archor[i*3+2] = nM;

				tmp1 = m_archorSort[i*2] ;
				tmp2 = m_archorSort[i*2+1] ;

				m_archorSort[i*2] = tmp2;
				m_archorSort[i*2+1] = tmp1;

				double dDAngle = atan2(ptA[0], ptA[2]);
				double dHAngle = atan2(ptA[1], sqrt(ptA[0] * ptA[0] + ptA[2] * ptA[2]));

				(*motstruct)[m_ncams * 6 + i * 3] = dDAngle;
				(*motstruct)[m_ncams * 6 + i * 3 + 1] = dHAngle;

				//double dwwDot = ptMain[0]*ptA[0] + ptMain[1]*ptA[1] + ptMain[2]*ptA[2];				
			}	
		}
	}
	double x,y,z;
	for (int i = 0; i < points.size();i++){
		x=points[i].xyz[0];
		y=points[i].xyz[1];
		z=points[i].xyz[2];
		
		points[i].world_depth[0]=atan2(x,z);
		points[i].world_depth[1]=atan2(y,sqrt(x*x+z*z));
		points[i].world_depth[2]=sqrt(x*x+y*y+z*z);

		points[i].world_inverse_depth[0]=points[i].world_depth[0];
		points[i].world_inverse_depth[1]=points[i].world_depth[1];
		points[i].world_inverse_depth[2]=1/points[i].world_depth[2];


		points[i].nM = m_archor[i * 3 + 1];
		points[i].nA = m_archor[i * 3 + 2];
		points[i].direction[0] = (*motstruct)[m_ncams * 6 + i * 3];
		points[i].direction[1] = (*motstruct)[m_ncams * 6 + i * 3 + 1];
		points[i].parallax = (*motstruct)[m_ncams * 6 + i * 3 + 2];

		int nM = points[i].nM;
		double dx = points[i].xyz[0] - cams[nM].camera_center[0];
		double dy = points[i].xyz[1] - cams[nM].camera_center[1];
		double dz = points[i].xyz[2] - cams[nM].camera_center[2];
		double d = sqrt(dx * dx + dy * dy + dz * dz);
		points[i].inverse_depth = 1 / d;
		points[i].depth = d;
		points[i].archored_xyz[0] = dx;
		points[i].archored_xyz[1] = dy;
		points[i].archored_xyz[2] = dz;

		points[i].archored_xy_inverse_z[0] = dx;
		points[i].archored_xy_inverse_z[1] = dy;
		points[i].archored_xy_inverse_z[2] = 1/dz;

		points[i].xy_inverse_z[0] = points[i].xyz[0];
		points[i].xy_inverse_z[1] = points[i].xyz[1];
		points[i].xy_inverse_z[2] = 1/points[i].xyz[2];
	}
	for (int i = 0; i < cams.size();i++){
		double R[9], q[4], axis_angle[3], t[3];
		double ey = cams[i].euler_angle[0];
		double ex = cams[i].euler_angle[1];
		double ez = cams[i].euler_angle[2];
		double c1 = cos(ey), c2 = cos(ex), c3 = cos(ez);
		double s1 = sin(ey), s2 = sin(ex), s3 = sin(ez);
		R[0]=c1*c3-s1*s2*s3;     R[1]=c2*s3;     R[2]=s1*c3+c1*s2*s3;
		R[3]=-c1*s3-s1*s2*c3;    R[4]=c2*c3;     R[5]=-s1*s3+c1*s2*c3;
		R[6]=-s1*c2;             R[7]=-s2;       R[8]=c1*c2;
		t[0] = -R[0] * cams[i].camera_center[0] - R[1] * cams[i].camera_center[1] - R[2] * cams[i].camera_center[2];
		t[1] = -R[3] * cams[i].camera_center[0] - R[4] * cams[i].camera_center[1] - R[5] * cams[i].camera_center[2];
		t[2] = -R[6] * cams[i].camera_center[0] - R[7] * cams[i].camera_center[1] - R[8] * cams[i].camera_center[2];
		cams[i].translation[0]=t[0];
		cams[i].translation[1]=t[1];
		cams[i].translation[2]=t[2];

		const double trace = R[0] + R[4] + R[8];
		if(trace >=0.0){
			double t = sqrt(trace + 1.0);
			q[0] = 0.5 * t;
			t = 0.5 / t;
			q[1] = (R[7] - R[5]) * t;
			q[2] = (R[2] - R[6]) * t;
			q[3] = (R[3] - R[1]) * t;
		}else{
			int i = 0;
			if(R[4]>R[0]){
				i = 1;
			}
			if(R[8]>R[4*i]){
				i = 2;
			}
			const int j = (i + 1) % 3;
			const int k = (j + 1) % 3;
			double t = sqrt(R[4 * i] - R[4 * j] - R[4 * k] + 1.0);
			q[i + 1] = 0.5 * t;
			t = 0.5 / t;
			q[0] = (R[3 * k + j] - R[3 * j + k]) * t;
			q[j + 1] = (R[3 * j + i] + R[3 * i + j]) * t;
			q[k + 1] = (R[3 * k + i] + R[3 * i + k]) * t;
		}
		cams[i].quat[0] = q[0];
		cams[i].quat[1] = q[1];
		cams[i].quat[2] = q[2];
		cams[i].quat[3] = q[3];
		// printf("%f %f %f %f\n", q[0], q[1], q[2], q[3]);

		const double sin_squared_theta = q[1]*q[1]+q[2]*q[2]+q[3]*q[3];
		if(sin_squared_theta>0.0){
			const double sin_theta = sqrt(sin_squared_theta);
			const double &cos_theta = q[0];
			const double two_theta = 2.0 * ((cos_theta < 0.0) ? atan2(-sin_theta, -cos_theta) : atan2(sin_theta, cos_theta));
			const double k = two_theta / sin_theta;
			axis_angle[0] = q[1] * k;
			axis_angle[1] = q[2] * k;
			axis_angle[2] = q[3] * k;
		}else{
			const double k = 2.0;
			axis_angle[0] = q[1] * k;
			axis_angle[1] = q[2] * k;
			axis_angle[2] = q[3] * k;
		}
		cams[i].axis_angle[0] = axis_angle[0];
		cams[i].axis_angle[1] = axis_angle[1];
		cams[i].axis_angle[2] = axis_angle[2];

		// 
		// Eigen::Vector3d physical_t(cams[i].translation[0], cams[i].translation[1], cams[i].translation[2]);
		// Eigen::Vector3d rv(cams[i].axis_angle[0], cams[i].axis_angle[1], cams[i].axis_angle[2]);
		// double theta = rv.norm();
		// Eigen::Matrix3d rv_hat;
		// rv_hat << 0, -rv.z(), rv.y(),
		// 		rv.z(), 0, -rv.x(),
		// 		-rv.y(), rv.x(), 0;
		// Eigen::Matrix3d J;
		// if(theta<1e-6){
		// 	J = Eigen::Matrix3d::Identity();
		// }else{
		// 	J = Eigen::Matrix3d::Identity() + (1 - cos(theta)) / (theta * theta) * rv_hat + (theta - sin(theta)) / (theta * theta * theta) * rv_hat * rv_hat;
		// }
		// Eigen::Vector3d rho = J.inverse() * physical_t;

		// cams[i].se3[0] = rho[0];
		// cams[i].se3[1] = rho[1];
		// cams[i].se3[2] = rho[2];
		cams[i].se3[0] = cams[i].quat[0];
		cams[i].se3[1] = cams[i].quat[1];
		cams[i].se3[2] = cams[i].quat[2];
		cams[i].se3[3] = cams[i].quat[3];
		cams[i].se3[4] = cams[i].translation[0];
		cams[i].se3[5] = cams[i].translation[1];
		cams[i].se3[6] = cams[i].translation[2];

		// for (int j = 0; j < 6; ++j){
		// 	printf("%f ", cams[i].se3[j]);
		// }
		// printf("\n");

		// double data[6] = {
		// 	rho[0], rho[1], rho[2],
		// 	cams[i].axis_angle[0],   cams[i].axis_angle[1],   cams[i].axis_angle[2]
		// };
		// Eigen::Map<const Eigen::Matrix<double, 6, 1>> se3_vec(data);
		// Sophus::SE3<double> T = Sophus::SE3<double>::exp(se3_vec);
		// Eigen::Matrix<double, 6, 1> xi = T.log();
		// for (int j = 0; j < 6; ++j){
		// 	cams[i].se3[j] = xi[j];
		// 	printf("%f ", cams[i].se3[j]);
		// }
		// printf("\n");
		// printf("%f %f %f\n", T.translation()[0], T.translation()[1], T.translation()[2]);
		// printf("%f %f %f\n", physical_t[0], physical_t[1], physical_t[2]);
		// cams[i].se3[0] = cams[i].translation[0];
		// cams[i].se3[1] = cams[i].translation[1];
		// cams[i].se3[2] = cams[i].translation[2];
	}
	

	if (m_bProvideXYZ)
		free(m_XYZ);

	fpXYZ = NULL;
}


bool PBA::pba_initializeMainArchor( 
	double* imgpts,	
	double* camera,	
	double* K,		
	double* feature,
	int nP,			
	int FID,		
	double* KR )	
{
	//solve  KRX = x
	Vector3d x;
	if (m_bProvideXYZ)
	{
		x(0) = m_XYZ[FID*3] - *(camera + nP*6+3);
		x(1) = m_XYZ[FID*3+1] - *(camera + nP*6+4);
		x(2) = m_XYZ[FID*3+2] - *(camera + nP*6+5);
	}
	
	else
	{
		double *ptr = m_KR + nP*9;	
		Matrix3d  A;
		A << ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7], ptr[8];
		
		double matx[3];				
		matx[0] = imgpts[0];
		matx[1] = imgpts[1];
		matx[2] = 1;
		
		Vector3d  b(matx);
		
		x = A.colPivHouseholderQr().solve(b);
	}

	double* pKR = KR + nP*9;
	double t = pKR[6]*x(0) + pKR[7]*x(1) + pKR[8]*x(2);	

	//compute azimuth and elevation angle
	double dDAngle = atan2( x(0), x(2) );				
	double dHAngle = atan2( x(1), sqrt(x(0)*x(0)+ x(2)*x(2)) );

	feature[0] = dDAngle;
	feature[1] = dHAngle;
	feature[2] = 0;

	if ( t < 0 )
		return true;
	else
		return false;
}

bool PBA::pba_initializeAssoArchor( 
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int nMI,
	int nAI,
	int FID,
	bool bLast )
{
	int nM = photo[nMI];                           
	int nA = photo[nAI];                         

	Vector3d  xM, xA;

	if (m_bProvideXYZ)
	{
		xM[0] = m_XYZ[FID*3]   - *(camera + nM*6+3);
		xM[1] = m_XYZ[FID*3+1] - *(camera + nM*6+4);
		xM[2] = m_XYZ[FID*3+2] - *(camera + nM*6+5);

		xA[0] = m_XYZ[FID*3]   - *(camera + nA*6+3);
		xA[1] = m_XYZ[FID*3+1] - *(camera + nA*6+4);
		xA[2] = m_XYZ[FID*3+2] - *(camera + nA*6+5);
	}
	else
	{
		//Main anchor ray
		double *ptr1 = m_KR + nM*9;
		Matrix3d  AM;	
		AM << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

		double matxM[3];
		matxM[0] = *(imgpts+2*nMI);
		matxM[1] = *(imgpts+2*nMI+1);
		matxM[2] = 1;

		Vector3d  bM(matxM);
		xM = AM.colPivHouseholderQr().solve(bM);			

		//Associate archor ray
		double *ptr2 = m_KR + nA*9;
		Matrix3d  AA;	
		AA << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

		double matxA[3];
		matxA[0] = *(imgpts+2*nAI);
		matxA[1] = *(imgpts+2*nAI+1);
		matxA[2] = 1;

		Vector3d  bA(matxA);
		xA = AA.colPivHouseholderQr().solve(bA);			
	}
		
	//Parallax Angle
	double dDot = xM(0)*xA(0) + xM(1)*xA(1) + xM(2)*xA(2);			

	double dDisM = sqrt( xM(0)*xM(0)+xM(1)*xM(1)+xM(2)*xM(2) );		
	double dDisA = sqrt( xA(0)*xA(0)+xA(1)*xA(1)+xA(2)*xA(2) );		

	if (dDot/(dDisM*dDisA)>1)			
		feature[2] = 0;
	else if (dDot/(dDisM*dDisA)<-1)		
		feature[2] = PI;
	else
	{
		double dw = acos( dDot/(dDisM*dDisA) );
		feature[2] = dw;
	}

	
	double pti2k[3];
	pti2k[0] = *(camera + nA*6+3) - *(camera + nM*6+3);    
	pti2k[1] = *(camera + nA*6+4) - *(camera + nM*6+4);    
	pti2k[2] = *(camera + nA*6+5) - *(camera + nM*6+5);    

	double dDot1 = xM[0]*pti2k[0] + xM[1]*pti2k[1] + xM[2]*pti2k[2];
	double dDisi2k = sqrt( pti2k[0]*pti2k[0] + pti2k[1]*pti2k[1] + pti2k[2]*pti2k[2] );
	double tmp = dDot1/(dDisM*dDisi2k);
	double dW2;
	if (tmp > 1)						
		dW2 = 0;
	if ( tmp < -1)						
		dW2 = PI;
	else
		dW2  = acos( tmp );

	return true;
}

bool PBA::pba_initializeOtheArchors( 
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int* archorSort,
	int nfeacout,
	int nOI,
	int FID )
{
	static int i = 0;
	double dw = feature[2];                   
	double dwNew;
	double dmaxw = dw;
	int   nNewI = 0;
	bool bAdjust = false;
	double dDot,dDisM,dDisA;

	if ( dw < MAXARCHOR   )
	{
		//current archor vector 
		int nO = photo[nOI];

		Vector3d  xO;
		if ( m_bProvideXYZ)
		{
			xO(0) = m_XYZ[FID*3] - *(camera + nO*6 + 3);
			xO(1) = m_XYZ[FID*3+1]-*(camera + nO*6 + 4);
			xO(2) = m_XYZ[FID*3+2]-*(camera + nO*6 + 5);
		}
		else
		{
			double *ptr1 = m_KR + nO*9;
			Matrix3d  AO;	
			AO << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

			double matxO[3];
			matxO[0] = *(imgpts+nOI*2);
			matxO[1] = *(imgpts+nOI*2+1);
			matxO[2] = 1;

			Vector3d  bO(matxO);
			xO = AO.colPivHouseholderQr().solve(bO);
		}

		double dDAngle = atan2( xO(0), xO(2) );
		double dHAngle = atan2( xO(1), sqrt(xO(0)*xO(0)+ xO(2)*xO(2)) );

		for ( i = 0; i < nfeacout; i++ )
		{
			//Main Archor Vector
			int nM = photo[i];                            
			Vector3d  xM;

			if (m_bProvideXYZ)
			{
				xM(0) = m_XYZ[FID*3]  -*(camera + nM*6+3);
				xM(1) = m_XYZ[FID*3+1]-*(camera + nM*6+4);
				xM(2) = m_XYZ[FID*3+2]-*(camera + nM*6+5);
			}
			else
			{
				double *ptr2 = m_KR + nM*9;
				Matrix3d  AM;	
				AM << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

				double matxM[3];
				matxM[0] = *(imgpts+i*2);
				matxM[1] = *(imgpts+i*2+1);
				matxM[2] = 1;

				Vector3d  bM(matxM);
				xM = AM.colPivHouseholderQr().solve(bM);
			}

			//Parallax angle between current archor and main archor
			dDot = xM(0)*xO(0) + xM(1)*xO(1) + xM(2)*xO(2);
			dDisM = sqrt( xM(0)*xM(0)+xM(1)*xM(1)+xM(2)*xM(2) );
			dDisA = sqrt( xO(0)*xO(0)+xO(1)*xO(1)+xO(2)*xO(2) );

			if( dDot/(dDisM*dDisA) > 1 )
				dwNew = 0;
			else if(dDot/(dDisM*dDisA)<-1)
				dwNew = PI;
			else
				dwNew = acos( dDot/(dDisM*dDisA) );

			if ( dwNew > dmaxw )
			{
				dmaxw = dwNew;
				archorSort[0] = nOI;   
				archorSort[1] = i;     
				feature[0] = dDAngle;
				feature[1] = dHAngle;
				feature[2] = dmaxw;
				bAdjust = true;
			}
		}
	}
	return bAdjust;
}



