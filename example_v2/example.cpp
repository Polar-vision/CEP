#include "BAExporter_v2.h"
#include "dataPath.h"
#include <vector>
#include <map>
#include <string>
using namespace std;

int main(int argc, char* argv[] )
{
	for(int i=0;i<9;i++){
		// if(i!=0&&i!=7&&i!=8){
		// 	continue;
		// }
		printf("hello ba!\n");
		// argv[1] = const_cast<char*>(dt);

		/*Parameterization of image point*/
		imagepointtype iptype;
		for(int j=0;j<1;j++){
			switch(j){
				case 0:iptype=uv;break;
				case 1:iptype=light_cone;break;
			}

			/*Parameterization of object point*/
			objectpointtype optype;
			switch(i){
				case 0:optype=xyz;break;
				case 1:optype=xy_inverse_z;break;
				case 2:optype=depth;break;
				case 3:optype=inverse_depth;break;
				case 4:optype=archored_xyz;break;
				case 5:optype=archored_xy_inverse_z;break;
				case 6:optype=archored_depth;break;
				case 7:optype=archored_inverse_depth;break;
				case 8:optype=parallax;break;
			}
			//zero archor
			// optype=xyz;
			// optype=xy_inverse_z;
			// optype=depth;
			// optype=inverse_depth;
			//one archor
			// optype=archored_xyz;
			// optype=archored_xy_inverse_z;
			// optype=archored_depth;
			// optype=archored_inverse_depth;
			//two archors
			// optype=parallax;
			/*Parameterization of 3d rotation*/
			rotation3dtype r3dtype;
			r3dtype=euler_angle;
			// r3dtype=axis_angle;
			// r3dtype=quaternion;
			
			parametertype paramtype;
			paramtype=rotation_translation_landmark;
			// paramtype=rotation_translation;
			// paramtype=rotation_landmark;
			// paramtype=translation_landmark;
			// paramtype=rotation;
			// paramtype=translation;
			// paramtype=landmark;

			manifoldtype manitype;
			// manitype=lie;
			// manitype=quaternion_manifold;
			// manitype=line_manifold;
			// manitype=sphere_manifold;
			// manitype=euclidean_manifold;
			manitype=none;

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

			}

			string pCheck = string(object_point_type)+"_"+rotation_3d_type+"_"+image_point_type;
			string originalPath(dt);
			size_t pos = originalPath.find_last_of("/\\");
			string parentPath = (pos != string::npos) ? originalPath.substr(0, pos) : "";
			string pP = parentPath + "/";
			string p1 = pP + "Cam.txt";//被噪声污染后的外参
			string p2 = pP + "Feature.txt";
			string p3 = pP + "XYZ.txt";//重新三角化后的物点
			string p4 = pP + "cal.txt";
			string pReport = "-report.txt";
			string pPose = "-FinalPose.txt";
			string p3D = "-Final3D.ply";
			string p5, p6, p7, pInit3D;
			pInit3D = pP + "XYZ.ply";

			p5 = pP + pCheck + pReport;
			p6 = pP + pCheck + pPose;
			p7 = pP + pCheck + p3D;

			char* szCam = const_cast<char*>(p1.c_str());
			char* szFea = const_cast<char*>(p2.c_str());
			//char* szXYZ = NULL;
			char* szXYZ = const_cast<char*>(p3.c_str());
			char* szCalib = const_cast<char*>(p4.c_str());
			char* szReport = const_cast<char*>(p5.c_str()); 
			char* szPose = const_cast<char*>(p6.c_str()); 
			char* sz3D = const_cast<char*>(p7.c_str()); 

			BAExporter BA;
			BA.ba_run(szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D, optype,r3dtype,iptype,paramtype,manitype);
		}
	}


	return 0;
}

