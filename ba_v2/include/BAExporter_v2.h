#include <cstddef>
#include <iostream>
#ifndef _BA_H
#define _BA_H

#ifdef BA_EXPORTS
    #define BAapi __declspec(dllexport)
#elif defined(BA_STATIC)
    #define BAapi  // 静态库无修饰符
#else
    #define BAapi __declspec(dllimport)
#endif

typedef enum ObjectPointType{
	//zero archors
	xyz=0,
	xy_inverse_z=1,
	depth=2,
	inverse_depth=3,
	//one archors
	archored_xyz=4,
	archored_xy_inverse_z=5,
	archored_depth=6,
	archored_inverse_depth=7,
	//two archors
	parallax=8
}objectpointtype;

typedef enum Rotation3DType{
	euler_angle=0,
	axis_angle=1,
	quaternion=2
}rotation3dtype;

typedef enum ImagePointType{
	uv=0,
	light_cone=1
}imagepointtype;

typedef enum ParameterType{
	rotation_translation_landmark=0,
	rotation_landmark=1,
	translation_landmark=2,
	rotation_translation=3,
	landmark=4,
	rotation=5,
	translation=6
}parametertype;

typedef enum ManifoldType{
	lie=0,
	quaternion_manifold=1,
	euclidean_manifold=2,
	sphere_manifold=3,
	line_manifold=4,
	none=5
}manifoldtype;

class IBA;//PBA和SBA的基类
class BAapi BAExporter
{
public:
	bool ba_run(char *szCam = NULL,
				char *szFea = NULL,
				char *szXYZ = NULL,
				char *szCalib = NULL,
				char *szReport = NULL,
				char *szPose = NULL,
				char *sz3D = NULL,
				objectpointtype optype = xyz,
				rotation3dtype r3dtype = euler_angle,
				imagepointtype iptype = uv,
				parametertype paramtype = rotation_translation_landmark,
				manifoldtype manitype = none);

	bool ba_initialize( char* szCamera, char* szFeature, char* szCalib =  NULL, char* szXYZ = NULL );


	BAExporter();
	~BAExporter();

private:
	IBA* ptr;
};

#endif