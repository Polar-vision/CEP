#include "BAExporter_v2.h"
#include "PBAImp_v2.h"

BAExporter::BAExporter(){
	ptr = new PBA;
}

BAExporter::~BAExporter(){
	delete ptr;
}

bool BAExporter::ba_run(char* szCam, char* szFea, char* szXYZ, char* szCalib, char* szReport,
	char* szPose, char* sz3D, objectpointtype optype, rotation3dtype r3dtype,imagepointtype iptype,
						parametertype paramtype,manifoldtype manitype){
	return ptr->ba_run(szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D, optype, r3dtype,iptype,paramtype,manitype);
}

bool BAExporter::ba_initialize(char* szCamera, char* szFeature, char* szCalib, char* szXYZ){
	return	ptr->ba_initialize(szCamera, szFeature, szCalib, szXYZ);
}
