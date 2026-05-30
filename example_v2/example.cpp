#include "BAExporter_v2.h"
#include "dataPath.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <map>
#include <string>
using namespace std;
namespace fs = std::filesystem;

struct RunRecord {
	int index = 0;
	string anchor_mode;
	string method_name;
	string code_name;
	string report_path;
	string convergence_path;
	string pose_path;
	string points_path;
	BAResult result;
};

const char* objectPointCodeName(objectpointtype optype)
{
	switch (optype)
	{
	case xyz: return "xyz";
	case xy_inverse_z: return "xy_inverse_z";
	case depth: return "depth";
	case inverse_depth: return "inverse_depth";
	case archored_xyz: return "archored_xyz";
	case archored_xy_inverse_z: return "archored_xy_inverse_z";
	case archored_depth: return "archored_depth";
	case archored_inverse_depth: return "archored_inverse_depth";
	case parallax: return "parallax";
	default: return "unknown";
	}
}

const char* comparisonMethodName(objectpointtype optype)
{
	switch (optype)
	{
	case xyz: return "xyz";
	case xy_inverse_z: return "inverse_depth";
	case depth: return "spherical";
	case inverse_depth: return "inverse_distance";
	case archored_xyz: return "anchored_xyz";
	case archored_xy_inverse_z: return "anchored_inverse_depth";
	case archored_depth: return "anchored_spherical";
	case archored_inverse_depth: return "anchored_inverse_distance";
	case parallax: return "parallax";
	default: return "unknown";
	}
}

const char* anchorModeName(objectpointtype optype)
{
	switch (optype)
	{
	case xyz:
	case xy_inverse_z:
	case depth:
	case inverse_depth:
		return "zero_anchor";
	case archored_xyz:
	case archored_xy_inverse_z:
	case archored_depth:
	case archored_inverse_depth:
		return "single_anchor";
	case parallax:
		return "double_anchor";
	default:
		return "unknown";
	}
}

const char* terminationName(int termination_type)
{
	switch (termination_type)
	{
	case 0: return "CONVERGENCE";
	case 1: return "NO_CONVERGENCE";
	case 2: return "FAILURE";
	case 3: return "USER_SUCCESS";
	case 4: return "USER_FAILURE";
	default: return "UNKNOWN";
	}
}

string csvEscape(const string& value)
{
	string escaped = value;
	size_t pos = 0;
	while ((pos = escaped.find('"', pos)) != string::npos)
	{
		escaped.insert(pos, 1, '"');
		pos += 2;
	}
	return "\"" + escaped + "\"";
}

double costDropPercent(const BAResult& result)
{
	if (result.initial_cost <= 0.0)
	{
		return 0.0;
	}
	return 100.0 * (result.initial_cost - result.final_cost) / result.initial_cost;
}

string formatDouble(double value, int precision = 6)
{
	ostringstream oss;
	oss << setprecision(precision) << value;
	return oss.str();
}

void writeComparisonCsv(const fs::path& csv_path, const vector<RunRecord>& records)
{
	ofstream csv(csv_path);
	csv << "index,anchor_mode,method,code_name,rotation,image,parameter,manifold,"
		<< "num_cameras,num_points,num_observations,iterations,successful_steps,"
		<< "unsuccessful_steps,num_linear_solves,initial_cost,final_cost,"
		<< "initial_rms_px,final_rms_px,cost_drop_percent,total_time_sec,"
		<< "minimizer_time_sec,linear_solver_time_sec,residual_evaluation_time_sec,"
		<< "jacobian_evaluation_time_sec,termination,report_file,convergence_file,"
		<< "pose_file,points_file\n";

	csv << setprecision(17);
	for (const auto& record : records)
	{
		const BAResult& r = record.result;
		csv << record.index << ","
			<< csvEscape(record.anchor_mode) << ","
			<< csvEscape(record.method_name) << ","
			<< csvEscape(record.code_name) << ","
			<< csvEscape("euler_angle") << ","
			<< csvEscape("uv") << ","
			<< csvEscape("rotation_translation_landmark") << ","
			<< csvEscape("none") << ","
			<< r.num_cameras << ","
			<< r.num_points << ","
			<< r.num_observations << ","
			<< r.num_iterations << ","
			<< r.num_successful_steps << ","
			<< r.num_unsuccessful_steps << ","
			<< r.num_linear_solves << ","
			<< r.initial_cost << ","
			<< r.final_cost << ","
			<< r.initial_rms << ","
			<< r.final_rms << ","
			<< costDropPercent(r) << ","
			<< r.total_time_sec << ","
			<< r.minimizer_time_sec << ","
			<< r.linear_solver_time_sec << ","
			<< r.residual_evaluation_time_sec << ","
			<< r.jacobian_evaluation_time_sec << ","
			<< csvEscape(terminationName(r.termination_type)) << ","
			<< csvEscape(fs::path(record.report_path).filename().string()) << ","
			<< csvEscape(fs::path(record.convergence_path).filename().string()) << ","
			<< csvEscape(fs::path(record.pose_path).filename().string()) << ","
			<< csvEscape(fs::path(record.points_path).filename().string()) << "\n";
	}
}

void writeComparisonMarkdown(const fs::path& md_path, const vector<RunRecord>& records)
{
	ofstream md(md_path);
	md << "# BA Method Comparison\n\n";
	md << "| # | Anchor | Method | Code | Iter | Time(s) | Initial RMS(px) | Final RMS(px) | Cost Drop(%) | Termination |\n";
	md << "| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |\n";

	for (const auto& record : records)
	{
		const BAResult& r = record.result;
		md << "| " << record.index
			<< " | " << record.anchor_mode
			<< " | " << record.method_name
			<< " | " << record.code_name
			<< " | " << r.num_iterations
			<< " | " << formatDouble(r.total_time_sec)
			<< " | " << formatDouble(r.initial_rms)
			<< " | " << formatDouble(r.final_rms)
			<< " | " << formatDouble(costDropPercent(r))
			<< " | " << terminationName(r.termination_type)
			<< " |\n";
	}

	md << "\nFull per-method Ceres reports are written next to this table as `*-report.txt` files.\n";
}

void writeComparisonTables(const fs::path& dataset_dir, const vector<RunRecord>& records)
{
	const fs::path csv_path = dataset_dir / "BA-comparison.csv";
	const fs::path md_path = dataset_dir / "BA-comparison.md";
	writeComparisonCsv(csv_path, records);
	writeComparisonMarkdown(md_path, records);
	printf("BA comparison CSV: %s\n", csv_path.string().c_str());
	printf("BA comparison Markdown: %s\n", md_path.string().c_str());
}

int main(int argc, char* argv[] )
{
	fs::path executableDir = fs::absolute(argv[0]).parent_path();
	fs::path dataPath = fs::path(dt);
	if (dataPath.is_relative()) {
		dataPath = executableDir / dataPath;
	}
	dataPath = fs::weakly_canonical(dataPath);
	if (!fs::exists(dataPath)) {
		printf("Missing data file: %s\n", dataPath.string().c_str());
		return 1;
	}

	const fs::path datasetDir = dataPath.parent_path();
	vector<RunRecord> records;

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
			case none:
				manifold_type="none";break;

			}

			string pCheck = string(object_point_type)+"_"+rotation_3d_type+"_"+image_point_type;
			string parentPath = datasetDir.string();
			string pP = parentPath + "/";
			string p1 = pP + "Cam.txt";// noisy initial camera poses
			string p2 = pP + "Feature.txt";
			string p3 = pP + "XYZ.txt";// retriangulated object points
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

			BAResult result;
			BAExporter BA;
			BA.ba_run(szCam, szFea, szXYZ, szCalib, szReport, szPose, sz3D, optype,r3dtype,iptype,paramtype,manitype,&result);

			RunRecord record;
			record.index = i;
			record.anchor_mode = anchorModeName(optype);
			record.method_name = comparisonMethodName(optype);
			record.code_name = objectPointCodeName(optype);
			record.report_path = p5;
			record.convergence_path = pP + "convergence_" + pCheck + ".txt";
			record.pose_path = p6;
			record.points_path = p7;
			record.result = result;
			records.push_back(record);
		}
	}

	writeComparisonTables(datasetDir, records);

	return 0;
}

