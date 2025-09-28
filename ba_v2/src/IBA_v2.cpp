
#include "IBA_v2.h"
#include <stdio.h>

//�����ļ��з�ע���е�����
int IBA::findNcameras(FILE* fp)
{
	int lineno, ncams, ch;

	lineno = ncams = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#') { /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);

		SKIP_LINE(fp);
		++lineno;
		if (ferror(fp))
		{
			fprintf(stderr, "findNcameras(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		++ncams;
	}
	return ncams;
}

IBA::IBA(void)
{
}

IBA::~IBA(void)
{
}

int IBA::countNDoubles(FILE* fp)
{
	int lineno, ch, np, i;
	char buf[MAXSTRLEN], * s;
	double dummy;

	lineno = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) return 0;

		ungetc(ch, fp);
		++lineno;
		if (!fgets(buf, MAXSTRLEN - 1, fp)) { /* read the line found... */
			fprintf(stderr, "countNDoubles(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		/* ...and count the number of doubles it has */
		for (np = i = 0, s = buf; 1; ++np, s += i) {
			ch = sscanf_s(s, "%lf%n", &dummy, &i);
			if (ch == 0 || ch == EOF) break;
		}

		rewind(fp);
		return np;
	}
	return 0; // should not reach this point
}

int IBA::skipNDoubles(FILE* fp, int nvals)
{
	int i;
	int j;

	for (i = 0; i < nvals; ++i)
	{
		j = fscanf_s(fp, "%*f");
		if (j == EOF) return EOF;

		if (ferror(fp)) return EOF - 1;
	}

	return nvals;
}

void IBA::readNpointsAndNprojections(FILE* fp, int* n3Dpts, int pnp, int* nprojs, int mnp)
{
	int nfirst, lineno, npts, nframes, ch, n;

	/* #parameters for the first line */
	nfirst = countNDoubles(fp);

	*n3Dpts = *nprojs = lineno = npts = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);
		++lineno;
		//skipNDoubles(fp, pnp);
		n = readNInts(fp, &nframes, 1);
		
		if (n != 1)
			exit(1);

		//printf("%d ", nframes);

		SKIP_LINE(fp);
		*nprojs += nframes;
		++npts;
	}

	*n3Dpts = npts;
}

void IBA::ba_readCablibration(FILE* fp, double* K)
{
	int n = fscanf_s(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf\n", &K[0], &K[3], &K[6], &K[1], &K[4], &K[7], &K[2], &K[5], &K[8]);

	if (n != 9)
	{
		fprintf(stderr, "BA error: Format of Calibaration is wrong\n");
		exit(1);
	}
}

void IBA::ba_readCameraPose(FILE* fp, double* params, int* m_v)
{
	int n, num, lineno = 0;
	double* tofilter;
	double* pPrams = params;
	int* pm_C = m_v;
	int jjg = 0;
	//the number of element per line is 8, it represents that focal length vary, or it is constant
	num = countNDoubles(fp);
	if (num == 8)
	{
		m_bFocal = true;
		m_K = (double*)malloc(m_ncams * 2 * sizeof(double));
		tofilter = (double*)malloc(8 * sizeof(double));
	}
	else if (num == 7) {
		tofilter = (double*)malloc(7 * sizeof(double));
	}
	else
		tofilter = (double*)malloc(6 * sizeof(double));

	while (!feof(fp))
	{
		if (num == 6) {
			n = readNDoubles(fp, tofilter, 6);
			if (n == -1) {
				//printf("%d %d\n", lineno, jjg);
				break;
			}
			Camera cam;
			cam.euler_angle[0] = tofilter[0];
			cam.euler_angle[1] = tofilter[1];
			cam.euler_angle[2] = tofilter[2];
			cam.camera_center[0] = tofilter[3];
			cam.camera_center[1] = tofilter[4];
			cam.camera_center[2] = tofilter[5];
			cam.camidx = 1;
			cams.push_back(cam);

			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2];
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5];
			*pm_C = 1;
			//m_v[lineno++] = 1;
		}
		if (num == 7){
			n = readNDoubles(fp, tofilter, 7);
			if (n == -1) {
				//printf("%d %d\n", lineno, jjg);
				break;
			}

			Camera cam;
			cam.euler_angle[0] = tofilter[0];
			cam.euler_angle[1] = tofilter[1];
			cam.euler_angle[2] = tofilter[2];
			cam.camera_center[0] = tofilter[3];
			cam.camera_center[1] = tofilter[4];
			cam.camera_center[2] = tofilter[5];
			cam.camidx = static_cast<int>(tofilter[6]);
			cams.push_back(cam);

			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2];
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5];
			*pm_C = static_cast<int>(tofilter[6]);
			++jjg;
		}
		if (num == 8){
			n = readNDoubles(fp, tofilter, 8);
			if (n == -1) {
				//printf("%d %d\n", lineno, jjg);
				break;
			}
			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2];
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5];

			m_K[lineno * 2] = tofilter[6];
			m_K[lineno * 2 + 1] = tofilter[7];
		}
			
		pPrams += 6;
		if (num == 7 || num == 6) {
			pm_C += 1;
		}
		++lineno;
	}
	if (tofilter != NULL) {
		free(tofilter);
		tofilter = NULL;
	}
}

int IBA::readNInts(FILE* fp, int* vals, int nvals)
{
	int i;
	int n, j;

	for (i = n = 0; i < nvals; ++i) {//��һ�Σ�nvals����1
		j = fscanf_s(fp, "%d", vals + i);//��һ�Σ���ȡ��һ����Ҳ����1����ά�������2DͶӰ����
		if (j == EOF) return EOF;//��һ�Σ�һ�㲻���ǣ�����j=1����ʾ�ɹ�

		if (j != 1 || ferror(fp)) return EOF - 1;//��һ�Σ�һ�㶼�ǣ�j=1����ʾ�ɹ�

		n += j;//��һ����Ϊn=0�����Ծ���j��Ҳ���ǳɹ���־j=1
	}//֮��Ļ���fscanfÿ��һ����ô�ͻ��ƶ�ָ�룬Ȼ��%d��������ո񡢻��з��������հ��ַ���ֱ��������%d����������
	//����������EOF����j!=0��
	return n;
}

int IBA::readNDoubles(FILE* fp, double* vals, int nvals)
{
	int i;
	int n, j;

	for (i = n = 0; i < nvals; ++i)//nvals��2����ȡxy����
	{
		j = fscanf_s(fp, "%lf", vals + i);//���￪ʼ������float��Ҳ��������ֵ
		if (j == EOF) return EOF;

		if (j != 1 || ferror(fp)) return EOF - 1;

		n += j;//����֮���2��Ϊ����2��
	}

	return n;//����2
}

void IBA::ba_readCameraPoseration(char* fname, double* ical)
{
	FILE* fp = nullptr;
	int  ch = EOF;
	fopen_s(&fp, fname, "r");
	if (fp == nullptr)
	{
		fprintf(stderr, "BA: Cannot open calbration file %s, exiting\n", fname);
		return;
	}

	int s = 0;
	for (int i = 0; i < nc_; i++) {
		int num = fscanf_s(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", \
			& ical[9 * i + 0], &ical[9 * i + 3], &ical[9 * i + 6], \
			& ical[9 * i + 1], &ical[9 * i + 4], &ical[9 * i + 7], \
			& ical[9 * i + 2], &ical[9 * i + 5], &ical[9 * i + 8]);
		Intrinsic intr;
		intr.fx = ical[9 * i + 0];
		intr.fy = ical[9 * i + 4];
		intr.cx = ical[9 * i + 6];
		intr.cy = ical[9 * i + 7];
		intrs.push_back(intr);
		s += num;
	}
	if (s != 9 * nc_)
	{
		fprintf(stderr, "BA error: Format of Calibration file is wrong");
		return;
	}

	fclose(fp);
	fp = NULL;
}

void IBA::ba_updateKR(double* KR, double* K, double* p)
{
	if (!m_bFocal)
	{
		int i = 0;
		double* ptAngle;
		double* pKR, * pK;
		double matR[9];//
		double matRG[9], matRB[9], matRA[9];
		//for (i = 0; i < m_ncams; i++) {
		//	printf("%d\n", m_V[i]);
		//}
		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;//ָ
			/*phi omega kappaϵͳ*/
			//matR=matRG*matRB*matRA
			//ptAngle=[kappa,phi,omega]
			double ey = ptAngle[0];
			double ex = ptAngle[1];
			double ez = ptAngle[2];
			double c1 = cos(ey);   double c2 = cos(ex);   double c3 = cos(ez);
			double s1 = sin(ey);   double s2 = sin(ex);   double s3 = sin(ez);
			matR[0]=c1*c3-s1*s2*s3;     matR[1]=c2*s3;     matR[2]=s1*c3+c1*s2*s3;
			matR[3]=-c1*s3-s1*s2*c3;    matR[4]=c2*c3;     matR[5]=-s1*s3+c1*s2*c3;
			matR[6]=-s1*c2;             matR[7]=-s2;       matR[8]=c1*c2;
			// matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			// matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			// matR[2] = -sin(ptAngle[1]);
			// matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			// matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			// matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			// matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			// matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			// matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			//omega
			matRG[0] = 1;		matRG[1] = 0;				matRG[2] = 0;
			matRG[3] = 0;		matRG[4] = cos(ptAngle[2]);	matRG[5] = sin(ptAngle[2]);
			matRG[6] = 0;		matRG[7] = -sin(ptAngle[2]);	matRG[8] = cos(ptAngle[2]);

			//phi
			matRB[0] = cos(ptAngle[1]);		matRB[1] = 0;		matRB[2] = -sin(ptAngle[1]);
			matRB[3] = 0;					matRB[4] = 1;		matRB[5] = 0;
			matRB[6] = sin(ptAngle[1]);		matRB[7] = 0;		matRB[8] = cos(ptAngle[1]);

			//kappa
			matRA[0] = cos(ptAngle[0]);		matRA[1] = sin(ptAngle[0]);			matRA[2] = 0;
			matRA[3] = -sin(ptAngle[0]);	matRA[4] = cos(ptAngle[0]);			matRA[5] = 0;
			matRA[6] = 0;					matRA[7] = 0;						matRA[8] = 1;

			//pKR=KR*matR
			pKR = KR + i * 9;
			pK = K + (m_V[i] - 1) * 9;//
			//printf("%d\n", m_V[i]);
			//printf("%f %f %f\n", pK[0], pK[3], pK[6]);
			//printf("%f %f %f\n", pK[1], pK[4], pK[7]);
			//printf("%f %f %f\n", pK[2], pK[5], pK[8]);
			pKR[0] = pK[0] * matR[0] + pK[3] * matR[3] + pK[6] * matR[6];
			pKR[1] = pK[0] * matR[1] + pK[3] * matR[4] + pK[6] * matR[7];
			pKR[2] = pK[0] * matR[2] + pK[3] * matR[5] + pK[6] * matR[8];
			pKR[3] = pK[1] * matR[0] + pK[4] * matR[3] + pK[7] * matR[6];
			pKR[4] = pK[1] * matR[1] + pK[4] * matR[4] + pK[7] * matR[7];
			pKR[5] = pK[1] * matR[2] + pK[4] * matR[5] + pK[7] * matR[8];
			pKR[6] = pK[2] * matR[0] + pK[5] * matR[3] + pK[8] * matR[6];
			pKR[7] = pK[2] * matR[1] + pK[5] * matR[4] + pK[8] * matR[7];
			pKR[8] = pK[2] * matR[2] + pK[5] * matR[5] + pK[8] * matR[8];
		}
	}
}

void IBA::readNpointsAndNprojectionsFromProj(FILE* fp, int& n3Dpts, int& nprojs)
{
	int nfirst, lineno, npts, nframes, ch, n;
	nprojs = 0;
	n3Dpts = 0;
	npts = 0;

	/* #parameters for the first line */
	nfirst = countNDoubles(fp);

	//*n3Dpts=*nprojs=lineno=npts=0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);
		++lineno;
		//skipNDoubles(fp, pnp);
		n = readNInts(fp, &nframes, 1);
		if (n != 1)
		{
			fprintf(stderr, "readNpointsAndNprojections(): error reading input file, line %d: "
				"expecting number of frames for 3D point\n", lineno);
			exit(1);
		}

		SKIP_LINE(fp);
		nprojs += nframes;
		++npts;
	}

	n3Dpts = npts;
}

void IBA::readPointProjections(FILE* fp, double* imgpts, int* photo, int* imgptsSum, int n3Dpts, int n2Dprojs)
{
	int nframes, ch, lineno, ptno, frameno, n;
	int i;
	int nproj2D = 0;

	lineno = ptno = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			lineno++;

			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);

		n = readNInts(fp, &nframes, 1);  /* read in number of image projections */
		if (n != 1)
		{
			fprintf(stderr, "sba_readProjectionAndInitilizeFeature(): error reading input file, line %d:\n"
				"expecting number of frames for 3D point\n", lineno);
			exit(1);
		}

		imgptsSum[ptno] = nframes;

		for (i = 0; i < nframes; ++i)
		{
			n = readNInts(fp, &frameno, 1); /* read in frame number... */

			photo[nproj2D] = frameno;

			n += readNDoubles(fp, imgpts + nproj2D * 2, 2); /* ...and image projection */

			nproj2D++;
		}
		fscanf_s(fp, "\n"); // consume trailing newline

		lineno++;
		ptno++;
	}
}
void IBA::readImagePts(const char* szProj, double** imgpts, int** photo, int** imgptsSum, int& n3Dpts, int& n2Dprojs)
{
	FILE* fpp = nullptr;
	fopen_s(&fpp, szProj, "r");
	if (fpp == nullptr) {
		fprintf(stderr, "cannot open file %s, exiting\n", szProj);
		exit(1);
	}
	readNpointsAndNprojectionsFromProj(fpp, n3Dpts, n2Dprojs);

	*imgpts = (double*)malloc(n2Dprojs * 2 * sizeof(double));
	if (*imgpts == NULL) {
		fprintf(stderr, "memory allocation for 'imgpts' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	*photo = (int*)malloc(n2Dprojs * sizeof(int));
	if (*photo == NULL)
	{
		fprintf(stderr, "memory allocation for 'struct' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	*imgptsSum = (int*)malloc(n3Dpts * sizeof(int));
	if (*imgptsSum == NULL)
	{
		fprintf(stderr, "memory allocation for 'struct' failed in readInitialSBAEstimate()\n");
		exit(1);
	}

	rewind(fpp);
	readPointProjections(fpp, *imgpts, *photo, *imgptsSum, n3Dpts, n2Dprojs);

	fclose(fpp);
}


