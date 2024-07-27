#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pro_fold11.h"

static inline double norm2(double u[3])
	/*returns the norm of a 3 dimensional vector*/
	{
	return u[0]*u[0]+u[1]*u[1]+u[2]*u[2];
	}
	
	
static inline double dot(double u[3],double v[3])
	/*returns the dot product of two 3-D vectors*/
	{
	return u[0]*v[0]+u[1]*v[1]+u[2]*v[2];
	}
	
	
double normalize(double u[3])
	/*normalizes u and changes its values to unit vector along u*/
	{
	int k;
	double n;
	
	n=sqrt(norm2(u));
	
	for(k=0;k<3;++k)
		u[k]/=n;
	
	return n;
	}
	
	
double perpnormalize(double v[3],double u[3])
	/*projects v on u and returns the norm of projected v or v_perp but
	 it also changes the values of array v to unit vector along  v_perp
	and it also normalizes u in the process changing u to u cap*/
	{
	int k;
	double d;
	
	d=dot(v,u);
	
	for(k=0;k<3;++k)
		v[k]-=d*u[k];
		
	return normalize(v);
	}
	
	
void mult(double a[3][3],double x[3],double ax[3])
	/*multiplies matrix a (3*3) with a 3-D vector x and stores it in ax vector*/
	{
	int i,k;
	
	for(i=0;i<3;++i)
		{
		ax[i]=0.;
		for(k=0;k<3;++k)
			ax[i]+=a[i][k]*x[k];
		}
	}
	
	
void cross(double u[3],double v[3],double uv[3])
	/*stores the cross product of u and v in uv*/
	{
	uv[0]=u[1]*v[2]-u[2]*v[1];
	uv[1]=u[2]*v[0]-u[0]*v[2];
	uv[2]=u[0]*v[1]-u[1]*v[0];
	}
	
	
int axes(double sym[3][3],double u[3],double v[3])
	/* to initialize u and v, u is chosen as the row of sym
	with the highest norm and v as the row of second highest 
	norm. After that, in an iteration, u is updated with sym*u
	and v is updated with v perp to u. If at any point, norm of
	u is less than v then they are exchanged with each other.*/
	{
	int i,j,k,iter;
	int r,iu[8]={2,1,3,1,2,3,0,0},iv[8]={1,2,3,0,0,3,2,1};
	double n[3],su[3],sv[3],e,f,tmp,dv;
	double DVMAX = 0.00000000000000000000000001;

	for(i=0;i<3;++i)
		n[i]=norm2(sym[i]);
		
	r=4*(n[0]>n[1])+2*(n[0]>n[2])+(n[1]>n[2]);
		
	for(k=0;k<3;++k)
		{
		u[k]=sym[iu[r]][k];
		v[k]=sym[iv[r]][k];
		}
		
	e=normalize(u);
	f=perpnormalize(v,u);

	for(iter=0;iter<totiter;++iter)
		{
		mult(sym,u,su);
		mult(sym,v,sv);
		
		e=normalize(su);
		f=normalize(sv);
		
		if(e<f)
			{
			tmp=e;
			e=f;
			f=tmp;
			
			for(k=0;k<3;++k)
				{
				tmp=su[k];
				su[k]=sv[k];
				sv[k]=tmp;
				}
			}
			
		for(k=0;k<3;++k)
			u[k]=su[k];
		
		perpnormalize(sv,u);
		
		dv=sqrt(1.-dot(sv,v));
		
		for(k=0;k<3;++k)
			v[k]=sv[k];
		
		//if(dv<DVMAX)
		//	break;
		}
		
	return iter;
	}
	
	
void rotproj(double in[3][3],double out[3][3])
/*out gives orthogonal matrix which maximises tr(out*in). To determine out, 
kabsch-Umeyama algorithm is used, which uses the left and right singular
vectors. To determine the left singular vectors, in*in is diagonalized whose 
eigenvectors(stored in rows of rot1) are the left singular vectors of in.
To determine the right singular vectors, rot1*in gives the right singular vectors
multiplied with the corresponding singular values of in. So, they are normalized and 
stored in rot2. To fix the determinant of rot2 as 1, the third singular vector is 
determined by taking the cross product of the first two. Hence,
out would be rot1^T*rot2 */
	{
	int i,j,k;
	double sym[3][3],rot1[3][3],rot2[3][3];
	
	for(i=0;i<3;++i)
	for(j=0;j<=i;++j)
		{
		sym[i][j]=0.;
		for(k=0;k<3;++k)
			sym[i][j]+=in[i][k]*in[j][k];
		}
		
	sym[0][1]=sym[1][0];
	sym[0][2]=sym[2][0];
	sym[1][2]=sym[2][1];
	
	axes(sym,rot1[0],rot1[1]);
	
	cross(rot1[0],rot1[1],rot1[2]);
	
	for(i=0;i<2;++i)
	for(j=0;j<3;++j)
		{
		rot2[i][j]=0.;
		for(k=0;k<3;++k)
			rot2[i][j]+=rot1[i][k]*in[k][j];
		}
		
	normalize(rot2[0]);
	normalize(rot2[1]);
	
	cross(rot2[0],rot2[1],rot2[2]);
	
	for(i=0;i<3;++i)
	for(j=0;j<3;++j)
		{
		out[i][j]=0.;
		for(k=0;k<3;++k)
			out[i][j]+=rot1[k][i]*rot2[k][j];
		}
	}
	
double RMSD(int n, double (*motif1)[3],double (*motif2)[3]){
	int i,k;
	double row_norm = 0.;
	for(i=0;i<n;++i){
		for(k=0;k<3;++k){
			row_norm += (motif1[i][k]-motif2[i][k])*(motif1[i][k]-motif2[i][k]);
			}
		}
	return sqrt(row_norm/n);	
	}
	
double NormDiff(int n, double (*motif1)[3],double (*motif2)[3]){
	int i,k;
	double row_norm = 0.;
	for(i=0;i<n;++i){
		for(k=0;k<3;++k){
			row_norm += (motif1[i][k]-motif2[i][k])*(motif1[i][k]-motif2[i][k]);
			}
		}
	return row_norm;	
	}

void center(int n, double (*x)[3],double (*xc)[3],double c[3])
	/*n is number of datapoints. x contains the cooradinates
	of those datapoints. the mean is calculated first and then
	it is subtracted from the coordinates which gives xc centered
	at the origin.*/
	{
	int i,k;
	
	for(k=0;k<3;++k)
		{
		c[k]=0.;
		
		for(i=0;i<n;++i)
			c[k]+=x[i][k];
			
		c[k]/=n;
		
		for(i=0;i<n;++i)
			xc[i][k]=x[i][k]-c[k];
		}
	}
	
	
void motifproj(int n, double (*motif)[3],double (*in)[3],double (*out)[3])
	{
	int k,l,i;
	double c[3],rot[3][3],rotmotif[3];

	center(n,in,out,c);//c contains the center of dataset in
	
	//* out is the dataset centered at origin on which we are projecting the motif *//
	// rot is the [3][3] matrix M^T from KU paper. Q is motif.transpose,
	// P is out.transpose. M = Q(P^T)
	for(k=0;k<3;++k)
	for(l=0;l<3;++l)
		{
		rot[k][l]=0.;
		
		for(i=0;i<n;++i)
			rot[k][l]+=out[i][k]*motif[i][l];
		}
	
	rotproj(rot,rot);
	
	for(i=0;i<n;++i)
		{
		mult(rot,motif[i],rotmotif);
		
		for(k=0;k<3;++k)
			out[i][k]=rotmotif[k]+c[k];
		}
	}
	
/*	
int main(int argc,char* argv[])
	{
	char *motiffile,*coordfile;
	FILE *fp;
	int i,k;
	double motif[6][3],coord[6][3],coordproj[6][3];
	
	if(argc==3)
		{
		motiffile=argv[1];
		coordfile=argv[2];
		}
	else
		{
		printf("expected two arguments: motiffile, coordfile\n");
		return 1;
		}
		
	fp=fopen(motiffile,"r");
	
	for(i=0;i<6;++i)
	for(k=0;k<3;++k)
		fscanf(fp,"%lf",&motif[i][k]);
		
	fclose(fp);
	
	fp=fopen(coordfile,"r");
	
	for(i=0;i<6;++i)
	for(k=0;k<3;++k)
		fscanf(fp,"%lf",&coord[i][k]);
		
	fclose(fp);
		
	motifproj(6,motif,coord,coordproj);
	
	for(i=0;i<6;++i)
		{
		for(k=0;k<3;++k)
			printf("%12.6lf",coordproj[i][k]);
		
		printf("\n");
		}
		
	return 0;
	}
*/
