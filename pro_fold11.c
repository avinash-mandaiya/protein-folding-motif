#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "pro_fold11.h"

char *ProID;
int ProLen,*RT,SizeWM,AA[91],lenAA[21],(*ATR[21]),NumHPAA,*HPAA,NumSPAA,*SPAA,iter,(*MatchArray)[2],NumSol,CHstep,totiter,totHBB,totHBS,NumSPatoms;
short int BBMotifNum[21][21][21][21][11], SolMotifNum[21][21][21], RotMotifNum[3][21];
double (*BBMotifs[21][21][21][21][11])[20][3], (**SolMotifs[21][21][21])[3], (**RotMotifs[3][21])[3];
double (**BBStore[2])[20][3], (***SolStore[2])[3];
double ATDMIN[4][4], ASDMIN[4], OODMIN, (**ActualCoor)[3], (*ActualSol)[3];

double (**BB[2])[20][3], (**BBA[2])[20][3], (**BBR[2])[20][3], (**BBB[2])[20][3];
double (**RS)[3], (**RSA)[3], (**RSR)[3], (**RSB)[3];
double (***SB[2])[3], (***SBA[2])[3], (***SBR[2])[3], (***SBB[2])[3];
double (****VB)[3], (****VBA)[3], (****VBB)[3], (****VBR)[3];
double (**VO)[3], (**VOA)[3], (**VOR)[3], (**VOB)[3];
double (***VS)[2][3], (***VSA)[2][3], (***VSB)[2][3], (***VSR)[2][3];
double (**atom)[3], (*Solatom)[3];
double tBBerr, tVXerr, tSBerr, tRSerr, toterr;
double beta, **etaBB[2], ****etaVB, ***etaVS, **etaVO, *etaSB[2], *etaRS, epsilon, *CAerr, fracSP, fracOS;
double **WeightsBB[2], **WeightsSB[2], (*Wrank)[2],(*WrankB)[3],(*WrankS)[3];

char errfile[50],statsfile[50],solfile[50],ProFile[50],CAfile[50];

int **SBMatch,**BBMatch;


int compare_error_asc_G( const void *pa, const void *pb)
        {    
	const double *a = pa;
	const double *b = pb;

        if ( *a < *b ) return -1;
        if ( *a > *b ) return +1;
        return 0;
        }


int compare_error_asc( const void *pa, const void *pb )
        {    
        const double (*a)[2] = pa;
        const double (*b)[2] = pb;
        if ( (*a)[1] < (*b)[1] ) return -1;
        if ( (*a)[1] > (*b)[1] ) return +1;
        return 0;
        }


int compare_error_asc2( const void *pa, const void *pb )
        {
        const double (*a)[3] = pa;
        const double (*b)[3] = pb;
        if ( (*a)[2] < (*b)[2] ) return -1;
        if ( (*a)[2] > (*b)[2] ) return +1;
        return 0;
        }


static inline double sq(double diff)
        {
        return diff*diff;
        }


void VEProj (double atom1[3], double atom2[3], double newatom1[3], double newatom2[3], double dmin)
        {
        int i;
        double dist;

        dist = 0.;
        for(i=0;i<3;++i)
                dist += sq(atom2[i]-atom1[i]);

        dist = sqrt(dist);

        if (dist < dmin)
                {
                for(i=0;i<3;++i)
                        {
                        newatom1[i] = atom1[i]+((dmin-dist)/(2*dist))*(atom1[i]-atom2[i]);
                        newatom2[i] = atom2[i]-((dmin-dist)/(2*dist))*(atom1[i]-atom2[i]);
                        }

                }
        else
                {
                for(i=0;i<3;++i)
                        {
                        newatom1[i] = atom1[i];
                        newatom2[i] = atom2[i];
                        }
                }
        }


void changeVar(double (**BBo[2])[20][3], double (***SBo[2])[3], double (****VBo)[3], double (***VSo)[2][3], double (**RSo)[3], double (**VOo)[3])
        {
        int i,j,k,m,n,p,q,l,c;

        for (i=0;i<ProLen-1;++i)
                for (j=0;j<ProLen-1;++j)
                        for (m=0;m<5;++m)
                                for (n=0;n<3;++n)
					for(c=0;c<2;++c)
						{
						BBo[c][i][j][m][n] = atom[i][m][n];
						BBo[c][i][j][m+5][n] = atom[i+1][m][n];
						BBo[c][i][j][m+10][n] = atom[j][m][n];
						BBo[c][i][j][m+15][n] = atom[j+1][m][n];
						}

        for (i=0;i<ProLen-2;++i)
                for (j=0;j<NumSol;++j)
			for(c=0;c<2;++c)
				{
				for (m=0;m<5;++m)
					for (n=0;n<3;++n)
						{
						SBo[c][i][j][m][n] = atom[i][m][n];
						SBo[c][i][j][m+5][n] = atom[i+1][m][n];
						SBo[c][i][j][m+10][n] = atom[i+2][m][n];
						}

				for (n=0;n<3;++n)
					SBo[c][i][j][15][n] = Solatom[j][n];

				if (lenAA[RT[i+1]] > 5)
					for (m=5;m<lenAA[RT[i+1]];++m)
						for (n=0;n<3;++n)
							SBo[c][i][j][m+11][n] = atom[i+1][m][n];
				}

        for(i=0;i<ProLen;++i)
                for(j=i+1;j<ProLen;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				for(n=0;n<lenAA[RT[j]];++n)
					for (k=0;k<3;++k)
						{
						VBo[i][j][m][n][k] = atom[i][m][k];
						VBo[j][i][m][n][k] = atom[j][n][k];
						}

        for(i=0;i<ProLen;++i)
                for(j=0;j<NumSol;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				for (k=0;k<3;++k)
					{
					VSo[i][j][m][0][k] = atom[i][m][k];
					VSo[i][j][m][1][k] = Solatom[j][k];
					}

        for(i=0;i<NumSol;++i)
                for(j=0;j<NumSol;++j)
			for (k=0;k<3;++k)
				VOo[i][j][k] = Solatom[i][k];

        for (i=0;i<ProLen;++i)
		{
		if (i == 0)
			{
			for (k=0;k<3;++k)
				{
				RSo[i][0][k] = atom[i+1][1][k];
				RSo[i][1][k] = atom[i+1][0][k];
				}
			q = 2;
			}

		else if (i == ProLen-1)
			{
			for (k=0;k<3;++k)
				{
				RSo[i][0][k] = atom[i-1][1][k];
				RSo[i][1][k] = atom[i-1][2][k];
				RSo[i][2][k] = atom[i-1][3][k];
				}
			q = 3;
			}

		else 
			{
			for (k=0;k<3;++k)
				{
				RSo[i][0][k] = atom[i-1][1][k];
				RSo[i][1][k] = atom[i-1][2][k];
				RSo[i][2][k] = atom[i-1][3][k];
				RSo[i][3][k] = atom[i+1][1][k];
				RSo[i][4][k] = atom[i+1][0][k];
				}
			q = 5;
			}

		for (j=q;j<q+lenAA[RT[i]];++j)
			for (k=0;k<3;++k)
				RSo[i][j][k] = atom[i][j-q][k];

		}
        }


void chull(double (**points)[3], double (*pointsSol)[3])
        {
        int i,j,k,l,r,n,z,sign,NearPlane,NumFaces,FCChull[100][3], SPatom_count, atom_num;
        double vec1[3],vec2[3],perp[3],dp1,dp2, projVecSol[NumSol][3], NearDisSol[NumSol][2];
        double LagM, projVec[NumSPAA][3], NearDis[NumSPAA][6];

        NumFaces = 0;


        for(i=0;i<NumHPAA;++i)
                for(j=i+1;j<NumHPAA;++j)
                        for(k=j+1;k<NumHPAA;++k)
                                {
                                for(n=0;n<3;++n)
                                        {
                                        vec1[n] = points[HPAA[i]][4][n] - points[HPAA[j]][4][n];
                                        vec2[n] = points[HPAA[k]][4][n] - points[HPAA[j]][4][n];
                                        }

                                perp[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
                                perp[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
                                perp[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];

                                r = 0;
                                while (r==i||r==j|r==k) {r++;}

                                if (r>=NumHPAA)
                                        printf("very few hydrophobic residues");


                                dp1 = 0.;
                                for(n=0;n<3;++n)
                                        dp1 += (perp[n]*(points[HPAA[r]][4][n] - points[HPAA[j]][4][n]));

                                for(l=0;l<NumHPAA;++l)
                                        {
                                        if (l==i||l==j||l==k||l==r)
                                                continue;

                                        dp2 = 0.;
                                        for(n=0;n<3;++n)
                                                {
                                                dp2 += perp[n]*(points[HPAA[l]][4][n] - points[HPAA[j]][4][n]);
                                                }

                                        if (dp1*dp2 < 0.00000000)
                                                {
                                                sign = 0;
                                                break;
                                                }

                                        else
                                                sign = 1;

                                        }

                                if (sign == 1)
                                        {
                                        FCChull[NumFaces][0] = HPAA[i];
                                        FCChull[NumFaces][1] = HPAA[j];
                                        FCChull[NumFaces][2] = HPAA[k];
                                        ++NumFaces;
                                        }
                                }

	double disFaces[NumFaces];


        //Convex Hull Constraint for Strongly Hydrophilic AA 

	SPatom_count = 0;

        for(i=0;i<NumSPAA;++i)
		for(atom_num=4;atom_num<lenAA[RT[SPAA[i]]];++atom_num)
			{
			if (ATR[RT[SPAA[i]]][atom_num] != 1 && ATR[RT[SPAA[i]]][atom_num] != 2)			
				continue;


			NearDis[SPatom_count][0] = 1000.;
			NearDis[SPatom_count][1] = SPAA[i];
			NearDis[SPatom_count][2] = atom_num;

			for(j=0;j<NumFaces;++j)
				{
				k = 0;
				while ( (HPAA[k]==FCChull[j][0]) || (HPAA[k]==FCChull[j][1]) || (HPAA[k]==FCChull[j][2]))
					k++;

				for(n=0;n<3;++n)
					{
					vec1[n] = points[FCChull[j][0]][4][n] - points[FCChull[j][1]][4][n];
					vec2[n] = points[FCChull[j][2]][4][n] - points[FCChull[j][1]][4][n];
					}

				perp[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
				perp[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
				perp[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];

				dp1 = 0.;
				dp2 = 0.;

				for(n=0;n<3;++n)
					{
					dp1 += perp[n]*(points[HPAA[k]][4][n] - points[FCChull[j][1]][4][n]);
					dp2 += perp[n]*(points[SPAA[i]][atom_num][n] - points[FCChull[j][1]][4][n]);
					}

				if (dp1*dp2 < 0.0)
					{
					sign = 0;
					NearDis[SPatom_count][0] = 0.;

					NearDis[SPatom_count][3] = points[SPAA[i]][atom_num][0];
					NearDis[SPatom_count][4] = points[SPAA[i]][atom_num][1];
					NearDis[SPatom_count][5] = points[SPAA[i]][atom_num][2];
	
					break;
					}
				else
					{
					sign = 1;

					disFaces[j] = 0.;
					for(n=0;n<3;++n)
						disFaces[j] += perp[n]*(points[SPAA[i]][atom_num][n]-points[FCChull[j][1]][4][n]);

					LagM = -1.*disFaces[j]/(sq(perp[0])+sq(perp[1])+sq(perp[2]));

					disFaces[j] = fabs(disFaces[j])/sqrt(sq(perp[0])+sq(perp[1])+sq(perp[2]));

					if(disFaces[j] < NearDis[SPatom_count][0])
						{
						NearDis[SPatom_count][0] = disFaces[j];

						NearDis[SPatom_count][3] = points[SPAA[i]][atom_num][0] + perp[0]*LagM;
						NearDis[SPatom_count][4] = points[SPAA[i]][atom_num][1] + perp[1]*LagM;
						NearDis[SPatom_count][5] = points[SPAA[i]][atom_num][2] + perp[2]*LagM;
						}
					}
				}

			SPatom_count++;
			}

        qsort(NearDis, SPatom_count, 6*sizeof(double), compare_error_asc_G);

        z = (int) (fracSP*SPatom_count+0.5);

        for (i=0;i<z;++i)
                {
                j = (int) (NearDis[i][1]+0.5);
		k = (int) (NearDis[i][2]+0.5);

                for(n=0;n<3;++n)
                        points[j][k][n] = NearDis[i][3+n];
		}


        //Convex Hull Constraint for Solvents

        for(i=0;i<NumSol;++i)
                {
                NearDisSol[i][0] = i;
                NearDisSol[i][1] = 1000.;

                for(j=0;j<NumFaces;++j)
                        {
                        k = 0;
                        while ( (HPAA[k]==FCChull[j][0]) || (HPAA[k]==FCChull[j][1]) || (HPAA[k]==FCChull[j][2]))
                                k++;

                        for(n=0;n<3;++n)
                                {
                                vec1[n] = points[FCChull[j][0]][4][n] - points[FCChull[j][1]][4][n];
                                vec2[n] = points[FCChull[j][2]][4][n] - points[FCChull[j][1]][4][n];
                                }

                        perp[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1];
                        perp[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2];
                        perp[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0];

                        dp1 = 0.;
                        dp2 = 0.;

                        for(n=0;n<3;++n)
                                {
                                dp1 += perp[n]*(points[HPAA[k]][4][n] - points[FCChull[j][1]][4][n]);
                                dp2 += perp[n]*(pointsSol[i][n] - points[FCChull[j][1]][4][n]);
                                }

                        if (dp1*dp2 < 0.0000000001)
                                {
                                sign = 0;
                                NearDisSol[i][1] = 0.;

                                projVecSol[i][0] = pointsSol[i][0];
                                projVecSol[i][1] = pointsSol[i][1];
                                projVecSol[i][2] = pointsSol[i][2];

                                break;
                                }
                        else
                                {
                                sign = 1;

                                disFaces[j] = 0.;
                                for(n=0;n<3;++n)
                                        disFaces[j] += perp[n]*(pointsSol[i][n]-points[FCChull[j][1]][4][n]);

                                LagM = -1.*disFaces[j]/(sq(perp[0])+sq(perp[1])+sq(perp[2]));

                                disFaces[j] = fabs(disFaces[j])/sqrt(sq(perp[0])+sq(perp[1])+sq(perp[2]));

                                if(disFaces[j] < NearDisSol[i][1])
                                        {
                                        NearDisSol[i][1] = disFaces[j];

                                        projVecSol[i][0] = pointsSol[i][0] + perp[0]*LagM;
                                        projVecSol[i][1] = pointsSol[i][1] + perp[1]*LagM;
                                        projVecSol[i][2] = pointsSol[i][2] + perp[2]*LagM;
                                        }
                                }
                        }
                }

        qsort(NearDisSol, NumSol, 2*sizeof(double), compare_error_asc);

        z = (int) (fracOS*NumSol);

        for (i=0;i<z;++i)
                {
                j = (int) (NearDisSol[i][0]+0.5);
                for(n=0;n<3;++n)
                        pointsSol[j][n] = projVecSol[j][n];

                }
        } 


void inertia(double (**points)[3], double (*pointsSol)[3], double cut)
        {
        int i,j,k,m,n,min_axis,max_axis,totatoms;
        double IT[3][3];
        double COM[3],ax[3][3],tmp[3],eigenV[3],minEV = 10000000.,maxEV = 0.,lambda,dp;

        for(m=0;m<3;++m)
                {
                COM[m] = 0.;

                for(n=0;n<3;++n)
                        IT[m][n] = 0.;
                }

	totatoms = 0;
        for(i=0;i<ProLen;++i)
                for(j=0;j<lenAA[RT[i]];++j)
			{
			++totatoms;
                        for(n=0;n<3;++n)
                                COM[n] += points[i][j][n];
			}	

        for(n=0;n<3;++n)
                COM[n] /= totatoms;

        for(j=0;j<ProLen;++j)
                for(i=0;i<lenAA[RT[j]];++i)
                        {
                        IT[0][0] += sq(points[j][i][1] - COM[1]) + sq(points[j][i][2] - COM[2]);
                        IT[1][1] += sq(points[j][i][0] - COM[0]) + sq(points[j][i][2] - COM[2]);
                        IT[2][2] += sq(points[j][i][1] - COM[1]) + sq(points[j][i][0] - COM[0]);
                        IT[1][0] -= (points[j][i][0] - COM[0])*(points[j][i][1] - COM[1]);
                        IT[2][0] -= (points[j][i][0] - COM[0])*(points[j][i][2] - COM[2]);
                        IT[1][2] -= (points[j][i][2] - COM[2])*(points[j][i][1] - COM[1]);
                        }

        IT[2][1] = IT[1][2];
        IT[0][1] = IT[1][0];
        IT[0][2] = IT[2][0];

	totiter = 50;
        n = axes(IT,ax[0],ax[1]);
        normalize(ax[0]);
        normalize(ax[1]);
        cross(ax[0],ax[1],ax[2]);
        normalize(ax[2]);

        mult(IT,ax[0],tmp);
        eigenV[0] = tmp[0]*ax[0][0]+tmp[1]*ax[0][1]+tmp[2]*ax[0][2];

        mult(IT,ax[1],tmp);
        eigenV[1] = tmp[0]*ax[1][0]+tmp[1]*ax[1][1]+tmp[2]*ax[1][2];

        mult(IT,ax[2],tmp);
        eigenV[2] = tmp[0]*ax[2][0]+tmp[1]*ax[2][1]+tmp[2]*ax[2][2];

        for(m=0;m<3;++m)
                {
                if (eigenV[m]<minEV)
                        {
                        minEV = eigenV[m];
                        min_axis = m;
                        }

                if (eigenV[m]>maxEV)
                        {
                        maxEV = eigenV[m];
                        max_axis = m;
                        }
                }

        lambda = sqrt(maxEV/(totatoms))/cut;


        if (lambda > 1.)
                {
                for(j=0;j<ProLen;++j)
                        for(i=0;i<lenAA[RT[j]];++i)
                                {
                                dp = 0.;
                                for(n=0;n<3;++n)
                                        dp += (points[j][i][n]-COM[n])*ax[max_axis][n];

                                for(n=0;n<3;++n)
                                        {
                                        tmp[n] = COM[n]+dp*ax[max_axis][n];
                                        points[j][i][n] = ((points[j][i][n]-tmp[n])/lambda) + tmp[n];
                                        }
                                }

                for(i=0;i<NumSol;++i)
                        {
                        dp = 0.;
                        for(n=0;n<3;++n)
                                dp += (pointsSol[i][n]-COM[n])*ax[max_axis][n];

                        for(n=0;n<3;++n)
                                {
                                tmp[n] = COM[n]+dp*ax[max_axis][n];
                                pointsSol[i][n] = ((pointsSol[i][n]-tmp[n])/lambda) + tmp[n];
                                }
                        }

                }
        }


void projA(double (**BBo[2])[20][3], double (***SBo[2])[3], double (****VBo)[3], double (***VSo)[2][3], double (**RSo)[3], double (**VOo)[3])
        {
        int i,j,k,m,n,z,p,q,r,c,CountHBB,CountHBS,resdiff,totMotifs,max_motifs,num_atoms,motif_type;
        double tempBB[20][3], tempSol[25][3], tempRS[19][3];
        double tempweight, curweight;


        // Volume Exclusion Constraint

        for(i=0;i<ProLen;++i)
		{
                for(j=i+2;j<ProLen;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				for(n=0;n<lenAA[RT[j]];++n)
                                        VEProj( VBo[i][j][m][n], VBo[j][i][m][n], VBA[i][j][m][n], VBA[j][i][m][n], ATDMIN[ATR[RT[i]][m]][ATR[RT[j]][n]]);

		j = i+1;

		if (lenAA[RT[i]] > 5 && lenAA[RT[j]] > 5 && j!=ProLen)
			{
			for(m=5;m<lenAA[RT[i]];++m)
				for(n=5;n<lenAA[RT[j]];++n)
                                        VEProj( VBo[i][j][m][n], VBo[j][i][m][n], VBA[i][j][m][n], VBA[j][i][m][n], ATDMIN[ATR[RT[i]][m]][ATR[RT[j]][n]]);
			}
		}

        for(i=0;i<ProLen;++i)
                for(j=0;j<NumSol;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				VEProj( VSo[i][j][m][0], VSo[i][j][m][1], VSA[i][j][m][0], VSA[i][j][m][1], ASDMIN[ATR[RT[i]][m]]);

        for(i=0;i<NumSol;++i)
                for(j=i+1;j<NumSol;++j)
			VEProj( VOo[i][j], VOo[j][i], VOA[i][j], VOA[j][i], OODMIN);

        // Hydrogen Bonding Constraint

	max_motifs = 5;

        // step 1: get weights

        for (i=0;i<ProLen-1;++i)
                for (j=0;j<ProLen-1;++j)
                        {
                        for(m=0;m<20;++m)
                                for(n=0;n<3;++n)
					for(c=0;c<2;++c)
						BBA[c][i][j][m][n] = BBo[c][i][j][m][n];

                        WeightsBB[0][i][j] = 100000000000.;
                        WeightsBB[1][i][j] = 100000000000.;

                        resdiff = i-j;

                        if (abs(resdiff) == 0)
                                {
                                WrankB[i*(ProLen-1)+j][0] = i;
                                WrankB[i*(ProLen-1)+j][1] = j;
                                WrankB[i*(ProLen-1)+j][2] = WeightsBB[1][i][j];
                                continue;
                                }

                        else
                                {
                                if (abs(resdiff) > 5)
                                        resdiff = 0;

                                totMotifs = BBMotifNum[RT[i]][RT[i+1]][RT[j]][RT[j+1]][resdiff+5];

                                if (totMotifs == 0)
                                        {
                                        WrankB[i*(ProLen-1)+j][0] = i;
                                        WrankB[i*(ProLen-1)+j][1] = j;
                                        WrankB[i*(ProLen-1)+j][2] = WeightsBB[1][i][j];

                                        continue;
                                        }

				totMotifs = (totMotifs < max_motifs) ? totMotifs : max_motifs;

				for(c=0;c<2;++c)
					for(z=0;z<totMotifs;++z)
						{
						totiter = 10;
						motifproj(20,BBMotifs[RT[i]][RT[i+1]][RT[j]][RT[j+1]][resdiff+5][z],BBo[c][i][j],tempBB);
						tempweight = RMSD(20,BBo[c][i][j],tempBB);
						if (tempweight < WeightsBB[c][i][j])
							{
							WeightsBB[c][i][j] = tempweight;

							for(m=0;m<20;++m)
								for(n=0;n<3;++n)
									BBStore[c][i][j][m][n] = tempBB[m][n];
							}
						}
                                }

                        WrankB[i*(ProLen-1)+j][0] = i;
                        WrankB[i*(ProLen-1)+j][1] = j;
                        WrankB[i*(ProLen-1)+j][2] = etaBB[1][i][j]*WeightsBB[1][i][j];
                        }


        for (i=0;i<ProLen-2;++i)
                for (j=0;j<NumSol;++j)
                        {
			num_atoms = (lenAA[RT[i+1]] < 5) ? 16 : 11+lenAA[RT[i+1]];

                        for(m=0;m<num_atoms;++m)
                                for(n=0;n<3;++n)
					for(c=0;c<2;++c)
						SBA[c][i][j][m][n] = SBo[c][i][j][m][n];

                        WeightsSB[0][i][j] = 100000000000.;
                        WeightsSB[1][i][j] = 100000000000.;

                        totMotifs = SolMotifNum[RT[i]][RT[i+1]][RT[i+2]];

			for(c=0;c<2;++c)
				for(z=0;z<totMotifs;++z)
					{
					totiter = 10;
					motifproj(num_atoms,SolMotifs[RT[i]][RT[i+1]][RT[i+2]][z],SBo[c][i][j],tempSol);

					tempweight = RMSD(num_atoms,SBo[c][i][j],tempSol);
					if (tempweight < WeightsSB[c][i][j])
						{
						WeightsSB[c][i][j] = tempweight;

						for(m=0;m<num_atoms;++m)
							for(n=0;n<3;++n)
								SolStore[c][i][j][m][n] = tempSol[m][n];
						}
					}

                        WrankS[i*(NumSol)+j][0] = i;
                        WrankS[i*(NumSol)+j][1] = j;
                        WrankS[i*(NumSol)+j][2] = etaSB[1][i]*WeightsSB[1][i][j];
                        }

	for (i=0;i<ProLen-1;++i)
		{
	
		for (j=0;j<ProLen-1;++j)
			{
			Wrank[j][0] = j;	
			Wrank[j][1] = etaBB[0][i][j]*WeightsBB[0][i][j];
			}

		SizeWM = ProLen-1;

		if (i!=0)
			{
			for (j=0;j<NumSol;++j)
				{
				Wrank[SizeWM+j][0] = SizeWM+j;
				Wrank[SizeWM+j][1] = etaSB[0][i-1]*WeightsSB[0][i-1][j];
				}

			SizeWM += NumSol; 
			}

		if (i!=ProLen-2)
			{
			for (j=0;j<NumSol;++j)
				{
				Wrank[SizeWM+j][0] = SizeWM+j;
				Wrank[SizeWM+j][1] = etaSB[0][i]*WeightsSB[0][i][j];
				}

			SizeWM += NumSol; 
			}

                qsort(Wrank, SizeWM, 2*sizeof(double), compare_error_asc);

                MatchArray[i][0] = (int) (Wrank[0][0]+.5);
                MatchArray[i][1] = (int) (Wrank[1][0]+.5);

                for(z=0;z<2;++z)
                        {
                        if (z==1 && ((RT[i]==10 || RT[i+1]==10) || (RT[i]==11 || RT[i+1]==11)))
                                {
                                MatchArray[i][1] = -1;
                                continue;
                                }

                        if (MatchArray[i][z] < ProLen-1)
                                {
                                for(m=0;m<20;++m)
                                        for(n=0;n<3;++n)
                                                BBA[0][i][MatchArray[i][z]][m][n] = BBStore[0][i][MatchArray[i][z]][m][n];
                                }

                        else
                                {
				if (i==0)
					{
					p = i;
					q = MatchArray[i][z]-ProLen+1;
					}

				else if (i==ProLen-2)
					{
					p = i-1;
					q = MatchArray[i][z]-ProLen+1;
					}

				else
					{
					if (MatchArray[i][z] < ProLen-1+NumSol)
						{
						p = i-1;
						q = MatchArray[i][z]-ProLen+1;
						}
					
					else
						{
						p = i;
						q = MatchArray[i][z]-ProLen+1-NumSol;
						}
					}


				num_atoms = (lenAA[RT[p+1]] < 5) ? 16 : 11+lenAA[RT[p+1]];

                                for(m=0;m<num_atoms;++m)
                                        for(n=0;n<3;++n)
                                                SBA[0][p][q][m][n] = SolStore[0][p][q][m][n];
                                }
                        }
                }


        qsort(WrankB, (ProLen-1)*(ProLen-1), 3*sizeof(double), compare_error_asc2);

        CountHBB = 0;
        CountHBS = 0;

        r = 0;

        while (CountHBB < totHBB)
                {
                p = (int) (WrankB[r][0]+.5);
                q = (int) (WrankB[r][1]+.5);

		for(m=0;m<20;++m)
			for(n=0;n<3;++n)
				BBA[1][p][q][m][n] = BBStore[1][p][q][m][n];

		CountHBB++;

                r++;
                }
/*
	for(i=0;i<ProLen-1;++i)
		for(j=0;j<NumSol;++j)
			{
			if (SBMatch[j][i] == 1)
				{
				for(m=0;m<11;++m)
					for(n=0;n<3;++n)
						SBA[i][j][m][n] = SolStore[i][j][m][n];
				}
			else
				{
				for(m=0;m<11;++m)
					for(n=0;n<3;++n)
						SBA[i][j][m][n] = SBo[i][j][m][n];

				}
			}
*/

/*
        qsort(WrankS, (ProLen-2)*(NumSol), 3*sizeof(double), compare_error_asc2);

        r = 0;


        while (CountHBS < totHBS)
                {
                p = (int) (WrankS[r][0]+.5);
                q = (int) (WrankS[r][1]+.5);
		
		num_atoms = (lenAA[RT[p+1]] < 5) ? 16 : 11+lenAA[RT[p+1]];

		for(m=0;m<num_atoms;++m)
			for(n=0;n<3;++n)
				SBA[1][p][q][m][n] = SolStore[1][p][q][m][n];

		CountHBS++;

                r++;
                }

	int count_motifs_rot[21];
	for (i=0;i<21;++i)
		count_motifs_rot[i] = 0;

*/

	// Rotamer Constraint

	for(i=0;i<ProLen;++i)
		{
		if (i==0)
			{
			motif_type = 2;
			num_atoms = 2+lenAA[RT[i]];
			}

		else if (i==ProLen-1)
			{
			motif_type = 1;
			num_atoms = 3+lenAA[RT[i]];
			}

		else 
			{
			motif_type = 0;
			num_atoms = 5+lenAA[RT[i]];
			}

		curweight = 1000000.;

		for(z=0;z<RotMotifNum[motif_type][RT[i]];++z)
			{
			//z = count_motifs_rot[RT[i]];
			//++count_motifs_rot[RT[i]];

			totiter = 50;
			motifproj(num_atoms,RotMotifs[motif_type][RT[i]][z],RSo[i],tempRS);
			tempweight = RMSD(num_atoms,RSo[i],tempRS);
			
			if (tempweight < curweight)
				{
				curweight = tempweight;

				for(m=0;m<num_atoms;++m)
					for(n=0;n<3;++n)
						RSA[i][m][n] = tempRS[m][n];
				}

			}

		//if (iter == 0)
		//	printf("%d %d %d %d %lf %lf\n",i,motif_type,RT[i],RotMotifNum[motif_type][RT[i]], curweight, etaRS[i]);
		}
        }


void reflect(double (**BBo[2])[20][3], double (***SBo[2])[3], double (****VBo)[3], double (***VSo)[2][3], double (**RSo)[3], double (**VOo)[3])
        {
        int i,j,n,p,q,c,num_atoms;

        for(i=0;i<ProLen;++i)
                for(j=i+1;j<ProLen;++j)
			for(p=0;p<lenAA[RT[i]];++p)
				for(q=0;q<lenAA[RT[j]];++q)
					for(n=0;n<3;++n)
						{
						VBR[i][j][p][q][n] = 2.*VBo[i][j][p][q][n] - VB[i][j][p][q][n];
						VBR[j][i][p][q][n] = 2.*VBo[j][i][p][q][n] - VB[j][i][p][q][n];
						}

        for(i=0;i<ProLen;++i)
                for(j=0;j<NumSol;++j)
			for(p=0;p<lenAA[RT[i]];++p)
				for(n=0;n<3;++n)
					{
					VSR[i][j][p][0][n] = 2.*VSo[i][j][p][0][n] - VS[i][j][p][0][n];
					VSR[i][j][p][1][n] = 2.*VSo[i][j][p][1][n] - VS[i][j][p][1][n];
					}

	for(i=0;i<NumSol;++i)
		for(j=i+1;j<NumSol;++j)
			for(n=0;n<3;++n)
				{
				VOR[i][j][n] = 2.*VOo[i][j][n] - VO[i][j][n];
				VOR[j][i][n] = 2.*VOo[j][i][n] - VO[j][i][n];
				}

        for (i=0;i<ProLen-1;++i)
                for (j=0;j<ProLen-1;++j)
                        for(p=0;p<20;++p)
                                for(n=0;n<3;++n)
					for(c=0;c<2;++c)
						BBR[c][i][j][p][n] = 2.*BBo[c][i][j][p][n] - BB[c][i][j][p][n];

        for (i=0;i<ProLen-2;++i)
                for (j=0;j<NumSol;++j)
			{
			num_atoms = (lenAA[RT[i+1]] < 5) ? 16 : 11+lenAA[RT[i+1]];
                        for(p=0;p<num_atoms;++p)
                                for(n=0;n<3;++n)
					for(c=0;c<2;++c)
						SBR[c][i][j][p][n] = 2.*SBo[c][i][j][p][n] - SB[c][i][j][p][n];
			}

	for(i=0;i<ProLen;++i)
		{
		if (i==0)
			num_atoms = 2+lenAA[RT[i]];

		else if (i==ProLen-1)
			num_atoms = 3+lenAA[RT[i]];

		else 
			num_atoms = 5+lenAA[RT[i]];

		for(j=0;j<num_atoms;++j)
			for(n=0;n<3;++n)
				RSR[i][j][n] = 2.*RSo[i][j][n] - RS[i][j][n];
		}
        }


void projB(double (**BBo[2])[20][3], double (***SBo[2])[3], double (****VBo)[3], double (***VSo)[2][3], double (**RSo)[3], double (**VOo)[3])
        {
        int i,j,k,m,n,p,q,l,c,num_atoms;
        double arbd,NumVar[ProLen][14],NumVarSol[NumSol];

        // Step 1: get atom positions from BB,VS,SB,VS,RS

        for (i=0;i<ProLen;++i)
		{
		num_atoms = (lenAA[RT[i]] == 4) ? 5 : lenAA[RT[i]]; 
                for (j=0;j<num_atoms;++j)
			{
			NumVar[i][j] = 0.;

                        for (n=0;n<3;++n)
                                atom[i][j][n] = 0.;
			}
		}

        for (j=0;j<NumSol;++j)
		{
		NumVarSol[j] = 0.;

                for (n=0;n<3;++n)
                        Solatom[j][n] = 0.;
		}

        for (i=0;i<ProLen-1;++i)
                for (j=0;j<ProLen-1;++j)
			for(c=0;c<2;++c)
				{
				if (i==j)
					continue;

				for (k=0;k<5;++k)
					{
					NumVar[i][k] += sq(etaBB[c][i][j]);
					NumVar[i+1][k] += sq(etaBB[c][i][j]);
					NumVar[j][k] += sq(etaBB[c][i][j]);
					NumVar[j+1][k] += sq(etaBB[c][i][j]);

					for (n=0;n<3;++n)
						{
						atom[i][k][n] += sq(etaBB[c][i][j])*BBo[c][i][j][k][n];
						atom[i+1][k][n] += sq(etaBB[c][i][j])*BBo[c][i][j][k+5][n];
						atom[j][k][n] += sq(etaBB[c][i][j])*BBo[c][i][j][k+10][n];
						atom[j+1][k][n] += sq(etaBB[c][i][j])*BBo[c][i][j][k+15][n];
						}
					}
				}

        for (i=0;i<ProLen-2;++i)
                for (j=0;j<NumSol;++j)
			for(c=0;c<1;++c)
				{
				for (k=0;k<5;++k)
					{
					NumVar[i][k] += sq(etaSB[c][i]);
					NumVar[i+1][k] += sq(etaSB[c][i]);
					NumVar[i+2][k] += sq(etaSB[c][i]);

					for (n=0;n<3;++n)
						{
						atom[i][k][n] += sq(etaSB[c][i])*SBo[c][i][j][k][n];
						atom[i+1][k][n] += sq(etaSB[c][i])*SBo[c][i][j][k+5][n];
						atom[i+2][k][n] += sq(etaSB[c][i])*SBo[c][i][j][k+10][n];
						}
					}				


				NumVarSol[j] += sq(etaSB[c][i]);

				for (n=0;n<3;++n)
					Solatom[j][n] += sq(etaSB[c][i])*SBo[c][i][j][15][n];

				if (lenAA[RT[i+1]] > 5)
					for (k=5;k<lenAA[RT[i+1]];++k)
						{
						NumVar[i+1][k] += sq(etaSB[c][i]);
						for (n=0;n<3;++n)
							atom[i+1][k][n] += sq(etaSB[c][i])*SBo[c][i][j][k+11][n];
						}
				}

        for(i=0;i<ProLen;++i)
		{
                for(j=i+2;j<ProLen;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				for(n=0;n<lenAA[RT[j]];++n)
					{
					NumVar[i][m] += sq(etaVB[i][j][m][n]);
					NumVar[j][n] += sq(etaVB[i][j][m][n]);

					for (k=0;k<3;++k)
						{
						atom[i][m][k] += sq(etaVB[i][j][m][n])*VBo[i][j][m][n][k];
						atom[j][n][k] += sq(etaVB[i][j][m][n])*VBo[j][i][m][n][k];
						}
					}

		j = i+1;

		if (lenAA[RT[i]] > 5 && lenAA[RT[j]] > 5 && j != ProLen)
			{
			for(m=5;m<lenAA[RT[i]];++m)
				for(n=5;n<lenAA[RT[j]];++n)
					{
					NumVar[i][m] += sq(etaVB[i][j][m][n]);
					NumVar[j][n] += sq(etaVB[i][j][m][n]);

					for (k=0;k<3;++k)
						{
						atom[i][m][k] += sq(etaVB[i][j][m][n])*VBo[i][j][m][n][k];
						atom[j][n][k] += sq(etaVB[i][j][m][n])*VBo[j][i][m][n][k];
						}
					}
			}
		}

        for(i=0;i<ProLen;++i)
                for(j=0;j<NumSol;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				{
				NumVar[i][m] += sq(etaVS[i][j][m]);
				NumVarSol[j] += sq(etaVS[i][j][m]);

				for (k=0;k<3;++k)
					{
					atom[i][m][k] += sq(etaVS[i][j][m])*VSo[i][j][m][0][k];
					Solatom[j][k] += sq(etaVS[i][j][m])*VSo[i][j][m][1][k];
					}
				}

        for(i=0;i<NumSol;++i)
                for(j=i+1;j<NumSol;++j)
			{
			NumVarSol[i] += sq(etaVO[i][j]);
			NumVarSol[j] += sq(etaVO[i][j]);

			for (k=0;k<3;++k)
				{
				Solatom[i][k] += sq(etaVO[i][j])*VOo[i][j][k];
				Solatom[j][k] += sq(etaVO[i][j])*VOo[j][i][k];
				}
			}


        for (i=0;i<ProLen;++i)
		{
		if (i == 0)
			{
			NumVar[i+1][1] += sq(etaRS[i]);
			NumVar[i+1][0] += sq(etaRS[i]);

			for (k=0;k<3;++k)
				{
				atom[i+1][1][k] += sq(etaRS[i])*RSo[i][0][k];
				atom[i+1][0][k] += sq(etaRS[i])*RSo[i][1][k];
				}
			q = 2;
			}

		else if (i == ProLen-1)
			{
			NumVar[i-1][1] += sq(etaRS[i]);
			NumVar[i-1][2] += sq(etaRS[i]);
			NumVar[i-1][3] += sq(etaRS[i]);

			for (k=0;k<3;++k)
				{
				atom[i-1][1][k] += sq(etaRS[i])*RSo[i][0][k];
				atom[i-1][2][k] += sq(etaRS[i])*RSo[i][1][k];
				atom[i-1][3][k] += sq(etaRS[i])*RSo[i][2][k];
				}
			q = 3;
			}

		else 
			{
			NumVar[i-1][1] += sq(etaRS[i]);
			NumVar[i-1][2] += sq(etaRS[i]);
			NumVar[i-1][3] += sq(etaRS[i]);
			NumVar[i+1][1] += sq(etaRS[i]);
			NumVar[i+1][0] += sq(etaRS[i]);

			for (k=0;k<3;++k)
				{
				atom[i-1][1][k] += sq(etaRS[i])*RSo[i][0][k];
				atom[i-1][2][k] += sq(etaRS[i])*RSo[i][1][k];
				atom[i-1][3][k] += sq(etaRS[i])*RSo[i][2][k];
				atom[i+1][1][k] += sq(etaRS[i])*RSo[i][3][k];
				atom[i+1][0][k] += sq(etaRS[i])*RSo[i][4][k];
				}
			q = 5;
			}

		for (j=q;j<q+lenAA[RT[i]];++j)
			{
			NumVar[i][j-q] += sq(etaRS[i]);

			for (k=0;k<3;++k)
				atom[i][j-q][k] += sq(etaRS[i])*RSo[i][j][k];
			}
		}

			
        for (i=0;i<ProLen;++i)
		{
		num_atoms = (lenAA[RT[i]] == 4) ? 5 : lenAA[RT[i]]; 
                for (j=0;j<num_atoms;++j)
                        for (n=0;n<3;++n)
                                atom[i][j][n] /= NumVar[i][j];
		}

        for (j=0;j<NumSol;++j)
                for (n=0;n<3;++n)
                        Solatom[j][n] /= NumVarSol[j];

        // Step 2: Convex Hull Constraint

        inertia(atom,Solatom,9.0);

        if(iter%CHstep==0)
                chull(atom, Solatom);

        // Step 3: get BB, VS, VB, RS, SB from atom positions

        changeVar(BBB, SBB, VBB, VSB, RSB, VOB);
        }


void getRMSCa()
        {
        int i,j,k,m,n,p,q,c;
        int NumVar[ProLen];
        FILE *fp;

        for (i=0;i<ProLen;++i)
                {
                CAerr[i] = 0.;
                NumVar[i] = 0;
                }


        for (i=0;i<ProLen-1;++i)
                for (j=0;j<ProLen-1;++j)
			if (i!=j)
				for (n=0;n<3;++n)
					for(c=0;c<2;++c)
						{
						CAerr[i] += sq(BBB[c][i][j][1][n] - BBA[c][i][j][1][n]);
						++NumVar[i];
						CAerr[i+1] += sq(BBB[c][i][j][6][n] - BBA[c][i][j][6][n]) ;
						++NumVar[i+1];
						CAerr[j] += sq(BBB[c][i][j][11][n] - BBA[c][i][j][11][n]);
						++NumVar[j];
						CAerr[j+1] += sq(BBB[c][i][j][16][n] - BBA[c][i][j][16][n]) ;
						++NumVar[j+1];
						}


        for (i=0;i<ProLen-2;++i)
                for (j=0;j<NumSol;++j)
                        for (n=0;n<3;++n)
				for(c=0;c<1;++c)
					{
					CAerr[i] += sq(SBB[c][i][j][1][n] - SBA[c][i][j][1][n]);
					++NumVar[i];
					CAerr[i+1] += sq(SBB[c][i][j][6][n] - SBA[c][i][j][6][n]) ;
					++NumVar[i+1];
					CAerr[i+2] += sq(SBB[c][i][j][11][n] - SBA[c][i][j][11][n]) ;
					++NumVar[i+2];
					}

        for (i=0;i<ProLen;++i)
                for (j=i+2;j<ProLen;++j)
                        {
                        for (k=0;k<3;++k)
                                {
				for(n=0;n<lenAA[RT[j]];++n)
					{
					CAerr[i] += sq(VBB[i][j][1][n][k] - VBA[i][j][1][n][k]);
					++NumVar[i];
					}

				for(m=0;m<lenAA[RT[i]];++m)
					{
					CAerr[j] += sq(VBB[j][i][m][1][k] - VBA[j][i][m][1][k]);
					++NumVar[j];
                                        }
                                }
                        }

        for (i=0;i<ProLen;++i)
                for (j=0;j<NumSol;++j)
                        for (k=0;k<3;++k)
				{
				CAerr[i] += sq(VSB[i][j][1][0][k] - VSA[i][j][1][0][k]);
				++NumVar[i];
				}

        for (i=0;i<ProLen;++i)
		for (k=0;k<3;++k)
			{
			if (i==0)
				{
				CAerr[i+1] += sq(RSB[i][0][k] - RSA[i][0][k]);
				++NumVar[i+1];
				CAerr[i] += sq(RSB[i][2][k] - RSA[i][2][k]);
				++NumVar[i];
				}

			else if (i==ProLen-1)
				{
				CAerr[i-1] += sq(RSB[i][0][k] - RSA[i][0][k]);
				++NumVar[i-1];
				CAerr[i] += sq(RSB[i][3][k] - RSA[i][3][k]);
				++NumVar[i];
				}
		
			else
				{
				CAerr[i-1] += sq(RSB[i][0][k] - RSA[i][0][k]);
				++NumVar[i-1];
				CAerr[i+1] += sq(RSB[i][3][k] - RSA[i][3][k]);
				++NumVar[i+1];
				CAerr[i] += sq(RSB[i][6][k] - RSA[i][6][k]);
				++NumVar[i];
				}
			}

        fp=fopen(CAfile,"a");

        for(j=0;j<ProLen;++j)
                {
                CAerr[j] /= NumVar[j];
                CAerr[j] = sqrt(CAerr[j]);
                fprintf(fp,"%lf ",CAerr[j]);
                }

        fprintf(fp,"\n");

        fclose(fp);
        }


void RRR()
        {
        int i,j,k,m,n,p,q,c, NumVXVar,NumSBVar,NumRSVar,totCons, VXCons,num_atoms;
        double diff,avgeta,temperr;

	projA(BB, SB, VB, VS, RS, VO);
	reflect(BBA, SBA, VBA, VSA, RSA, VOA);
	projB(BBR, SBR, VBR, VSR, RSR, VOR);

        tBBerr = 0.;
        tVXerr = 0.;
        tSBerr = 0.;
        tRSerr = 0.;
        toterr = 0.;

        avgeta = 0.;

        for (i=0;i<ProLen-1;++i)
                {
                for (j=0;j<ProLen-1;++j)
			for(c=0;c<2;++c)
				{
				if (i==j)
					continue;

				temperr = 0.;
				for(p=0;p<20;++p)
					for(n=0;n<3;++n)
						{
						diff = BBB[c][i][j][p][n] - BBA[c][i][j][p][n];
						BB[c][i][j][p][n] += beta*diff;
						temperr += sq(diff);
						}
				tBBerr += temperr;
				temperr /= 60.;
				etaBB[c][i][j] += epsilon*temperr;
				avgeta += etaBB[c][i][j];
				}
                }

        toterr += tBBerr;
        tBBerr /= ((ProLen-1)*(ProLen-2)*2*60.);
        tBBerr = sqrt(tBBerr);

	NumSBVar = 0;

	for(c=0;c<1;++c)
		for (i=0;i<ProLen-2;++i)
			{
			temperr = 0.;
			num_atoms = (lenAA[RT[i+1]] < 5) ? 16 : 11+lenAA[RT[i+1]];
			NumSBVar += 3*num_atoms*NumSol;

			for (j=0;j<NumSol;++j)
				for(p=0;p<num_atoms;++p)
					for(n=0;n<3;++n)
						{
						diff = SBB[c][i][j][p][n] - SBA[c][i][j][p][n];
						SB[c][i][j][p][n] += beta*diff;
						temperr += sq(diff);
						}

			tSBerr += temperr;
			temperr /= (num_atoms*3.*NumSol);
			etaSB[c][i] += epsilon*temperr;
			avgeta += etaSB[c][i];
			}

        toterr += tSBerr;
        tSBerr /= NumSBVar;
        tSBerr = sqrt(tSBerr);

        NumVXVar = 0;

        for(i=0;i<ProLen;++i)
		{
                for(j=i+2;j<ProLen;++j)
			{
			NumVXVar += lenAA[RT[i]]*lenAA[RT[j]];
			for(m=0;m<lenAA[RT[i]];++m)
				for(n=0;n<lenAA[RT[j]];++n)
					{
					temperr = 0.;
					for (k=0;k<3;++k)
						{
						diff = VBB[i][j][m][n][k] - VBA[i][j][m][n][k];
						VB[i][j][m][n][k] += beta*diff;
						temperr += sq(diff);

						diff = VBB[j][i][m][n][k] - VBA[j][i][m][n][k];
						VB[j][i][m][n][k] += beta*diff;
						temperr += sq(diff);
						}
					tVXerr += temperr;
					temperr /= 6.;
					etaVB[i][j][m][n] += epsilon*temperr;
					avgeta += etaVB[i][j][m][n];
					}
			}

		j = i+1;

		if (lenAA[RT[i]] > 5 && lenAA[RT[j]] > 5 && j != ProLen)
			{
			NumVXVar += (lenAA[RT[i]]-5)*(lenAA[RT[j]]-5);
			for(m=5;m<lenAA[RT[i]];++m)
				for(n=5;n<lenAA[RT[j]];++n)
					{
					temperr = 0.;
					for (k=0;k<3;++k)
						{
						diff = VBB[i][j][m][n][k] - VBA[i][j][m][n][k];
						VB[i][j][m][n][k] += beta*diff;
						temperr += sq(diff);

						diff = VBB[j][i][m][n][k] - VBA[j][i][m][n][k];
						VB[j][i][m][n][k] += beta*diff;
						temperr += sq(diff);
						}
					tVXerr += temperr;
					temperr /= 6.;
					etaVB[i][j][m][n] += epsilon*temperr;
					avgeta += etaVB[i][j][m][n];
					}
			}
		}

        for(i=0;i<ProLen;++i)
                for(j=0;j<NumSol;++j)
			{
			NumVXVar += lenAA[RT[i]];
			for(m=0;m<lenAA[RT[i]];++m)
				{
				temperr = 0.;
				for (k=0;k<3;++k)
					{
					diff = VSB[i][j][m][0][k] - VSA[i][j][m][0][k];
					VS[i][j][m][0][k] += beta*diff;
					temperr += sq(diff);

					diff = VSB[i][j][m][1][k] - VSA[i][j][m][1][k];
					VS[i][j][m][1][k] += beta*diff;
					temperr += sq(diff);
					}
				tVXerr += temperr;
				temperr /= 6.;
				etaVS[i][j][m] += epsilon*temperr;
				avgeta += etaVS[i][j][m];
				}
			}

        for(i=0;i<NumSol;++i)
                for(j=i+1;j<NumSol;++j)
			{
			temperr = 0.;
			NumVXVar += 1;
			for (k=0;k<3;++k)
				{
				diff = VOB[i][j][k] - VOA[i][j][k];
				VO[i][j][k] += beta*diff;
				temperr += sq(diff);

				diff = VOB[j][i][k] - VOA[j][i][k];
				VO[j][i][k] += beta*diff;
				temperr += sq(diff);
				}
			tVXerr += temperr;
			temperr /= 6.;
                        etaVO[i][j] += epsilon*temperr;
                        avgeta += etaVO[i][j];
			}

        toterr += tVXerr;
        tVXerr /= (NumVXVar*6.);
        tVXerr = sqrt(tVXerr);

	NumRSVar = 0;

        for (i=0;i<ProLen;++i)
		{
		if (i==0)
			num_atoms = 2+lenAA[RT[i]];

		else if (i==ProLen-1)
			num_atoms = 3+lenAA[RT[i]];

		else 
			num_atoms = 5+lenAA[RT[i]];

		temperr = 0.; 
		NumRSVar += num_atoms;
		for (j=0;j<num_atoms;++j)
			for (k=0;k<3;++k)
				{
				diff = RSB[i][j][k] - RSA[i][j][k];
				RS[i][j][k] += beta*diff;
				temperr += sq(diff);
				}
		tRSerr += temperr;
		temperr /= (num_atoms*3.);
		etaRS[i] += epsilon*temperr;
		avgeta += etaRS[i];

		}

        toterr += tRSerr;
        tRSerr /= (NumRSVar*3.);
        tRSerr = sqrt(tRSerr);

        toterr /= (((ProLen-1)*(ProLen-2)*2*60.)+(NumVXVar*6.)+(NumRSVar*3.)+NumSBVar);

        toterr = sqrt(toterr);

        totCons = (((ProLen-1)*(ProLen-2)*2.)+NumVXVar+(2*ProLen-2));
        avgeta /= totCons;


        for (i=0;i<ProLen-1;++i)
                for (j=0;j<ProLen-1;++j)
			for(c=0;c<2;++c)
				etaBB[c][i][j] /= avgeta;

        for (i=0;i<ProLen-2;++i)
		for(c=0;c<1;++c)
			etaSB[c][i] /= avgeta;

        for (i=0;i<ProLen;++i)
                for (j=i+1;j<ProLen;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				for(n=0;n<lenAA[RT[j]];++n)
					etaVB[i][j][m][n] /= avgeta;

        for (i=0;i<ProLen;++i)
                for (j=0;j<NumSol;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				etaVS[i][j][m] /= avgeta;

        for (i=0;i<NumSol;++i)
                for (j=0;j<NumSol;++j)
                        etaVO[i][j] /= avgeta;

        for (i=0;i<ProLen;++i)
		//etaRS[i] = 1.;
		etaRS[i] /= avgeta;
/*
	for(c=0;c<2;++c)
		for (i=0;i<ProLen-2;++i)
			{
			avgeta = 0.;
			for (j=0;j<NumSol;++j)
				avgeta += etaSB[c][i][j];
	
			avgeta /= NumSol; 

			for (j=0;j<NumSol;++j)
				etaSB[c][i][j] = avgeta;
			}
*/
        }


void printsol(char *solfile)
        {    
        FILE *fp; 
        int i,j,k;

        fp=fopen(solfile,"a");

        for(i=0;i<5;++i)
                for(j=0;j<ProLen;++j)
                        {    
                        for(k=0;k<3;++k)
                                fprintf(fp,"%lf ",atom[j][i][k]);
                        fprintf(fp,"\n");
                        }    

        for(i=0;i<NumSol;++i)
                {    
                for(k=0;k<3;++k)
                        fprintf(fp,"%lf ",Solatom[i][k]);
                fprintf(fp,"\n");
                }    

        fprintf(fp,"\n");

        fclose(fp);
        }


int solve(int maxiter,int iterstride,double stoperr)
        {
        FILE *fperr,*fpsol;
        int i,j;

        fperr=fopen(errfile,"w");

        for(iter=0;iter<=maxiter;++iter)
                {    
                RRR();

                if(iter%iterstride==0)
                        {    
                        fprintf(fperr,"%.6e\t%.6e %.6e %.6e %.6e\n",toterr,tBBerr,tSBerr,tVXerr,tRSerr);
                        printsol(solfile);
                        getRMSCa();
                        }    

                if(toterr<stoperr)
                        {    
                        fprintf(fperr,"%.6e\t%.6e %.6e %.6e %.6e\n",toterr,tBBerr,tSBerr,tVXerr,tRSerr);
                        printsol(solfile);
                        getRMSCa();
                        fclose(fperr);
                        return iter;
                        }
		}

        fclose(fperr);
/*
	if (iter%10000 == 0)
		{
		changeVar(BB, SB, VB, VS, RS, VO);

		for (i=0;i<ProLen-1;++i)
			for (j=0;j<ProLen-1;++j)
				etaBB[i][j] = 1.;

		for (i=0;i<ProLen-1;++i)
			for (j=0;j<NumSol;++j)
				etaSB[i][j] = 1.;

		for (i=0;i<ProLen;++i)
			for (j=i+1;j<ProLen;++j)
				etaVB[i][j] = 1.;

		for (i=0;i<ProLen;++i)
			for (j=0;j<NumSol;++j)
				etaVS[i][j] = 1.;

		for (i=0;i<NumSol;++i)
			for (j=0;j<NumSol;++j)
				etaVO[i][j] = 1.;

		for (i=0;i<ProLen;++i)
			etaRS[i] = 1.;
		}
*/

        return 0;
        }


void makevars()
        {
        int i,j,k,m,n,DorA,c,num_atoms;

	for(c=0;c<2;++c)
		{
		BB[c] = malloc((ProLen-1)*sizeof(double*));
		BBA[c] = malloc((ProLen-1)*sizeof(double*));
		BBR[c] = malloc((ProLen-1)*sizeof(double*));
		BBB[c] = malloc((ProLen-1)*sizeof(double*));

		for(i=0;i<(ProLen-1);i++)
			{
			BB[c][i] = malloc((ProLen-1)*sizeof(double*[20][3]));
			BBA[c][i] = malloc((ProLen-1)*sizeof(double*[20][3]));
			BBR[c][i] = malloc((ProLen-1)*sizeof(double*[20][3]));
			BBB[c][i] = malloc((ProLen-1)*sizeof(double*[20][3]));
			}
		}

	for(c=0;c<2;++c)
		{
		SB[c] = malloc((ProLen-2)*sizeof(double**));
		SBA[c] = malloc((ProLen-2)*sizeof(double**));
		SBR[c] = malloc((ProLen-2)*sizeof(double**));
		SBB[c] = malloc((ProLen-2)*sizeof(double**));

		for(i=0;i<(ProLen-2);i++)
			{
			SB[c][i] = malloc(NumSol*sizeof(double*));
			SBA[c][i] = malloc(NumSol*sizeof(double*));
			SBR[c][i] = malloc(NumSol*sizeof(double*));
			SBB[c][i] = malloc(NumSol*sizeof(double*));
	
			num_atoms = (lenAA[RT[i+1]] < 5) ? 16 : 11+lenAA[RT[i+1]];

			for(j=0;j<NumSol;++j)
				{
				SB[c][i][j] = malloc(num_atoms*sizeof(double*[3]));
				SBA[c][i][j] = malloc(num_atoms*sizeof(double*[3]));
				SBR[c][i][j] = malloc(num_atoms*sizeof(double*[3]));
				SBB[c][i][j] = malloc(num_atoms*sizeof(double*[3]));
				}
			}
		}

        RS = malloc(ProLen*sizeof(double*));
        RSA = malloc(ProLen*sizeof(double*));
        RSR = malloc(ProLen*sizeof(double*));
        RSB = malloc(ProLen*sizeof(double*));

        for(i=0;i<ProLen;i++)
                {
		if (i==0)
			num_atoms = 2+lenAA[RT[i]];

		else if (i==ProLen-1)
			num_atoms = 3+lenAA[RT[i]];

		else 
			num_atoms = 5+lenAA[RT[i]];

		RS[i] = malloc(num_atoms*sizeof(double*[3]));
		RSA[i] = malloc(num_atoms*sizeof(double*[3]));
		RSR[i] = malloc(num_atoms*sizeof(double*[3]));
		RSB[i] = malloc(num_atoms*sizeof(double*[3]));
                }

        VB = malloc(ProLen*sizeof(double***));
        VBA = malloc(ProLen*sizeof(double***));
        VBR = malloc(ProLen*sizeof(double***));
        VBB = malloc(ProLen*sizeof(double***));

        for(i=0; i<ProLen; ++i)
                {
                VB[i] = malloc(ProLen*sizeof(double**));
                VBA[i] = malloc(ProLen*sizeof(double**));
                VBR[i] = malloc(ProLen*sizeof(double**));
                VBB[i] = malloc(ProLen*sizeof(double**));
		}
		
        for(i=0; i<ProLen; ++i)
		for(j=i+1; j<ProLen; ++j)
			{
			VB[i][j] = malloc(lenAA[RT[i]]*sizeof(double*));
			VBA[i][j] = malloc(lenAA[RT[i]]*sizeof(double*));
			VBR[i][j] = malloc(lenAA[RT[i]]*sizeof(double*));
			VBB[i][j] = malloc(lenAA[RT[i]]*sizeof(double*));

			VB[j][i] = malloc(lenAA[RT[i]]*sizeof(double*));
			VBA[j][i] = malloc(lenAA[RT[i]]*sizeof(double*));
			VBR[j][i] = malloc(lenAA[RT[i]]*sizeof(double*));
			VBB[j][i] = malloc(lenAA[RT[i]]*sizeof(double*));

			for(k=0; k<lenAA[RT[i]]; ++k)
				{
				VB[i][j][k] = malloc(lenAA[RT[j]]*sizeof(double*[3]));
				VBA[i][j][k] = malloc(lenAA[RT[j]]*sizeof(double*[3]));
				VBR[i][j][k] = malloc(lenAA[RT[j]]*sizeof(double*[3]));
				VBB[i][j][k] = malloc(lenAA[RT[j]]*sizeof(double*[3]));

				VB[j][i][k] = malloc(lenAA[RT[j]]*sizeof(double*[3]));
				VBA[j][i][k] = malloc(lenAA[RT[j]]*sizeof(double*[3]));
				VBR[j][i][k] = malloc(lenAA[RT[j]]*sizeof(double*[3]));
				VBB[j][i][k] = malloc(lenAA[RT[j]]*sizeof(double*[3]));
				}
			}

	VS = malloc(ProLen*sizeof(double**));
	VSA = malloc(ProLen*sizeof(double**));
	VSR = malloc(ProLen*sizeof(double**));
	VSB = malloc(ProLen*sizeof(double**));

	for(j=0; j<ProLen; ++j)
		{
		VS[j] = malloc(NumSol*sizeof(double*));
		VSA[j] = malloc(NumSol*sizeof(double*));
		VSR[j] = malloc(NumSol*sizeof(double*));
		VSB[j] = malloc(NumSol*sizeof(double*));

		for(k=0; k<NumSol; ++k)
			{
			VS[j][k] = malloc(lenAA[RT[j]]*sizeof(double*[2][3]));
			VSA[j][k] = malloc(lenAA[RT[j]]*sizeof(double*[2][3]));
			VSR[j][k] = malloc(lenAA[RT[j]]*sizeof(double*[2][3]));
			VSB[j][k] = malloc(lenAA[RT[j]]*sizeof(double*[2][3]));
			}
		}

	VO = malloc(NumSol*sizeof(double*));
	VOA = malloc(NumSol*sizeof(double*));
	VOR = malloc(NumSol*sizeof(double*));
	VOB = malloc(NumSol*sizeof(double*));

	for(i=0;i<NumSol;++i)
		{
		VO[i] = malloc(NumSol*sizeof(double*[3]));
		VOA[i] = malloc(NumSol*sizeof(double*[3]));
		VOR[i] = malloc(NumSol*sizeof(double*[3]));
		VOB[i] = malloc(NumSol*sizeof(double*[3]));
		}

	atom = malloc(ProLen*sizeof(double*));
        for (i=0;i<ProLen;++i)
		{
		if (lenAA[RT[i]] > 4)
			atom[i] = malloc(lenAA[RT[i]]*sizeof(double*[3]));
		else
			atom[i] = malloc((lenAA[RT[i]]+1)*sizeof(double*[3]));
		}	

        Solatom = malloc(NumSol*sizeof(double*[3]));

        CAerr = malloc(ProLen*sizeof(double));

        etaBB[0] = malloc((ProLen-1)*sizeof(double*));
        etaBB[1] = malloc((ProLen-1)*sizeof(double*));

        etaSB[0] = malloc((ProLen-2)*sizeof(double));
        etaSB[1] = malloc((ProLen-2)*sizeof(double));

	for(c=0;c<2;++c)
		for(i=0;i<(ProLen-1);i++)
			etaBB[c][i] = malloc((ProLen-1)*sizeof(double));

        etaVB = malloc(ProLen*sizeof(double***));
        etaVS = malloc(ProLen*sizeof(double**));

        for (i=0;i<ProLen;++i)
                {
                etaVB[i] = malloc(ProLen*sizeof(double**));
                etaVS[i] = malloc(NumSol*sizeof(double*));

		for(j=i+1; j<ProLen; ++j)
			{
			etaVB[i][j] =  malloc(lenAA[RT[i]]*sizeof(double*));
			for(m=0;m<lenAA[RT[i]];++m)
				etaVB[i][j][m] =  malloc(lenAA[RT[j]]*sizeof(double));
			}

		for(j=0; j<NumSol; ++j)
			etaVS[i][j] =  malloc(lenAA[RT[i]]*sizeof(double));
                }

	etaVO = malloc(NumSol*sizeof(double*));

	for(i=0;i<NumSol;++i)
		etaVO[i] = malloc(NumSol*sizeof(double));

        etaRS = malloc(ProLen*sizeof(double));

	for(c=0;c<2;++c)
		{
		BBStore[c] = malloc((ProLen-1)*sizeof(double*));
		SolStore[c] = malloc((ProLen-2)*sizeof(double**));

		for (i=0;i<(ProLen-1);++i)
			BBStore[c][i] =  malloc((ProLen-1)*sizeof(double*[20][3]));

		for (i=0;i<(ProLen-2);++i)
			{
			SolStore[c][i] =  malloc((NumSol)*sizeof(double*));

			num_atoms = (lenAA[RT[i+1]] < 5) ? 16 : 11+lenAA[RT[i+1]];

			for(j=0;j<NumSol;++j)
				SolStore[c][i][j] = malloc(num_atoms*sizeof(double*[3]));
			}
		}
        }


int getprotein(char *protein_file)
        {
        int i,j,k,c,atom_num;
        FILE *fp;
        char buf[1024];
        char aa;

        fp=fopen(protein_file,"r");
        if(!fp)
                {
                printf("protein_file not found\n");
                return 0;
                }

        fscanf(fp,"%d%*[ ]%d%*[\n]",&ProLen,&NumSol);

        RT=malloc(ProLen*sizeof(int));

        fgets(buf, sizeof buf, fp);

        fclose(fp);

        NumHPAA = 0;
        NumSPAA = 0;

        for(i=0;i<ProLen;++i)
                {
                RT[i] = AA[buf[i]];
                if (RT[i] >= 16)
                        NumHPAA++;
                if (RT[i] <= 7 )
                        NumSPAA++;
                }

        HPAA = malloc(NumHPAA*sizeof(int));
        SPAA = malloc(NumSPAA*sizeof(int));

        j = 0;
        k = 0;

	NumSPatoms = 0;

        for(i=0;i<ProLen;++i)
                {
                if (RT[i] >= 16)
                        {
                        HPAA[j] = i;
                        j++;
                        }

                if (RT[i] <= 7)
                        {
                        SPAA[k] = i;
                        k++;
			
			for(atom_num=5;atom_num<lenAA[RT[i]];++atom_num)
				if ((ATR[RT[i]][atom_num] == 1) || (ATR[RT[i]][atom_num] == 2))
					NumSPatoms++;
                        }
                }

        Wrank = malloc((ProLen-1+(2*NumSol))*sizeof(double*[2]));
        WrankB = malloc(((ProLen-1)*(ProLen-1))*sizeof(double*[3]));
        WrankS = malloc(((ProLen-2)*(NumSol))*sizeof(double*[3]));

        MatchArray = malloc((ProLen-1)*sizeof(int*[2]));

        WeightsBB[0] = malloc((ProLen-1)*sizeof(double*));
        WeightsBB[1] = malloc((ProLen-1)*sizeof(double*));

	for(c=0;c<2;++c)
		for(i=0;i<ProLen-1;++i)
			WeightsBB[c][i] = malloc((ProLen-1)*sizeof(double));

        WeightsSB[0] = malloc((ProLen-2)*sizeof(double*));
        WeightsSB[1] = malloc((ProLen-2)*sizeof(double*));

	for(c=0;c<2;++c)
		for(i=0;i<ProLen-2;++i)
			WeightsSB[c][i] = malloc(NumSol*sizeof(double));

        return 1;
        }


int getMotifsBB(char *motif_file)
        {
        int i,j,lineNum,resdiff,m,n,z,r1,r2,r3,r4;
	short int CurrCount[21][21][21][21][11];
        double coor1,coor2,coor3;
        FILE *fp;
        char buf[1024],a1,a2,a3,a4;

        fp=fopen(motif_file,"r");
        if(!fp)
                {
                printf("motif_file not found\n");
                return 0;
                }

        for(r1=0;r1<21;r1++)
        for(r2=0;r2<21;r2++)
        for(r3=0;r3<21;r3++)
        for(r4=0;r4<21;r4++)
                for(j=0;j<11;j++)
                        {
                        BBMotifNum[r1][r2][r3][r4][j] = 0;
                        CurrCount[r1][r2][r3][r4][j] = -1;
                        }

        while(fgets(buf, sizeof buf, fp))
                {
                if (strlen(buf)<20 && buf[0] != '\n')
                        {
                        sscanf(buf, "%c%*[ ]%c%*[ ]%c%*[ ]%c%*[ ]%d",&a1,&a2,&a3,&a4,&resdiff);

                        r1 = AA[a1];
                        r2 = AA[a2];
                        r3 = AA[a3];
                        r4 = AA[a4];

                        ++BBMotifNum[r1][r2][r3][r4][resdiff];
                        }
                }

        for(r1=0;r1<21;r1++)
        for(r2=0;r2<21;r2++)
        for(r3=0;r3<21;r3++)
        for(r4=0;r4<21;r4++)
                for(j=0;j<11;j++)
                        {
                        if (BBMotifNum[r1][r2][r3][r4][j] != 0)
                                BBMotifs[r1][r2][r3][r4][j] = malloc((BBMotifNum[r1][r2][r3][r4][j])*sizeof(double*[20][3]));
                        }

        rewind(fp);

        i = 0;
        while(fgets(buf, sizeof buf, fp))
                {
                if (strlen(buf)<20 && buf[0] != '\n')
                        {
                        sscanf(buf, "%c%*[ ]%c%*[ ]%c%*[ ]%c%*[ ]%d",&a1,&a2,&a3,&a4,&resdiff);

                        r1 = AA[a1];
                        r2 = AA[a2];
                        r3 = AA[a3];
                        r4 = AA[a4];

                        z=++CurrCount[r1][r2][r3][r4][resdiff];

                        lineNum = 0;

                        continue;
                        }

                else if (strlen(buf) > 20)
                        {
                        sscanf(buf, "%lf%*[ \t ]%lf%*[ \t ]%lf", &coor1, &coor2, &coor3);

                        BBMotifs[r1][r2][r3][r4][resdiff][z][lineNum][0] = coor1;
                        BBMotifs[r1][r2][r3][r4][resdiff][z][lineNum][1] = coor2;
                        BBMotifs[r1][r2][r3][r4][resdiff][z][lineNum][2] = coor3;

                        ++lineNum;
                        }
                }

        fclose(fp);

        return 1;
        }


int getMotifsSol(char *motif_file)
        {
        int i,j,lineNum,m,n,z,r1,r2,r3,num_atoms;
	short int CurrCount[21][21][21];
        double coor1,coor2,coor3;
        FILE *fp;
        char buf[1024],a1,a2,a3;

        fp=fopen(motif_file,"r");
        if(!fp)
                {
                printf("motif_file not found\n");
                return 0;
                }

        for(r1=0;r1<21;r1++)
        for(r2=0;r2<21;r2++)
        for(r3=0;r3<21;r3++)
                {
                SolMotifNum[r1][r2][r3] = 0;
                CurrCount[r1][r2][r3] = -1;
                }

        while(fgets(buf, sizeof buf, fp))
                {
                if (strlen(buf)<20 && buf[0] != '\n')
                        {
                        sscanf(buf, "%c%*[ ]%c%*[ ]%c",&a1,&a2,&a3);

                        r1 = AA[a1];
                        r2 = AA[a2];
                        r3 = AA[a3];

                        ++SolMotifNum[r1][r2][r3];
                        }
                }
        for(r1=0;r1<21;r1++)
        for(r2=0;r2<21;r2++)
        for(r3=0;r3<21;r3++)
                {
                if (SolMotifNum[r1][r2][r3] != 0)
			{
                        SolMotifs[r1][r2][r3] = malloc((SolMotifNum[r1][r2][r3])*sizeof(double*));
			num_atoms = (lenAA[r2] < 5) ? 16 : 11+lenAA[r2];
			for(i=0;i<SolMotifNum[r1][r2][r3];++i)
				SolMotifs[r1][r2][r3][i] = malloc(num_atoms*sizeof(double*[3]));
		
			}
                }

        rewind(fp);

        while(fgets(buf, sizeof buf, fp))
                {
                if (strlen(buf)<20 && buf[0] != '\n')
                        {
                        sscanf(buf, "%c%*[ ]%c%*[ ]%c",&a1,&a2,&a3);

                        r1 = AA[a1];
                        r2 = AA[a2];
                        r3 = AA[a3];

                        z=++CurrCount[r1][r2][r3];

                        lineNum = 0;

                        continue;
                        }

                else if (strlen(buf) > 20)
                        {
                        sscanf(buf, "%lf%*[ \t ]%lf%*[ \t ]%lf", &coor1, &coor2, &coor3);

                        SolMotifs[r1][r2][r3][z][lineNum][0] = coor1;
                        SolMotifs[r1][r2][r3][z][lineNum][1] = coor2;
                        SolMotifs[r1][r2][r3][z][lineNum][2] = coor3;

                        ++lineNum;
                        }
                }
        fclose(fp);

        return 1;
        }


int getMotifsRot(char *motif_file)
        {
        int i,j,lineNum,m,n,z,r1,motif_type,num_atoms;
	short int CurrCount[3][21];
        double coor1,coor2,coor3;
        FILE *fp;
        char buf[1024],a1,a2;

        fp=fopen(motif_file,"r");
        if(!fp)
                {
                printf("motif_file not found\n");
                return 0;
                }

	for(i=0;i<3;i++)
        for(r1=0;r1<21;r1++)
                {
                RotMotifNum[i][r1] = 0;
                CurrCount[i][r1] = -1;
                }

        while(fgets(buf, sizeof buf, fp))
                {
                if (strlen(buf)<5 && buf[0] != '\n')
                        {
                        sscanf(buf, "%c%*[ ]%d",&a1,&motif_type);

                        r1 = AA[a1];

                        ++RotMotifNum[motif_type][r1];
                        }
                }

	for(i=0;i<3;i++)
        for(r1=0;r1<21;r1++)
                {
                if (RotMotifNum[i][r1] != 0)
			{
                        RotMotifs[i][r1] = malloc((RotMotifNum[i][r1])*sizeof(double*));

			for(z=0;z<RotMotifNum[i][r1];++z)
				{
				if (i==0)
					num_atoms = 5+lenAA[r1];
				else if (i==1)
					num_atoms = 3+lenAA[r1];
				else
					num_atoms = 2+lenAA[r1];
	
				RotMotifs[i][r1][z] =  malloc((num_atoms)*sizeof(double*[3]));
				}
			}
                }

        rewind(fp);

        while(fgets(buf, sizeof buf, fp))
                {
                if (strlen(buf)<5 && buf[0] != '\n')
                        {
                        sscanf(buf, "%c%*[ ]%d",&a1,&motif_type);

                        r1 = AA[a1];
                        z=++CurrCount[motif_type][r1];
			lineNum = 0;

                        continue;
                        }

                else if (strlen(buf) > 5)
                        {
                        sscanf(buf, "%lf%*[ \t ]%lf%*[ \t ]%lf", &coor1, &coor2, &coor3);

                        RotMotifs[motif_type][r1][z][lineNum][0] = coor1;
                        RotMotifs[motif_type][r1][z][lineNum][1] = coor2;
                        RotMotifs[motif_type][r1][z][lineNum][2] = coor3;

                        ++lineNum;
                        }
                }

        fclose(fp);

        return 1;
        }


int getVEXinfo()
        {
        int i,j,k;
        FILE *fp;
        char buf[1024];

        fp = fopen("VEX_3.txt","r");
        if(!fp)
                {
                printf("VEX.txt not found\n");
                return 0;
                }

        for(i=0;i<4;++i)
                for(j=0;j<4;++j)
                        fscanf(fp,"%lf%*[ \t ]%*[\n]",&ATDMIN[i][j]);

        for(i=0;i<4;++i)
		fscanf(fp,"%lf%*[ \t ]%*[\n]",&ASDMIN[i]);

	fscanf(fp,"%lf%*[\n]",&OODMIN);

        fclose(fp);

        return 1;

        }


void AAmapping()
        {
        int i,j;

        for(i=65;i<=90;++i)
                AA[i]=0;

        AA['R']=1;
        AA['N']=2;
        AA['D']=3;
        AA['Q']=4;
        AA['E']=5;
        AA['H']=6;
        AA['K']=7;

        AA['C']=8;
        AA['M']=9;
        AA['G']=10;
        AA['P']=11;
        AA['S']=12;
        AA['T']=13;
        AA['W']=14;
        AA['Y']=15;

        AA['I']=16;
        AA['L']=17;
        AA['V']=18;
        AA['F']=19;
        AA['A']=20;

	/* number of atoms in all the amino acids*/
	
	lenAA[1] = 11; 
	lenAA[2] = 8; 
	lenAA[3] = 8; 
	lenAA[4] = 9; 
	lenAA[5] = 9; 
	lenAA[6] = 10; 
	lenAA[7] = 9; 

	lenAA[8] = 6; 
	lenAA[9] = 8; 
	lenAA[10] = 4; 
	lenAA[11] = 7; 
	lenAA[12] = 6; 
	lenAA[13] = 7; 
	lenAA[14] = 14; 
	lenAA[15] = 12; 

	lenAA[16] = 8; 
	lenAA[17] = 8; 
	lenAA[18] = 7; 
	lenAA[19] = 11; 
	lenAA[20] = 5; 

	/* atom types in all the amino acids. Atoms "CNOS" are indexed by 0,1,2,3 respectively */

	for(i=1;i<21;++i)
		{
		ATR[i] = malloc(lenAA[i]*sizeof(int));

		for(j=0;j<lenAA[i];++j)
			ATR[i][j] = 0;
		
		ATR[i][0] = 1;
		ATR[i][3] = 2;
		}
	
	ATR[AA['C']][5] = 3;
	ATR[AA['W']][8] = 1;
	ATR[AA['D']][6] = 2;
	ATR[AA['D']][7] = 2;
	ATR[AA['N']][6] = 2;
	ATR[AA['N']][7] = 1;
	ATR[AA['Y']][11] = 2;
	ATR[AA['S']][5] = 2;
	ATR[AA['Q']][7] = 2;
	ATR[AA['Q']][8] = 1;
	ATR[AA['E']][7] = 2;
	ATR[AA['E']][8] = 2;
	ATR[AA['K']][8] = 1;
	ATR[AA['R']][7] = 1;
	ATR[AA['R']][9] = 1;
	ATR[AA['R']][10] = 1;
	ATR[AA['T']][5] = 2;
	ATR[AA['M']][6] = 3;
	ATR[AA['H']][6] = 1;
	ATR[AA['H']][9] = 1;
        }


double urand(double a, double b)
        {
        double num;
        num = (b-a)*(((double)rand())/RAND_MAX)+a;
        return num;
        }


void init()
        {
        int i,j,k,l,m,n,c,num_atoms;
        FILE *fp;
        double a;

        a = 20.;

	for(i=0;i<ProLen;++i)
		{
		num_atoms = (lenAA[RT[i]] == 4) ? 5 : lenAA[RT[i]];
		for (j=0;j<num_atoms;++j)
                        for(k=0;k<3;++k)
                                atom[i][j][k] = urand(-a,a);
		}

        for(j=0;j<NumSol;++j)
                for(k=0;k<3;++k)
                        Solatom[j][k] = urand(-a,a);

        changeVar(BB, SB, VB, VS, RS, VO);

        // initializing etas

	for (c=0;c<2;++c)
		for (i=0;i<ProLen-1;++i)
			for (j=0;j<ProLen-1;++j)
				etaBB[c][i][j] = 1.;

	for (c=0;c<1;++c)
		for (i=0;i<ProLen-2;++i)
			etaSB[c][i] = 1.;

        for (i=0;i<ProLen;++i)
                for (j=i+1;j<ProLen;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				for(n=0;n<lenAA[RT[j]];++n)
					etaVB[i][j][m][n] = 1.;

        for (i=0;i<ProLen;++i)
                for (j=0;j<NumSol;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				etaVS[i][j][m] = 1.;

        for (i=0;i<NumSol;++i)
                for (j=0;j<NumSol;++j)
                        etaVO[i][j] = 1.;

        for (i=0;i<ProLen;++i)
		etaRS[i] = 1.;

	// Known matching
/*
        SBMatch = malloc(NumSol*sizeof(int*));
	
	for (i=0;i<NumSol;++i)
		SBMatch[i] = malloc((ProLen-2)*sizeof(int));

        fp = fopen("match.txt","r");
        if(!fp)
                printf("match_file not found\n");

	for (i=0;i<NumSol;++i)
		for(j=0;j<(ProLen-2);++j)
			SBMatch[i][j] = 0;

	for (i=0;i<NumSol;++i)
		{
		fscanf(fp,"%d%*[\n]",&j);
		for(k=0;k<j;++k)
			{
			fscanf(fp,"%d%*[ ]%*[\n]",&l);
			SBMatch[i][l] = 1;
			}
		}

        fclose(fp);
*/
        }


void initSol(double lambda)
        {
        int i,j,k,l,m,n,c,num_atoms;
        FILE *fp;
        char coorfile[50];

	ActualCoor = malloc(ProLen*sizeof(double*));

        for(i=0;i<ProLen;++i)
		{
		num_atoms = (lenAA[RT[i]] == 4) ? 5 : lenAA[RT[i]];
		ActualCoor[i] = malloc(num_atoms*sizeof(double*[3]));
		}

        ActualSol = malloc(NumSol*sizeof(double*[3]));

        sprintf(coorfile,"%s_coor.txt",ProID);

        fp = fopen(coorfile,"r");
        if(!fp)
                printf("coor_file not found\n");

	for(i=0;i<ProLen;++i)
		{
		num_atoms = (lenAA[RT[i]] == 4) ? 5 : lenAA[RT[i]];
		for (j=0;j<num_atoms;++j)
                        for(k=0;k<3;++k)
                                fscanf(fp,"%lf%*[ ]%*[\n]",&ActualCoor[i][j][k]);
		}

        for(j=0;j<NumSol;++j)
                for(k=0;k<3;++k)
                        fscanf(fp,"%lf%*[ ]%*[\n]",&ActualSol[j][k]);

        fclose(fp);

	for(i=0;i<ProLen;++i)
		{
		num_atoms = (lenAA[RT[i]] == 4) ? 5 : lenAA[RT[i]];
		for (j=0;j<num_atoms;++j)
                        for(k=0;k<3;++k)
                                atom[i][j][k] = ActualCoor[i][j][k]+lambda*urand(-1.,1.);
		}

        for(j=0;j<NumSol;++j)
                for(k=0;k<3;++k)
			Solatom[j][k] = ActualSol[j][k]+lambda*urand(-1.,1.);

        changeVar(BB, SB, VB, VS, RS, VO);

        // initializing etas

	for (c=0;c<2;++c)
		for (i=0;i<ProLen-1;++i)
			for (j=0;j<ProLen-1;++j)
				etaBB[c][i][j] = 1.;

	for (c=0;c<2;++c)
		for (i=0;i<ProLen-2;++i)
			etaSB[c][i] = 1.;

        for (i=0;i<ProLen;++i)
                for (j=i+1;j<ProLen;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				for(n=0;n<lenAA[RT[j]];++n)
					etaVB[i][j][m][n] = 1.;

        for (i=0;i<ProLen;++i)
                for (j=0;j<NumSol;++j)
			for(m=0;m<lenAA[RT[i]];++m)
				etaVS[i][j][m] = 1.;

        for (i=0;i<NumSol;++i)
                for (j=0;j<NumSol;++j)
                        etaVO[i][j] = 1.;

        for (i=0;i<ProLen;++i)
		etaRS[i] = 1.;
		
	// Known matching
/*
        SBMatch = malloc(NumSol*sizeof(int*));
	
	for (i=0;i<NumSol;++i)
		SBMatch[i] = malloc((ProLen-2)*sizeof(int));

        fp = fopen("match.txt","r");
        if(!fp)
                printf("match_file not found\n");

	for (i=0;i<NumSol;++i)
		for(j=0;j<(ProLen-2);++j)
			SBMatch[i][j] = 0;

	for (i=0;i<NumSol;++i)
		{
		fscanf(fp,"%d%*[\n]",&j);
		for(k=0;k<j;++k)
			{
			fscanf(fp,"%d%*[ ]%*[\n]",&l);
			SBMatch[i][l] = 1;
			}
		}

        fclose(fp);
*/	
	}


int main(int argc,char* argv[])
        {
        char *id,MotifFileBB[20],MotifFileSol[20],MotifFileRot[20],buf[1024];
        int c,maxiter,iterstride,iter,i,j,m,n,z,p,q,r,seed;
        double deltat, stoperr, iterpersec, peta;
        FILE *fp;
        clock_t start;

        if(argc==14)
                {
                ProID = argv[1];
                fracSP = atof(argv[2]);
                fracOS = atof(argv[3]);
                totHBB = atoi(argv[4]);
                totHBS = atoi(argv[5]);
                CHstep = atoi(argv[6]);
                id = argv[7];
                beta = atof(argv[8]);
                maxiter = atoi(argv[9]);
                iterstride = atoi(argv[10]);
                stoperr = atof(argv[11]);
                epsilon = atof(argv[12]);
                seed = atoi(argv[13]);
                }
        else
                {
                printf("expected more/less arguments \n");
                return 1;
                }

        sprintf(errfile,"%s.err",id);
        sprintf(statsfile,"%s.stats",id);
        sprintf(solfile,"%s.sol",id);
        sprintf(CAfile,"%s.ca",id);

        AAmapping();
        sprintf(ProFile,"%s.txt",ProID);

        if(!getprotein(ProFile))
                return 1;

        fp=fopen(solfile,"w");
        fclose(fp);

        fp=fopen(CAfile,"w");
        fclose(fp);

        fp=fopen(statsfile,"w");
        for(c=0;c<argc;++c)
                fprintf(fp,"%s ",argv[c]);
        fprintf(fp,"\n\n");

        fprintf(fp,"Protein Size: %d\n\n",ProLen);
        fprintf(fp,"Residue Types:\n");

        for(c=0;c<ProLen;++c)
                fprintf(fp,"%d ",RT[c]);
        fprintf(fp,"\n\n");

        fclose(fp);

        sprintf(MotifFileBB,"%s_BBv1.txt",ProID);
        if(!getMotifsBB(MotifFileBB))
                return 1;

        sprintf(MotifFileSol,"%s_Solv1.txt",ProID);
        if(!getMotifsSol(MotifFileSol))
                return 1;

        sprintf(MotifFileRot,"%s_Rotv1.txt",ProID);
        if(!getMotifsRot(MotifFileRot))
                return 1;

        if(!getVEXinfo())
                return 1;

        makevars();
        srand(seed);

        init();

        //initSol(0.0);

/*
        changeVar(BB, SB, VB, VS, RS, VO);
	projA(BB, SB, VB, VS, RS, VO);
	reflect(BBA, SBA, VBA, VSA, RSA, VOA);
	projB(BBR, SBR, VBR, VSR, RSR, VOR);
        RRR();
        iter=solve(maxiter,iterstride,stoperr);

*/
        start=clock();
        iter=solve(maxiter,iterstride,stoperr);

        deltat=((double)(clock()-start))/CLOCKS_PER_SEC;

        fp=fopen(statsfile,"a");
        if(iter)
                {
                iterpersec=iter/deltat;
                fprintf(fp,"\nTime elapsed: %10.2lf \nNumber of iterations: %d \n",deltat,iter);
                }
        else
                {
                iterpersec=maxiter/deltat;
                fprintf(fp,"\nTime elapsed: %10.2lf \n",deltat);
                }

        fprintf(fp,"iterations/sec:%10.2lf\n",iterpersec);

        fprintf(fp,"\n\n Backbone-Backbone Matching:\n\n");

	for(r=0;r<totHBB;++r)
                {
                p = (int) (WrankB[r][0]+.5);
                q = (int) (WrankB[r][1]+.5);

		fprintf(fp,"%d %d %.6e %.6e\n",p,q,etaBB[1][p][q],WeightsBB[1][p][q]);
		}

        fclose(fp);

/*
        fprintf(fp,"\n\n Backbone-Solvent Matching:\n\n");

	for(r=0;r<totHBS;++r)
                {
                p = (int) (WrankS[r][0]+.5);
                q = (int) (WrankS[r][1]+.5);

		fprintf(fp,"%d %d %.6e %.6e\n",p,q,etaSB[1][p],WeightsSB[1][p][q]);
		}

        for(i=0;i<ProLen-1;++i)
                {

                if ( (RT[i]==10 || RT[i+1]==10) || (RT[i]==11 || RT[i+1]==11))
                        {
                        m = MatchArray[i][0];
                        peta = (m < (ProLen-1)) ? etaBB[0][i][m] : etaSB[0][i][m-ProLen+1];
                        fprintf(fp,"%d %d %.6e %.6e\n",i,m,peta,WeightsBB[0][i][m]);
                        }

                else
                        {
                        m = MatchArray[i][0];
                        n = MatchArray[i][1];

                        peta = (m < (ProLen-1)) ? etaBB[0][i][m] : etaSB[0][i][m-ProLen+1];
                        fprintf(fp,"%d %d %.6e %.6e\n",i,m,peta,WeightsBB[0][i][m]);

                        peta = (n < (ProLen-1)) ? etaBB[0][i][n] : etaSB[0][i][n-ProLen+1];
                        fprintf(fp,"%d %d %.6e %.6e\n",i,n,peta,WeightsBB[0][i][n]);
                        }
                }
*/

        return 0;
        }
