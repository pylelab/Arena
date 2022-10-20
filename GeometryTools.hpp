/*************************************************************************
Geometry Functions for Structural Bioinformatics
Author: Cao Yang
contact: cao@scu.edu.cn
Date:2008.7.11.
*************************************************************************/

#ifndef Geometryhpp
#define Geometryhpp

#include <vector>
#include <cmath>
#include <iostream>
#include <cstdlib> 
using namespace std;

const double PI=3.141592653589793;
const double Extra=1.0e-4;
const double UpMax=1.0e+10;
typedef vector<vector<double> > matrix;


void ShowMyvector(const vector<double> &cc)
{
   for(int i=0; i<cc.size(); ++i)
   {
      cerr<<cc[i]<<'\t';
   }
   cerr<<endl;
}

//Functions using vectors
inline void crossproduct(const vector<double> &c1, const vector<double> &c2, vector<double> &cc)
{
   cc[0] = c1[1] * c2[2] - c1[2] * c2[1];
   cc[1] = c1[2] * c2[0] - c1[0] * c2[2];
   cc[2] = c1[0] * c2[1] - c2[0] * c1[1];
}

inline double innerproduct(const vector<double> &c1, const vector<double> &c2)
{
   return c1[0] * c2[0] + c1[1] * c2[1] + c1[2] * c2[2];
}

inline double innerproduct_univ(const vector<double> &c1, const vector<double> &c2)
{
	double sum=0;
	for(int i=0; i<c1.size(); i++)
		sum += c1[i] * c2[i];
	return sum;
}

inline void vectorsum(const vector<double> &c1, const vector<double> &c2, vector<double> &cc)
{
   cc[0] = c1[0] + c2[0];
   cc[1] = c1[1] + c2[1];
   cc[2] = c1[2] + c2[2];
}

inline void vectorsum_univ(const vector<double> &c1, const vector<double> &c2, vector<double> &cc)
{
	for(int i=0; i<cc.size(); i++)
		cc[i] = c1[i] + c2[i];
}

inline void subtract(const vector<double> &c1, const vector<double> &c2, vector<double> &cc)
{
   cc[0] = c1[0] - c2[0];
   cc[1] = c1[1] - c2[1];
   cc[2] = c1[2] - c2[2];
}

inline void subtract_univ(const vector<double> &c1, const vector<double> &c2, vector<double> &cc)
{
	for(int i=0; i<cc.size(); i++)
		cc[i] = c1[i] - c2[i];
}

inline bool norm(const vector<double> &c, vector<double> &cc)
{
   double len = c[0]*c[0] + c[1]*c[1] +c[2]*c[2];
   if(len<Extra) 
   {
      cerr<<"Error in norm()! length~=0!\n";
      ShowMyvector(c);
      ShowMyvector(cc);
      return false;
   }
   len = 1/sqrt(len);
   cc[0] = c[0]*len;
   cc[1] = c[1]*len;
   cc[2] = c[2]*len;
   return true;
}

inline bool norm_univ(const vector<double> &c, vector<double> &cc)
{
	double len = 0;
	for(int i=0; i<c.size(); i++)
		len+=c[i]*c[i];
	if(len<Extra) 
	{
		cerr<<"Error in norm()! length~=0!\n";
		return false;
	}
	len = 1/sqrt(len);
	for(int i=0; i<c.size(); i++)
		cc[i] = c[i]*len;
	return true;
}

inline bool norm_warnless(vector<double> &c)
{
   double len = c[0]*c[0] + c[1]*c[1] +c[2]*c[2];
   if(len<Extra) 
   {
      return false;
   }
   len = 1/sqrt(len);
   c[0] *= len;
   c[1] *= len;
   c[2] *= len;
   return true;
}

inline bool norm_univ_warnless(vector<double> &c)
{
	double len = 0;
	for(int i=0; i<c.size(); i++)
		len+=c[i]*c[i];
	if(len<Extra) 
	{
		return false;
	}
	len = 1/sqrt(len);
	for(int i=0; i<c.size(); i++)
		c[i]*=len;
	return true;
}

inline double vectorlength(const vector<double> &c)
{
   return sqrt(c[0]*c[0]+ c[1]*c[1] + c[2]*c[2]);
}

inline double vectorlength_univ(const vector<double> &c)
{
	double len = 0;
	for(int i=0; i<c.size(); i++)
		len+=c[i]*c[i];
	return sqrt(len);
}

inline double vectorlength2(const vector<double> &c)
{
   return (c[0]*c[0]+ c[1]*c[1] + c[2]*c[2]);
}

inline double vectorlength2_univ(const vector<double> &c)
{
	double len = 0;
	for(int i=0; i<c.size(); i++)
		len+=c[i]*c[i];
	return len;
}

inline void multi(double coefficient, vector<double> &c, vector<double> &cc)
{
   cc[0] = c[0]*coefficient;
   cc[1] = c[1]*coefficient;
   cc[2] = c[2]*coefficient;
}

inline void multi_univ(double coefficient, vector<double> &c, vector<double> &cc)
{
	for(int i=0; i<c.size(); i++)
		cc[i] = c[i] * coefficient;
}

inline void multi(double coefficient, vector<double> &cc)
{
   cc[0] *= coefficient;
   cc[1] *= coefficient;
   cc[2] *= coefficient;
}

inline void multi_univ(double coefficient, vector<double> &cc)
{
	for(int i=0; i<cc.size(); i++)
		cc[i]*=coefficient;
}

inline double deg2rad(double deg)
{
   return deg*PI/180;
}

inline double rad2deg(double rad)
{
   return rad*180/PI;
}


inline double Points2Distance2(const vector<double> &c1, const vector<double> &c2)
{
   double a=(c1[0]-c2[0]);
   double b=(c1[1]-c2[1]);
   double c=(c1[2]-c2[2]);
   return a*a+b*b+c*c;
}

inline double Points2Distance(const vector<double> &c1, const vector<double> &c2)
{
   double a=(c1[0]-c2[0]);
   double b=(c1[1]-c2[1]);
   double c=(c1[2]-c2[2]);
   return sqrt(a*a+b*b+c*c);
}

//Return the angle of c1-c2-c3. Unit: radian
//Angle <c1c2c3 
inline double Points2Angle(const vector<double> &c1, const vector<double> &c2, const vector<double> &c3)
{
   double a, b, c, tmp1, tmp2, tmp3, alpha;
   
   tmp1=c1[0]-c2[0];
   tmp2=c1[1]-c2[1];
   tmp3=c1[2]-c2[2];
   a=tmp1*tmp1+tmp2*tmp2+ tmp3*tmp3;
   
   tmp1=c2[0]-c3[0];
   tmp2=c2[1]-c3[1];
   tmp3=c2[2]-c3[2];
   b=tmp1*tmp1+tmp2*tmp2+ tmp3*tmp3;
   
   tmp1=c3[0]-c1[0];
   tmp2=c3[1]-c1[1];
   tmp3=c3[2]-c1[2];
   c=tmp1*tmp1+tmp2*tmp2+ tmp3*tmp3;
   
   alpha=(a+b-c)/(2*sqrt(a*b));
   
   if(alpha>1+Extra || alpha<-1-Extra)
      cerr<<"Warning, double Points2Angle()\n";

   if (alpha>1) alpha=1;
   else if(alpha<-1) alpha=-1;
   return acos(alpha);
}

//Return the angle of <c1-c2,c3-c4>. Unit: radian
inline double Points4Angle(const vector<double> &c1, const vector<double> &c2, const vector<double> &c3, const vector<double> &c4)
{
    double a, b, a1, a2, a3, b1, b2, b3, alpha;
   
    a1=c2[0]-c1[0];
    a2=c2[1]-c1[1];
    a3=c2[2]-c1[2];
    a = a1*a1 + a2*a2 + a3*a3;
   
    b1=c4[0]-c3[0];
    b2=c4[1]-c3[1];
    b3=c4[2]-c3[2];
    b = b1*b1 + b2*b2 + b3*b3;
   
    alpha=(a1*b1+a2*b2+a3*b3)/sqrt(a*b);
   
    if(alpha>1+Extra || alpha<-1-Extra)
        cerr<<"Error, double Points2Angle()\n";
   
    if (alpha>1)      alpha=1;
    else if(alpha<-1) alpha=-1;
    return acos(alpha);
}

/************************Points to Dihedral************************/
//If use inline here there will be a error when combiling codes.
//Return the angle of c1-c2-c3-c4. Unit: radian
//points c1-c2-c3-c4 should be in two different pane.
//or there may be faults.
//by Yang Cao
//
//Change by chengxin: if they are on the same plane, return -2*PI
inline double Points2Dihedral(const vector<double> &c1, const vector<double> &c2, const vector<double> &c3, const vector<double> &c4)
{
   vector<double> vector1(3,0), vector2(3,0), vector3(3,0);
   subtract(c1, c2, vector1);
   subtract(c2, c3, vector2);
   subtract(c3, c4, vector3);
   
   vector<double> v1(3,0), v2(3,0);
   crossproduct(vector2, vector1, v1);
   crossproduct(vector3, vector2, v2);
   
   vector<double> v3(3,0), v4(3,0);
   if(!norm (v1, v3)) {cerr<<"Error in Points2Dihedral 1\n"<<endl; return -2*PI;}//exit(1);}
   if(!norm (v2, v4)) {cerr<<"Error in Points2Dihedral 2\n"<<endl; return -2*PI;}//exit(1);}
   
   double dihedral = innerproduct(v3, v4);
   
   if (dihedral>1 && dihedral<1+Extra)
   {
      //cerr<<"dihedral "<<dihedral<<" "<<acos(dihedral)<<"\n";
      dihedral=1;
   }
   else if(dihedral<-1 && dihedral>-1-Extra)
   {
      dihedral=-1;
   }
   else if(dihedral>1+Extra || dihedral<-1-Extra)
   {
      cerr<<"Error, double Points2Dihedral()\n";
      exit(0);
   }
   
   vector<double> v5(3,0);
   crossproduct(v4, v3, v5);
   double direction = innerproduct(v5, vector2);
   
   if (direction>0)
    {
       return  acos(dihedral);
    }  
    else
    {  
       return -acos(dihedral);
    }
   
}
//////////////////////////////////Matrix Tool/////////////////////////////////

//Method to initialize a matrix
inline void SetMatrix(matrix &sm, int m, int n)
{
   vector<double> tmp(n, 0);
   sm.assign(m, tmp);
}

inline void MatrixTimesMatrix(const matrix &a, const matrix &v, matrix &x, int m, int n, int l)
{
   int i=0, j=0, k=0;
   double sum=0;
   for(i=0; i<m; ++i)
   {
      for(k=0; k<l; ++k)
      {
         for(j=0, sum=0; j<n; ++j)
         {
            sum+=a[i][j]*v[j][k];
         }
         x[i][k]=sum;
      }
   }
}

inline bool TransVectorTimesVector(const vector<double> &trans, const vector<double> &vctor, matrix &mtx)
{
   if(trans.size()!=vctor.size())
   {
      cerr<<"Error in TransVectorTimesVector()"<<endl;
      return false;
   }
   int i=0, j=0;
   SetMatrix(mtx, trans.size(), vctor.size());
   for(i=0; i<trans.size(); ++i)
      for(j=0; j<vctor.size(); ++j)
         mtx[i][j]=trans[i]*vctor[j];
   
   return true;
}

inline bool MatrixTimesTransVector(const matrix &mtx, const vector<double> &tvt, vector<double> &vct)
{
   if(mtx[0].size()!=tvt.size())
   {
      cerr<<"Error in MatrixTimesTransVector()"<<endl;
      return false;
   }
   int i=0, j=0;
   vct.assign(mtx.size(), 0);
   for(i=0; i<mtx.size(); ++i)
      for(j=0; j<mtx[i].size(); ++j)
         vct[i]+=mtx[i][j]*tvt[j];
   
   return true;
}

inline void RealTimesMatrix(double at, const matrix &mx, matrix &mc)
{
   int i=0, j=0;
   for(i=0; i<mx.size(); ++i)
      for(j=0; j<mx[i].size(); ++j)
         mc[i][j]=at*mx[i][j];
}

inline bool MatrixAddMatrix(const matrix &ma, const matrix &mb, matrix &mc)
{
   if(ma.size()!=mb.size() || ma[0].size()!=mb[0].size()
   ||ma.size()!=mc.size() || ma[0].size()!=mc[0].size())
   {
      cerr<<"Error size! MatrixAddMatrix()"<<endl;
      return false;
   }
   int i=0, j=0;
   for(i=0; i<ma.size(); ++i)
      for(j=0; j<ma[i].size(); ++j)
         mc[i][j]=ma[i][j]+mb[i][j];
    
   return true;     
}

inline bool norm2(const vector<double> &c, vector<double> &cc)
{
   if(c.size()!=cc.size())
   {
      cerr<<"Error in norm() "<<endl;
      ShowMyvector(c);
      ShowMyvector(cc);
      return false;
   }
   int i=0; 
   double len=0;
   for(i=0; i<c.size(); ++i)
      len+=c[i]*c[i];
   len=sqrt(len);

   if(len<Extra) return false;
   for(i=0; i<c.size(); ++i)
      cc[i]=c[i]/len;
   
   return true;
   
}

inline bool norm2_warnless(vector<double> &c)
{
   int i=0; 
   double len=0;
   for(i=0; i<c.size(); ++i)
      len+=c[i]*c[i];
   len=1/sqrt(len);

   if(len<Extra) return false;
   for(i=0; i<c.size(); ++i)
      c[i]*=len;
   
   return true;
   
}

/************************Rotation Matrix****************************
* Give the Rotation Axis and the Rotation Angle
* Generate the Rotation Matrix
* Author: C.Y.
* Date 2006.11.30.
* BUG: romtx[0][2]=0 should be romtx[3][2]=0
*************************End***************************************/
inline bool RotationMatrixA(const vector<double> &axis, double angle, matrix &romtx)
{
   if(axis.size()!=3) 
   {
      cerr<<"Error 1 in RotationMatrixA()"<<endl;
      return false;
   }
   vector<double> ouc(3, 0);
   if(!norm(axis, ouc)) 
   {
      cerr<<"Error 2 in RotationMatrixA()"<<endl;
      return false;
   }
   int i=0, j=0;
   matrix su, u_i, unit, tmp;
   
   SetMatrix(su, 3, 3);
   su[0][0]=0; su[0][1]=-ouc[2]; su[0][2]=ouc[1]; 
   su[1][0]=ouc[2]; su[1][1]=0; su[1][2]=-ouc[0]; 
   su[2][0]=-ouc[1]; su[2][1]=ouc[0]; su[2][2]=0;
   
   RealTimesMatrix(sin(angle), su, su);
   /* The method below is trivial.
   matrix u_t, u_c;
   SetMatrix(u_t, 3, 1);
   SetMatrix(u_c, 1, 3);
   SetMatrix(u_i, 3, 3);
   u_t[0][0]=ouc[0]; u_t[1][0]=ouc[1]; u_t[2][0]=ouc[2];
   u_c[0][0]=ouc[0]; u_c[0][1]=ouc[1]; u_c[0][2]=ouc[2]; 
   MatrixTimesMatrix(u_t, u_c, u_i, 3, 1, 3);*/
   TransVectorTimesVector(ouc, ouc, u_i);
   
   SetMatrix(unit, 3, 3);
   SetMatrix(tmp, 3, 3);
   unit[0][0]=1; unit[1][1]=1; unit[2][2]=1;
   RealTimesMatrix(-1, u_i, tmp);
   MatrixAddMatrix(unit, tmp, tmp);
   RealTimesMatrix(cos(angle), tmp, tmp); 
   
   MatrixAddMatrix(u_i, tmp, tmp);
   MatrixAddMatrix(su, tmp, tmp);
   
   SetMatrix(romtx, 4, 4);
   for(i=0; i<tmp.size(); ++i)
      for(j=0; j<tmp[i].size(); ++j)
         romtx[i][j]=tmp[i][j];
   romtx[0][3]=0; romtx[1][3]=0; romtx[2][3]=0; 
   romtx[3][0]=0; romtx[3][1]=0; romtx[3][2]=0; 
   romtx[3][3]=1;
   
   return true;
}
/************************Rotation Matrix****************************
* Give the Rotation Axis and the Rotation Angle
* Generate the Rotation Matrix
* Author: C.Y.
* Date 2006.12.6.
*************************End***************************************/
inline bool RotationMatrixB(const vector<double> &axis, double angle, matrix &romtx)
{
   if(axis.size()!=3) 
   {
      cerr<<"Error 1 in RotationMatrixB()"<<endl;
      return false;
   }
   vector<double> ouc(3, 0);
   if(!norm(axis, ouc)) 
   {
      cerr<<"Error 2 in RotationMatrixB()"<<endl;
      return false;
   }
   double c=cos(angle);
   double s=sin(angle);
   double t=1-c;
   SetMatrix(romtx, 4, 4);
   romtx[0][0]=t*ouc[0]*ouc[0]+c;romtx[0][1]=t*ouc[0]*ouc[1]+s*ouc[2];
   romtx[0][2]=t*ouc[0]*ouc[2]-s*ouc[1];romtx[0][3]=0;
   
   romtx[1][0]=t*ouc[1]*ouc[0]-s*ouc[2];romtx[1][1]=t*ouc[1]*ouc[1]+c;
   romtx[1][2]=t*ouc[1]*ouc[2]+s*ouc[0];romtx[1][3]=0;
   
   romtx[2][0]=t*ouc[2]*ouc[0]+s*ouc[1];romtx[2][1]=t*ouc[2]*ouc[1]-s*ouc[0];
   romtx[2][2]=t*ouc[2]*ouc[2]+c;romtx[2][3]=0;
   
   romtx[3][0]=0;romtx[3][1]=0;romtx[3][2]=0;romtx[3][3]=1;
   
   return true;
}

/***************************CoordinateRotation********************
* Give the Rotation Axis and the Rotation Angle for a PointA
* Generate the coordinate of PointB after operation
* Author: C.Y.
* Date 2006.12.2.  
*****************************END**********************************/
bool CoordinateRotation(vector<double> &pointA, const vector<double> &axis, double angle, vector<double> &pointB)
{
   if(pointA.size()!=3)
   {
      cerr<<"Error in CoordinateRotation()"<<endl;
      return false;
   }
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   pointA.push_back(1);
   if(!MatrixTimesTransVector(rotmtx, pointA, pointB))
   { 
      pointA.pop_back(); 
      return false;
   }
   pointA.pop_back(); pointB.pop_back();
   return true;
}

/***************************CoordinateRotation********************
* For a PointA, Given the Coordinates of 2 points(axisA, axisB) 
* in line of Rotation Axis and the Rotation Angle 
* Generate the coordinate of PointB after operation
* Author: C.Y.
* Date 2006.12.2.  
*
*****************************END**********************************/
bool CoordinateRotation(const vector<double> &pointA, const vector<double> &axisA, const vector<double> &axisB, double angle, vector<double> &pointB)
{
   if(pointA.size()!=3 || axisA.size()!=3 || axisB.size()!=3)
   {
      cerr<<"Error in CoordinateRotation()"<<endl;
      return false;
   }
   
   vector<double> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   vector<double> point_A(4, 1);
   point_A[0]=pointA[0]-axisA[0];
   point_A[1]=pointA[1]-axisA[1];
   point_A[2]=pointA[2]-axisA[2];
   
   if(!MatrixTimesTransVector(rotmtx, point_A, pointB))
   { 
      point_A.pop_back(); 
      return false;
   }
   point_A.pop_back(); pointB.pop_back();
   
   pointB[0]+=axisA[0];
   pointB[1]+=axisA[1];
   pointB[2]+=axisA[2];
   
   return true;
}

/******************************************************************************
Coordinate Rotation for a group of points
CaoYang
2008.6.14.
******************************************************************************/
bool GroupRotation(const vector<double> &axisA, const vector<double> &axisB, double angle, vector<vector<double> > &pointB)
{
   if(axisA.size()!=3 || axisB.size()!=3)
   {
      cerr<<"Error in GroupRotation()"<<endl;
      return false;
   }
   
   vector<double> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   //vector<double> point_A(4, 0);
   double point_A[3];
   for(int i=0; i<pointB.size(); ++i)
   {
      point_A[0]=pointB[i][0]-axisA[0];
      point_A[1]=pointB[i][1]-axisA[1];
      point_A[2]=pointB[i][2]-axisA[2];

      pointB[i]=axisA;
      pointB[i][0]+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
      pointB[i][1]+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
      pointB[i][2]+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
   }
   return true;
}
//CoordinateRotation for a given group of points
bool GroupRotation(const vector<double> &axisA, const vector<double> &axisB, double angle, vector<vector<double> > &pointB, const vector<short>& index)
{
   if(axisA.size()!=3 || axisB.size()!=3)
   {
      cerr<<"Error in GroupRotation()"<<endl;
      return false;
   }
   
   vector<double> axis(3, 0);
   axis[0]=axisB[0]-axisA[0]; 
   axis[1]=axisB[1]-axisA[1];
   axis[2]=axisB[2]-axisA[2]; 
   
   matrix rotmtx;
   if(!RotationMatrixB(axis, deg2rad(angle), rotmtx)) return false;
   
   double point_A[3];
   int i, j, m, n;
   for(j=0; j<index.size(); j++)
   {
      i=index[j];
      point_A[0]=pointB[i][0]-axisA[0];
      point_A[1]=pointB[i][1]-axisA[1];
      point_A[2]=pointB[i][2]-axisA[2];

      pointB[i]=axisA;
      pointB[i][0]+=rotmtx[0][0]*point_A[0]+rotmtx[0][1]*point_A[1]+rotmtx[0][2]*point_A[2]+rotmtx[0][3];
      pointB[i][1]+=rotmtx[1][0]*point_A[0]+rotmtx[1][1]*point_A[1]+rotmtx[1][2]*point_A[2]+rotmtx[1][3];
      pointB[i][2]+=rotmtx[2][0]*point_A[0]+rotmtx[2][1]*point_A[1]+rotmtx[2][2]*point_A[2]+rotmtx[2][3];
   
   }
   return true;
}


/******************************************************************************
Coordinate Translation for a group of points
2008.6.14.
******************************************************************************/
bool GroupTranslation(const vector<double> &trans, vector<vector<double> > &pointB)
{
   if(trans.size()!=3)
   {
      cerr<<"Error in GroupTranslation\n";
      return false;
   }
   for(int i=0; i<pointB.size(); ++i)
   {
      pointB[i][0]+=trans[0];
      pointB[i][1]+=trans[1];
      pointB[i][2]+=trans[2];
   }
   return true;
}

#endif


