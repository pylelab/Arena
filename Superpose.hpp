/*********************************************************
Compute the best rotation matrix and translation vector that optimize 
overlap between two sets of coordinates 
Date: 20100118
*********************************************************/
#include <vector>
#include <cmath>
#include <cstdlib>

using namespace std;

/* vector functions using c arrays */

#define ROTATE(a,i,j,k,l) { g = a[i][j]; \
                            h = a[k][l]; \
                            a[i][j] = g-s*(h+g*tau); \
                            a[k][l] = h+s*(g-h*tau); }
/*   
 * jacobi3
 *
 *    computes eigenval and eigen_vec of a real 3x3
 * symmetric matrix. On output, elements of a that are above 
 * the diagonal are destroyed. d[1..3] returns the 
 * eigenval of a. v[1..3][1..3] is a matrix whose 
 * columns contain, on output, the normalized eigen_vec of
 * a. n_rot returns the number of Jacobi rotations that were required.
 */
int jacobi3(float a[3][3], float d[3], float v[3][3], int* n_rot)
{
  int count, k, i, j;
  float tresh, theta, tau, t, sum, s, h, g, c, b[3], z[3];

  /*Initialize v to the identity matrix.*/
  for (i=0; i<3; i++) 
  { 
    for (j=0; j<3; j++) 
      v[i][j] = 0.0;
    v[i][i] = 1.0;
  }

  /* Initialize b and d to the diagonal of a */
  for (i=0; i<3; i++) 
    b[i] = d[i] = a[i][i];

  /* z will accumulate terms */
  for (i=0; i<3; i++) 
    z[i] = 0.0; 
  
  *n_rot = 0;

  /* 50 tries */
  for (count=0; count<50; count++)     
  {

    /* sum off-diagonal elements */
    sum = 0.0;
    for (i=0; i<2; i++) 
    {
      for (j=i+1; j<3; j++)
         sum += fabs(a[i][j]);
    }

    /* if converged to machine underflow */
    if (sum == 0.0) 
      return(1);

    /* on 1st three sweeps... */
    if (count < 3) 
      tresh = sum * 0.2 / 9.0;    
    else       
      tresh = 0.0;      

    for (i=0; i<2; i++) 
    {
      for (j=i+1; j<3; j++) 
      {
        g = 100.0 * fabs(a[i][j]);

        /*  after four sweeps, skip the rotation if
         *   the off-diagonal element is small 
         */
        if ( count > 3  &&  fabs(d[i])+g == fabs(d[i])
              &&  fabs(d[j])+g == fabs(d[j]) ) 
        {
          a[i][j] = 0.0;
        } 
        else if (fabs(a[i][j]) > tresh) 
        {
          h = d[j] - d[i];
          
          if (fabs(h)+g == fabs(h))
          {
            t = a[i][j] / h;
          }
          else 
          {
            theta = 0.5 * h / (a[i][j]);
            t = 1.0 / ( fabs(theta) +
                        (float)sqrt(1.0 + theta*theta) );
            if (theta < 0.0) 
              t = -t;
          }
          
          c = 1.0 / (float) sqrt(1 + t*t);
          s = t * c;
          tau = s / (1.0 + c);
          h = t * a[i][j];

          z[i] -= h;
          z[j] += h;
          d[i] -= h;
          d[j] += h;

          a[i][j] = 0.0;

          for (k=0; k<=i-1; k++) 
            ROTATE(a, k, i, k, j)

          for (k=i+1; k<=j-1; k++) 
            ROTATE(a, i, k, k, j)

          for (k=j+1; k<3; k++) 
            ROTATE(a, i, k, j, k)

          for (k=0; k<3; k++) 
            ROTATE(v, k, i, k, j)

          ++(*n_rot);
        }
      }
    }

    for (i=0; i<3; i++) 
    {
      b[i] += z[i];
      d[i] = b[i];
      z[i] = 0.0;
    }
  }

  printf("Too many iterations in jacobi3\n");
  return (0);
}  

void normalize(float a[3])
{
  float  b = 1/sqrt((float)(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]));
  a[0] *= b;
  a[1] *= b;
  a[2] *= b;
}

float dot_prod(float a[3], float b[3])
{
  return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

static void cross(float a[3], float b[3], float c[3])
{
  a[0] = b[1]*c[2] - b[2]*c[1];
  a[1] = b[2]*c[0] - b[0]*c[2];
  a[2] = b[0]*c[1] - b[1]*c[0];
}

/* 
 * diagonalize_symmetric 
 *
 *    Diagonalize a 3x3 matrix & sort eigenval by size
 */
int diagonalize_symmetric(float matrix[3][3], 
                          float eigen_vec[3][3], 
                          float eigenval[3])
{
  int n_rot, i, j, k;
  float vec[3][3];
  float val; 
  
  if (!jacobi3(matrix, eigenval, vec, &n_rot)) 
  {
    printf("convergence failed\n");
    return (0);
  }

  /* sort solutions by eigenval */
  for (i=0; i<3; i++) 
  {
    k = i;
    val = eigenval[i];
    
    for (j=i+1; j<3; j++)
      if (eigenval[j] >= val)
      { 
        k = j;
        val = eigenval[k];
      }
       
    if (k != i) 
    {
      eigenval[k] = eigenval[i];
      eigenval[i] = val;
      for (j=0; j<3; j++) 
      {
        val = vec[j][i];
        vec[j][i] = vec[j][k];
        vec[j][k] = val;
      }
    }
  }

  /* transpose such that first index refers to solution index */
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      eigen_vec[i][j] = vec[j][i];

  return (1);
}


void VectorTimesMatrix(float *v1, const vector<vector<float> > &rotation_matrix,  float *v2)
{
	v2[0]=v1[0]*rotation_matrix[0][0]+v1[1]*rotation_matrix[1][0]+v1[2]*rotation_matrix[2][0];
	v2[1]=v1[0]*rotation_matrix[0][1]+v1[1]*rotation_matrix[1][1]+v1[2]*rotation_matrix[2][1];
	v2[2]=v1[0]*rotation_matrix[0][2]+v1[1]*rotation_matrix[1][2]+v1[2]*rotation_matrix[2][2];
}

void ChangeCoor(vector<float> &v1, const vector<vector<float> > &rotation_matrix, const vector<float> &translate_vec,vector<float> &v2)
{
	v2[0]=v1[0]*rotation_matrix[0][0]+v1[1]*rotation_matrix[1][0]+v1[2]*rotation_matrix[2][0]+translate_vec[0];
	v2[1]=v1[0]*rotation_matrix[0][1]+v1[1]*rotation_matrix[1][1]+v1[2]*rotation_matrix[2][1]+translate_vec[1];
	v2[2]=v1[0]*rotation_matrix[0][2]+v1[1]*rotation_matrix[1][2]+v1[2]*rotation_matrix[2][2]+translate_vec[2];
}

void RotateCoor(const vector<vector<float> > &target_stru,const vector<vector<float> > &refrence_stru, vector<vector<float> > &RotMatix, vector<float> &TranVect)
{
	int i, j, k, n;
	
	const int n_list=target_stru.size();
	
	//save two structure
	float ref_xlist[n_list][3],mov_xlist[n_list][3];
	for(i=0;i<n_list;i++) 
	{
		for(j=0;j<3;j++)
		{
			ref_xlist[i][j]=target_stru[i][j];
			mov_xlist[i][j]=refrence_stru[i][j];
		}
	}
	
	//translate two structure to the centre of mass, calculate correlation matrix
	float ref_com[3];
	float mov_com[3];
	float mov_to_ref[3];
      float R[3][3];
	
	/* calculate the centre of mass */
	for (i=0; i<3; i++)
	{ 
		mov_com[i] = 0.0;
		ref_com[i] = 0.0;
	}
	
	for (n=0; n<n_list; n++) 
	{
		for (i=0; i<3; i++)
		{ 
			mov_com[i] += mov_xlist[n][i];
			ref_com[i] += ref_xlist[n][i];
		}
	}
	
	for (i=0; i<3; i++)
	{
		mov_com[i] /= n_list;
		ref_com[i] /= n_list;
		mov_to_ref[i] = ref_com[i] - mov_com[i];
	}
	
	/* shift mov_xlist and ref_xlist to centre of mass */
	for (n=0; n<n_list; n++) 
	{
		for (i=0; i<3; i++)
		{ 
		  mov_xlist[n][i] -= mov_com[i];
		  ref_xlist[n][i] -= ref_com[i];
		}
	}
	
	/* initialize */
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++) R[i][j] = 0.0;
	}
	
	for (n=0; n<n_list; n++) 
	{
	/*
	 * correlation matrix R:   
	 *   R[i,j) = sum(over n): y(n,i) * x(n,j)  
	 *   where x(n) and y(n) are two vector sets   
	 */
		for (i=0; i<3; i++)
		{
		  	for (j=0; j<3; j++) R[i][j] += mov_xlist[n][i] * ref_xlist[n][j];
		}
	}
	
	//calculate rotation matrix and translate vector
	float Rt[3][3], RtR[3][3];
	float left_eigenvec[3][3], right_eigenvec[3][3], eigenval[3];
	float v[3];
	
	/* build Rt, transpose of R  */
	for (i=0; i<3; i++)
	{
		for (j=0; j<3; j++) Rt[i][j] = R[j][i];
	}
	
	/* make symmetric RtR = Rt X R */
	for (i=0; i<3; i++) 
	{
		for (j=0; j<3; j++)
		{
		  	RtR[i][j] = 0.0;
		  	for (k = 0; k<3; k++) RtR[i][j] += Rt[k][i] * R[j][k];
		}
	}
	
	diagonalize_symmetric(RtR, right_eigenvec, eigenval);
	
	
	/* right_eigenvec's should be an orthogonal system but could be left
	* or right-handed. Let's force into right-handed system.
	*/
	cross(&right_eigenvec[2][0], &right_eigenvec[0][0], &right_eigenvec[1][0]);
	
	/* From the Kabsch algorithm, the eigenvec's of RtR
	* are identical to the right_eigenvec's of R.
	* This means that left_eigenvec = R x right_eigenvec 
	*/
	for (i=0; i<3; i++) 
	{
		for (j=0; j<3; j++) left_eigenvec[i][j] = dot_prod(&right_eigenvec[i][0], &Rt[j][0]);
	}
	
	for (i=0; i<3; i++) normalize(&left_eigenvec[i][0]);
	
	/* 
	* Force left_eigenvec[2] to be orthogonal to the other vectors.
	* First check if the rotational matrices generated from the 
	* orthogonal eigenvectors are in a right-handed or left-handed
	* co-ordinate system - given by sigma. Sigma is needed to
	* resolve this ambiguity in calculating the RMSD.
	*/
	cross(v, &left_eigenvec[0][0], &left_eigenvec[1][0]);
	
	for (i=0; i<3; i++) left_eigenvec[2][i] = v[i]; 
	
	
	TranVect.assign(3, 0.);
	RotMatix.assign(3, TranVect);
	/* calc optimal rotation matrix that minimises residual */
	for (i=0;i<3; i++)
	{
		for (j=0; j<3; j++) 
		{
			for (k=0; k<3; k++) 
			      RotMatix[i][j] += left_eigenvec[k][i] * right_eigenvec[k][j];
		}
	}
	float temp[3];
	VectorTimesMatrix(ref_com, RotMatix, temp);
	TranVect[0]=refrence_stru[0][0]-mov_xlist[0][0]-temp[0];
	TranVect[1]=refrence_stru[0][1]-mov_xlist[0][1]-temp[1];
	TranVect[2]=refrence_stru[0][2]-mov_xlist[0][2]-temp[2];
	
}

