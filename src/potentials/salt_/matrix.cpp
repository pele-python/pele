/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <math.h>
#include <stdlib.h>
#include "matrix.h"

#define PITIMES2 2*3.141592654

namespace votca { namespace tools {

void matrix::RandomRotation()
{
    matrix &M=(*this);
    double theta = drand48() * PITIMES2; /* Rotation about the pole (Z).      */
    double phi   = drand48() * PITIMES2; /* For direction of pole deflection. */
    double z     = drand48() * 2.0;      /* For magnitude of pole deflection. */

    /* Compute a vector V used for distributing points over the sphere  */
    /* via the reflection I - V Transpose(V).  This formulation of V    */
    /* will guarantee that if x[1] and x[2] are uniformly distributed,  */
    /* the reflected points will be uniform on the sphere.  Note that V */
    /* has length sqrt(2) to eliminate the 2 in the Householder matrix. */

    double r  = sqrt( z );
    double Vx = sin( phi ) * r;
    double Vy = cos( phi ) * r;
    double Vz = sqrt( 2.0 - z );    

    /* Compute the row vector S = Transpose(V) * R, where R is a simple */
    /* rotation by theta about the z-axis.  No need to compute Sz since */
    /* it's just Vz.                                                    */

    double st = sin( theta );
    double ct = cos( theta );
    double Sx = Vx * ct - Vy * st;
    double Sy = Vx * st + Vy * ct;

    /* Construct the rotation matrix  ( V Transpose(V) - I ) R, which   */
    /* is equivalent to V S - R.                                        */

    M[0][0] = Vx * Sx - ct;
    M[0][1] = Vx * Sy - st;
    M[0][2] = Vx * Vz;

    M[1][0] = Vy * Sx + st;
    M[1][1] = Vy * Sy - ct;
    M[1][2] = Vy * Vz;

    M[2][0] = Vz * Sx;
    M[2][1] = Vz * Sy;
    M[2][2] = 1.0 - z;   /* This equals Vz * Vz - 1.0 */   
}
        
int cjcbi(matrix &a, matrix &v, double eps, int jt)
{ 
    int n=3;
    int i,j,p,q,l;
    double fm,cn,sn,omega,x,y,d;
    l=1;
    
    v.UnitMatrix();
    
    while (true) {
        fm=0.0;
        for (i=1; i<=n-1; i++)
            for (j=0; j<=i-1; j++) {
                d=fabs(a[i][j]);
                if ((i!=j)&&(d>fm))
                    { fm=d; p=i; q=j;}
            }
        if (fm<eps)  return(1);
        if (l>jt)  return(-1);
        l=l+1;
        x=-a[p][q]; y=(a[q][q]-a[p][p])/2.0;
        omega=x/sqrt(x*x+y*y);
        if (y<0.0) omega=-omega;
        sn=1.0+sqrt(1.0-omega*omega);
        sn=omega/sqrt(2.0*sn);
        cn=sqrt(1.0-sn*sn);
        fm=a[p][p];
        a[p][p]=fm*cn*cn+a[q][q]*sn*sn+a[p][q]*omega;
        a[q][q]=fm*sn*sn+a[q][q]*cn*cn-a[p][q]*omega;
        a[p][q]=0.0; a[q][p]=0.0;
        for (j=0; j<=n-1; j++)
        if ((j!=p)&&(j!=q))
          { //u=p*n+j; w=q*n+j;
            fm=a[p][j];
            a[p][j]=fm*cn+a[q][j]*sn;
            a[q][j]=-fm*sn+a[q][j]*cn;
          }
        for (i=0; i<=n-1; i++)
          if ((i!=p)&&(i!=q))
            { //u=i*n+p; w=i*n+q;
              fm=a[i][p];
              a[i][p]=fm*cn+a[i][q]*sn;
              a[i][q]=-fm*sn+a[i][q]*cn;
            }
        for (i=0; i<=n-1; i++)
          { //u=i*n+p; w=i*n+q;
            fm=v[i][p];
            v[i][p]=fm*cn+v[i][q]*sn;
            v[i][q]=-fm*sn+v[i][q]*cn;
          }
      }
    return(1);
}
using namespace std;
void matrix::SolveEigensystem(eigensystem_t &out)
{
    matrix m(*this);
    matrix v;
    cjcbi(m, v);

    for(int i=0; i<3; ++i) {
        out.eigenvalues[i] = m[i][i];
        out.eigenvecs[i] = vec(v[0][i], v[1][i], v[2][i]);
        out.eigenvecs[i].normalize();
    }

    //cout << v << endl;
    // sort by eigenvalues
    if(out.eigenvalues[0] > out.eigenvalues[1]) {
        std::swap(out.eigenvalues[0], out.eigenvalues[1]);
        std::swap(out.eigenvecs[0], out.eigenvecs[1]);
    }
    if(out.eigenvalues[1] > out.eigenvalues[2]) {
        std::swap(out.eigenvalues[1], out.eigenvalues[2]);
        std::swap(out.eigenvecs[1], out.eigenvecs[2]);
    }
    if(out.eigenvalues[0] > out.eigenvalues[1]) {
        std::swap(out.eigenvalues[0], out.eigenvalues[1]);
        std::swap(out.eigenvecs[0], out.eigenvecs[1]);
    }    
}

void matrix::Invert()
{
    matrix mi;
    matrix &m=*this;
    mi[0][0] = m[1][1]*m[2][2] - m[1][2]*m[2][1];
    mi[1][0] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
    mi[2][0] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
    mi[0][1] = m[0][2]*m[2][1] - m[0][1]*m[2][2];
    mi[1][1] = m[0][0]*m[2][2] - m[0][2]*m[2][0];
    mi[2][1] = m[0][1]*m[2][0] - m[0][0]*m[2][1];
    mi[0][2] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
    mi[1][2] = m[0][2]*m[1][0] - m[0][0]*m[1][2];
    mi[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];

///calc the determinant, with the diagonal
    double D;

    D = m[0][0]*mi[0][0] + m[0][1]*mi[1][0] + m[0][2]*mi[2][0];
    
    mi /= D;
    *this = mi;
}

}}
