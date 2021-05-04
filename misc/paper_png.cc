#include "write_png.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>

int main() {

    FILE *fp;
    double *rs;
    double *gs;
    double *bs;
    // Size of ouptut image
    const int m=256,n=256,mn=m*n;
    
    rs=(double*) malloc(sizeof(double)*mn);
    gs=(double*) malloc(sizeof(double)*mn);
    bs=(double*) malloc(sizeof(double)*mn);
    
    fp=fopen("../hw5/ftest.out/r.1","rb");
    fread(rs,sizeof(double),mn,fp);
    fclose(fp);

    fp=fopen("../hw5/ftest.out/g.1","rb");
    fread(gs,sizeof(double),mn,fp);
    fclose(fp);

    fp=fopen("../hw5/ftest.out/b.1","rb");
    fread(bs,sizeof(double),mn,fp);
    fclose(fp);

   
    // Allocate memory for the three color channels per gridpoint
    double *z=new double[3*mn],*zp=z,fr,t,x,y;

    // Loop over the pixels and set the colors according to a radial pattern
    int idx;
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            idx=j+i*m;
            *(zp++)=rs[j+i*m];
            // *(zp++)=1.;
            *(zp++)=1.;
            *(zp++)=gs[j+i*m];
            // *(zp++)=bs[j+i*m];
            printf("r %g g %g b %g\n",rs[idx],gs[idx],bs[idx]);
        }
    }

    // Call routine to write PNG file, and free the dynamically allocated
    // memory
    write_png("test.png",m,n,z,0.,1.);
    delete [] z;
    delete [] rs;
    delete [] gs;
    delete [] bs;

}




