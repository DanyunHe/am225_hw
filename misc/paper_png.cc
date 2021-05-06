#include "write_png.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>

int main() {

    FILE *fp;
    float *rs;
    float *gs;
    float *bs;
    // Size of ouptut image
    const int m=1025,n=1025,mn=m*n;
    
    rs=(float*) malloc(sizeof(float)*mn);
    gs=(float*) malloc(sizeof(float)*mn);
    bs=(float*) malloc(sizeof(float)*mn);
    
    fp=fopen("../hw5/ftest.out/r.19","rb");
    fread(rs,sizeof(float),mn,fp);
    fclose(fp);

    fp=fopen("../hw5/ftest.out/g.19","rb");
    fread(gs,sizeof(float),mn,fp);
    fclose(fp);

    fp=fopen("../hw5/ftest.out/b.19","rb");
    fread(bs,sizeof(float),mn,fp);
    fclose(fp);

   
    // Allocate memory for the three color channels per gridpoint
    double *z=new double[3*mn],*zp=z,fr,t,x,y;

    // Loop over the pixels and set the colors according to a radial pattern
    int idx;
    for(int j=0;j<n;j++) {
        for(int i=0;i<m;i++) {
            
                idx=j+i*m;
                *(zp++)=rs[idx];
                // *(zp++)=1.;
                // *(zp++)=1.;
                *(zp++)=gs[idx];
                *(zp++)=bs[idx];

            
            // if(i==j) printf("r %g g %g b %g\n",rs[idx],gs[idx],bs[idx]);
        }
    }

    // Call routine to write PNG file, and free the dynamically allocated
    // memory
    write_png("paper19.png",m,n,z,0.,1.);
    delete [] z;
    delete [] rs;
    delete [] gs;
    delete [] bs;

}




