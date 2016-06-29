#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "ppm.h"
#include "deconv_util.h"

#define EPS 1.0e-3


int power2_no_less(int in)
{
    int out = 1;
    while(out < in) out *= 2;
    return out;
}


void load_kernel(char *filename, int w, int h, double ***fwd, double ***bwd)
{
    int i, j, x, y, kw, kh;
    float sum, max;
    unsigned char *input_image;
    int *ip;
    float **psf;
    double *work, **f, **b;
    
    fprintf(stderr, "loading kernel\n");
    if((input_image = load_pgm(filename, &kw, &kh)) == NULL) exit(1);
    fprintf(stderr, "kernel size: %d x %d\n", kw, kh);
    psf = malloc_float2d(kw, kh);
    max = 0;
    for(i = 0; i < kw; i++)
    for(j = 0; j < kh; j++)
    {
        psf[i][j] = (float)input_image[kw * j + i] / 255.0;
        if(max < psf[i][j]) max = psf[i][j];
    }
    free(input_image);
    sum = 0;
    for(i = 0; i < kw; i++)
    for(j = 0; j < kh; j++)
    {
        //if(psf[i][j] < max / 15.0) psf[i][j] = 0; // threshold
        //else sum += psf[i][j];
        sum += psf[i][j];
    }
    for(i = 0; i < kw; i++)
        for(j = 0; j < kh; j++) psf[i][j] /= sum;


    f = malloc_double2d(w, h * 2);
    b = malloc_double2d(w, h * 2);

    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        x = (i + kw/2) % w;
        y = (j + kh/2) % h;
        if(x < kw && y < kh) f[i][j*2] = (double)psf[x][y];
        else f[i][j*2] = 0;
        f[i][j*2+1] = 0;
    }

    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        x = (i + kw/2) % w;
        y = (j + kh/2) % h;
        if(x < kw && y < kh) b[i][j*2] = (double)psf[kw-1-x][kh-1-y];
        else b[i][j*2] = 0;
        b[i][j*2+1] = 0;
    }

    work = malloc_double1d((w > h ? w : h) / 2);
    ip = malloc_int1d(2 + (int)sqrt((w > h ? w : h) + 0.5));
    ip[0] = 0;
    cdft2d(w, h * 2, -1, f, NULL, ip, work);
    cdft2d(w, h * 2, -1, b, NULL, ip, work);

    free_float2d(psf);
    free_double1d(work);
    free_int1d(ip);
    *fwd = f;
    *bwd = b;
}


void enlarge(int ow, int oh, float **org, int w, int h, float **img)
{
    int i, j, x, y;
    
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        x = i - (w - ow) / 2;
        y = j - (h - oh) / 2;
        if(x < 0) x = 0; else if(x >= ow) x = ow - 1;
        if(y < 0) y = 0; else if(y >= oh) y = oh - 1;
        img[i][j] = org[x][y];
    }
}

void crop(int w, int h, float **img, int org_w, int org_h, float **org)
{
    int i, j, x, y;

    for(i = 0; i < org_w; i++)
    for(j = 0; j < org_h; j++)
    {
        x = i + (w - org_w) / 2;
        y = j + (h - org_h) / 2;
        org[i][j] = img[x][y];
    }
}

void simple_edge_taper(int w, int h, int org_w, int org_h, float **img, double **fwd, int taper_width)
{
    int i, j, x, y, min;
    float frac;
    double re, im;
    int *ip;
    double *work, **blr;

    if(taper_width <= 0) return;
    
    blr = malloc_double2d(w, h * 2);
    work = malloc_double1d((w > h ? w : h) / 2);
    ip = malloc_int1d(2 + (int)sqrt((w > h ? w : h) + 0.5));
    ip[0] = 0;
    
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        blr[i][j*2] = (double)img[i][j];
        blr[i][j*2+1] = 0;
    }
    cdft2d(w, h * 2, -1, blr, NULL, ip, work);
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        re = blr[i][j*2] * fwd[i][j*2] - blr[i][j*2+1] * fwd[i][j*2+1];
        im = blr[i][j*2] * fwd[i][j*2+1] + blr[i][j*2+1] * fwd[i][j*2];
        blr[i][j*2] = re / (double)(w * h);
        blr[i][j*2+1] = im / (double)(w * h);
    }
    cdft2d(w, h * 2, 1, blr, NULL, ip, work);
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        x = i - (w - org_w) / 2;
        y = j - (h - org_h) / 2;
        if(x < 0 || y < 0 || x >= org_w || y >= org_h) frac = 0;
        else
        {
            min = w + h;
            if(min > x) min = x;
            if(min > y) min = y;
            if(min > (org_w - x)) min = org_w - x;
            if(min > (org_h - y)) min = org_h - y;
            
            if(min > taper_width) frac = 1.0;
            else frac = (float)min / (float)taper_width;
        }
        img[i][j] = img[i][j] * frac + blr[i][j*2] * (1.0 - frac);
    }
    
    free_double2d(blr);
    free_double1d(work);
    free_int1d(ip);
}


void rich_lucy(int w, int h, float **x, float **org, int niter, double **fwd, double **bwd)
{
    int i, j, k;
    double re, im;
    int *ip;
    double *work, **blr;
    
    work = malloc_double1d((w > h ? w : h) / 2);
    ip = malloc_int1d(2 + (int)sqrt((w > h ? w : h) + 0.5));
    ip[0] = 0;

    blr = malloc_double2d(w, h * 2);

    for(k = 0; k < niter; k++)
    {
        for(i = 0; i < w; i++)
        for(j = 0; j < h; j++)
        {
            blr[i][j*2] = (double)x[i][j];
            blr[i][j*2+1] = 0;
        }
        cdft2d(w, h * 2, -1, blr, NULL, ip, work);
        for(i = 0; i < w; i++)
        for(j = 0; j < h; j++)
        {
            re = blr[i][j*2] * fwd[i][j*2] - blr[i][j*2+1] * fwd[i][j*2+1];
            im = blr[i][j*2] * fwd[i][j*2+1] + blr[i][j*2+1] * fwd[i][j*2];
            blr[i][j*2] = re / (double)(w * h);
            blr[i][j*2+1] = im / (double)(w * h);
        }
        cdft2d(w, h * 2, 1, blr, NULL, ip, work);

        for(i = 0; i < w; i++)
        for(j = 0; j < h; j++)
        {
            blr[i][j*2] = (double)org[i][j] / (blr[i][j*2] + EPS);
            blr[i][j*2+1] = 0;
        }
        cdft2d(w, h * 2, -1, blr, NULL, ip, work);
        for(i = 0; i < w; i++)
        for(j = 0; j < h; j++)
        {
            re = blr[i][j*2] * bwd[i][j*2] - blr[i][j*2+1] * bwd[i][j*2+1];
            im = blr[i][j*2] * bwd[i][j*2+1] + blr[i][j*2+1] * bwd[i][j*2];
            blr[i][j*2] = re / (double)(w * h);
            blr[i][j*2+1] = im / (double)(w * h);
        }
        cdft2d(w, h * 2, 1, blr, NULL, ip, work);
            
        for(i = 0; i < w; i++)
            for(j = 0; j < h; j++) x[i][j] *= (float)blr[i][j*2];
    }

    free_double1d(work);
    free_double2d(blr);
    free_int1d(ip);
}
