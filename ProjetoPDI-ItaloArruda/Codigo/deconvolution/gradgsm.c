#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "util.h"
#include "deconv_util.h"

#define NUM_INNER_LOOP 10
#define EPS            (1.0 / (255.0 * 255.0))


static float beta;


#define SCALE 0.085
#define POWER 0.3

static double heavy_tailed_prior(float theta)
{
    // Jeffreys
//    return 1.0 / sqrt(theta * theta * 0.5 + EPS);
    // Generalized Gaussian
    return -pow(theta * theta + EPS, POWER / 2.0) / SCALE;
}

static float expect_heavy_tailed_prior(float theta)
{
    // Jeffreys
//    return 1.0 / (theta * theta * 0.5 + EPS);
    // Garrote
//    float t2 = theta * theta;
//    return 6.0 / (t2 + (float)(fabs(theta) * sqrt(t2 + 12.0 * noise)));
    // Laplace
//    return beta / (float)(sqrt(noise) * fabs(theta));
//    return beta / (float)(sqrt(noise * (theta * theta + 0.1)));
//    return 1.0 / (theta * theta);
    // Generalized Gaussian
    return POWER / SCALE / (float)pow(theta * theta + EPS, (2.0 - POWER) / 2.0);
}

static float compute_log_posterior(int w, int h, float **x, float **org, double **fwd, double **blr, int *ip, double *work, float noise);


void grad_gsm(int w, int h, float **x, float **org, float converge, int niter, double **fwd, double **bwd, float noise)
{
    int i, j, m, n, first;
    float lambda1, lambdan;
    float alpha, beta0, rho_opt;
    float posterior_prev, posterior;
    double re, im;
    int *ip;
    float **y, **d, **f, **x_prev, **x_next, **cr;
    double *work, **blr, **scl;

    y = malloc_float2d(w, h);
    f = malloc_float2d(w, h);
    d = malloc_float2d(w, h);
    x_prev = malloc_float2d(w, h);
    x_next = malloc_float2d(w, h);
    cr = malloc_float2d(w, h);
    blr = malloc_double2d(w, h * 2);
    scl = malloc_double2d(w, h);
    work = malloc_double1d((w > h ? w : h) / 2);
    ip = malloc_int1d(2 + (int)sqrt((w > h ? w : h) + 0.5));
    ip[0] = 0;

    
    lambdan = 1.0;
    lambda1 = 0;
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        re = fwd[i][j*2] * bwd[i][j*2] - fwd[i][j*2+1] * bwd[i][j*2+1];
        //im = fwd[i][j*2] * bwd[i][j*2+1] + fwd[i][j*2+1] * bwd[i][j*2]; // this should be zero
        scl[i][j] = re / (double)(w * h);
        lambda1 += (float)(re * re);
    }
    lambda1 /= (float)(w * h);
    //printf("lambda1: %f\n", lambda1);
    rho_opt = (1.0 - (float)sqrt(lambda1 / lambdan)) / (1.0 + (float)sqrt(lambda1 / lambdan));
    alpha = rho_opt * rho_opt + 1.0;
    beta = 2.0 * alpha / (lambda1 + lambdan);
    beta0 = beta / alpha;
    

    // prepare y' = H^T y
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        blr[i][j*2] = (double)org[i][j]; // observation
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
        for(j = 0; j < h; j++) y[i][j] = (float)blr[i][j*2];

    
    posterior = compute_log_posterior(w, h, x, org, fwd, blr, ip, work, noise);
    
    n = 0;
    do
    {
        posterior_prev = posterior;
        
        for(i = 0; i < w; i++)
            for(j = 0; j < h; j++) d[i][j] = expect_heavy_tailed_prior(x[i][j]);

        first = 1;
        
        for(m = 0; m < NUM_INNER_LOOP; m++)
        {
            //printf("outer %d, inner %d\n", n, m);

            for(i = 0; i < w; i++)
                for(j = 0; j < h; j++) f[i][j] = x[i][j];
            
            for(i = 0; i < w; i++)
            for(j = 0; j < h; j++)
            {
                blr[i][j*2] = f[i][j];
                blr[i][j*2+1] = 0;
            }
            cdft2d(w, h * 2, -1, blr, NULL, ip, work);
            for(i = 0; i < w; i++)
            for(j = 0; j < h; j++)
            {
                blr[i][j*2] *= scl[i][j];//fwd[i][j*2] / (double)(w * h);
                blr[i][j*2+1] *= scl[i][j];//fwd[i][j*2] / (double)(w * h);
            }
            cdft2d(w, h * 2, 1, blr, NULL, ip, work);
            for(i = 0; i < w; i++)
                for(j = 0; j < h; j++) f[i][j] = (float)blr[i][j*2];
            
            for(i = 0; i < w; i++)
                for(j = 0; j < h; j++) cr[i][j] = x[i][j] + (f[i][j] - x[i][j] - y[i][j]) / (1.0 + noise * d[i][j]);

            if(first)
            {
                for(i = 0; i < w; i++)
                    for(j = 0; j < h; j++) x_next[i][j] = x[i][j] - beta0 * cr[i][j];
                first = 0;
            }
            else
            {
                for(i = 0; i < w; i++)
                    for(j = 0; j < h; j++) x_next[i][j] = alpha * x[i][j] + (1.0 - alpha) * x_prev[i][j] - beta * cr[i][j];
            }

            for(i = 0; i < w; i++)
                for(j = 0; j < h; j++) x_prev[i][j] = x[i][j];
            
            for(i = 0; i < w; i++)
                for(j = 0; j < h; j++) x[i][j] = x_next[i][j];
        }

        posterior = compute_log_posterior(w, h, x, org, fwd, blr, ip, work, noise);
        //printf("%f, %f, %f\n", posterior, posterior_prev, (posterior - posterior_prev) / -posterior);

        n++; printf("%d iterations\n", n);
    }
    while(niter ? n < niter : posterior - posterior_prev > converge * fabs(posterior));

    free_float2d(y);
    free_float2d(f);
    free_float2d(d);
    free_float2d(x_prev);
    free_float2d(x_next);
    free_float2d(cr);
    free_double2d(blr);
    free_double2d(scl);
    free_double1d(work);
    free_int1d(ip);
}


static float compute_log_posterior(int w, int h, float **x, float **org, double **fwd, double **blr, int *ip, double *work, float noise)
{
    int i, j;
    double re, im, p = 0;

    for(i = 0; i < w; i++)
        for(j = 0; j < h; j++) p += heavy_tailed_prior(x[i][j]);
    p *= 2.0 * noise;
    
    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        blr[i][j*2] = x[i][j];
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
        for(j = 0; j < h; j++) p -= (org[i][j] - blr[i][j*2]) * (org[i][j] - blr[i][j*2]);

    return (float)p;
}
