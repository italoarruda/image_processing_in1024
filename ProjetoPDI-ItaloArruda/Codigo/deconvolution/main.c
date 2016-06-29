#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "ppm.h"
#include "deconv_util.h"


void grad_gsm(int w, int h, float **x, float **org, float converge, int niter, double **fwd, double **bwd, float noise);
void reconstruct_from_gradient(int w, int h, int lvl,
                               double ***xx, double ***bb,
                               float **img,
                               float **grad_x, float **grad_y,
                               float **border,
                               float avg);

static int multigrid_level;
static double ***xx, ***bb;
static float **border_color;
static float border_avg;


static void finite_dif_x(int w, int h, float **img, float **dif);
static void finite_dif_y(int w, int h, float **img, float **dif);

static void set_border(int w, int h, float **img);


static void print_help()
{
    fprintf(stderr, "main in.pgm/ppm psf.pgm prior [options]\n");
    fprintf(stderr, "-e int  : edge taper size\n");
    fprintf(stderr, "-g      : undo gamma correction\n");
    fprintf(stderr, "-c float: convergence threshold\n");
    fprintf(stderr, "-n int  : maximum number of iterations\n");
    fprintf(stderr, "-t      : two-stage iteration\n");
    fprintf(stderr, "-r int  : post-process RL iterations\n");
    exit(0);
}

int	main(int argc, char *argv[])
{
    int i, j, k, n, w, h, ow, oh;
    int taper, gamma, grayscale, num_comp, niter_rl, niter_gsm, two_stage;
    float noise, converge, sum;
    unsigned char *input_image;
    float ***org, **img, **img_x, **img_y, **srs, **srs_x, **srs_y;
    double **fwd, **bwd;

    
    if(argc < 4) print_help();
    
    noise = (float)atof(argv[3]);
    noise = noise * noise;
    
    n = 4; taper = 0; gamma = 0; niter_rl = 0; niter_gsm = 0; two_stage = 0; converge = 0.01;
    while(n < argc)
    {
        if(strcmp(argv[n], "-g") == 0) { gamma = 1; n++; }
        else if(strcmp(argv[n], "-t") == 0) { two_stage = 1; n++; }
        else if(strcmp(argv[n], "-e") == 0 && n+1 < argc) { taper = atoi(argv[n+1]); n+=2; }
        else if(strcmp(argv[n], "-c") == 0 && n+1 < argc) { converge = atoff(argv[n+1]); n+=2; }
        else if(strcmp(argv[n], "-r") == 0 && n+1 < argc) { niter_rl = atoi(argv[n+1]); n+=2; }
        else if(strcmp(argv[n], "-n") == 0 && n+1 < argc) { niter_gsm = atoi(argv[n+1]); n+=2; }
        else print_help();
    }

    
    printf("loading image\n");
    {
        if(strcmp(&argv[1][strlen(argv[1])-4], ".pgm") == 0) grayscale = 1;
        else grayscale = 0;
    
        if(grayscale)
        {
            if((input_image = load_pgm(argv[1], &ow, &oh)) == NULL) exit(1);
            num_comp = 1;
        }
        else
        {
            if((input_image = load_ppm(argv[1], &ow, &oh)) == NULL) exit(1);
            num_comp = 3;
        }

        org = malloc_float3d(num_comp, ow, oh);
        for(k = 0; k < num_comp; k++)
            for(i = 0; i < ow; i++)
                for(j = 0; j < oh; j++) org[k][i][j] = (float)input_image[(j * ow + i) * num_comp + k] / 255.0;
        free(input_image);
    }
	printf("image size: %d x %d x %d\n", ow, oh, num_comp);
   
    if(gamma)
        for(i = 0; i < num_comp * ow * oh; i++) org[0][0][i] = inverse_gamma(org[0][0][i]);

    
    w = power2_no_less(ow);
    h = power2_no_less(oh);
    printf("enlarged to %d x %d\n", w, h);
    
    img = malloc_float2d(w, h);
    img_x = malloc_float2d(w, h);
    img_y = malloc_float2d(w, h);
    srs = malloc_float2d(w, h);
    srs_x = malloc_float2d(w, h);
    srs_y = malloc_float2d(w, h);

    
    load_kernel(argv[2], w, h, &fwd, &bwd);


    for(k = 0; k < num_comp; k++)
    {
        enlarge(ow, oh, org[k], w, h, img);
        simple_edge_taper(w, h, ow, oh, img, fwd, taper);

        set_border(w, h, img);

        for(i = 0; i < w * h; i++) srs[0][i] = img[0][i];

        printf("%d component X\n", k + 1);
        finite_dif_x(w, h, img, img_x);
        finite_dif_x(w, h, srs, srs_x);
        if(two_stage) grad_gsm(w, h, srs_x, img_x, converge, niter_gsm, fwd, bwd, noise * 2.0);
        grad_gsm(w, h, srs_x, img_x, converge, niter_gsm, fwd, bwd, noise);
        
        printf("%d component Y\n", k + 1);
        finite_dif_y(w, h, img, img_y);
        finite_dif_y(w, h, srs, srs_y);
        if(two_stage) grad_gsm(w, h, srs_y, img_y, converge, niter_gsm, fwd, bwd, noise * 2.0);
        grad_gsm(w, h, srs_y, img_y, converge, niter_gsm, fwd, bwd, noise);

        reconstruct_from_gradient(w, h, multigrid_level, xx, bb, srs, srs_x, srs_y, border_color, border_avg);

        sum = 0;
        for(i = 0; i < w * h; i++) sum += img[0][i] - srs[0][i];
        sum /= (float)(w * h);
        for(i = 0; i < w * h; i++)
        {
            srs[0][i] += sum;
            if(niter_rl > 0 && srs[0][i] < 1.0/255.0) srs[0][i] = 1.0/255.0;
        }

        if(niter_rl > 0) rich_lucy(w, h, srs, img, niter_rl, fwd, bwd);

        crop(w, h, srs, ow, oh, org[k]);
    }

    if(gamma)
        for(i = 0; i < num_comp * ow * oh; i++) org[0][0][i] = gamma_correct(org[0][0][i]);
    
    if(grayscale) store_pgm_normalized_grayscale("output.pgm", ow, oh, org[0]);
    else store_ppm_normalized_rgb("output.ppm", ow, oh, org);

    free_float3d(org);
    free_float2d(img);
    free_float2d(img_x);
    free_float2d(img_y);
    free_float2d(srs);
    free_float2d(srs_x);
    free_float2d(srs_y);
    free_double2d(fwd);
    free_double2d(bwd);
    
    return 0;
}


static void finite_dif_x(int w, int h, float **img, float **dif)
{
    int i, j;

    for(i = 0; i < w; i++)
        for(j = 0; j < h; j++)
            dif[i][j] = img[(i+1)%w][j] - img[i][j];
}
static void finite_dif_y(int w, int h, float **img, float **dif)
{
    int i, j;

    for(i = 0; i < w; i++)
        for(j = 0; j < h; j++)
            dif[i][j] = img[i][(j+1)%h] - img[i][j];
}


static void set_border(int w, int h, float **img)
{
    int i, j, k, tw, th;
    
	multigrid_level = 0;
    tw = w; th = h;
    while(tw > 1 && th > 1) { tw /= 2; th /= 2; multigrid_level++; }
    printf("min{%d, %d} = 2^%d\n", w, h, multigrid_level);
    multigrid_level++;

    xx = (double ***)malloc(multigrid_level * sizeof(double **));
    bb = (double ***)malloc(multigrid_level * sizeof(double **));
    
    tw = w;
	th = h;
    for(k = 0; k < multigrid_level; k++)
    {
        xx[k] = malloc_double2d(tw, th);
        bb[k] = malloc_double2d(tw, th);
		for(i = 0; i < tw; i++)
		for(j = 0; j < th; j++)
		{
			xx[k][i][j] = 0;
			bb[k][i][j] = 0;
		}
        tw = (tw + 1) / 2;
        th = (th + 1) / 2;
    }

    border_color = malloc_float2d(4, w > h ? w : h);

    border_avg = 0;
    for(j = 0; j < h; j++)
    {
        border_color[0][j] = img[0][j]; // left
        border_color[1][j] = img[w-1][j]; // right
        border_avg += img[0][j] + img[w-1][j];
    }
    for(i = 0; i < w; i++)
    {
        border_color[2][i] = img[i][0]; // top
        border_color[3][i] = img[i][h-1]; // bottom
        border_avg += img[i][0] + img[i][h-1];
    }
    border_avg /= (float)(2 * (w + h));

    // offset from average
    for(j = 0; j < h; j++)
    {
        border_color[0][j] -= border_avg;
        border_color[1][j] -= border_avg;
    }
    for(i = 0; i < w; i++)
    {
        border_color[2][i] -= border_avg;
        border_color[3][i] -= border_avg;
    }
}

