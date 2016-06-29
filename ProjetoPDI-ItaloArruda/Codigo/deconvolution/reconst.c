#include <stdio.h>
#include <math.h>

#define EPS 1.0e-8
#define EPS_SCALE 1
#define BORDER 0.5

static void set_equation_coef(float **gx, float **gy, int w, int h, double **b, float **bdr);
static double residual(int w, int h, double **x, double **b);
static void v_cycle(int w, int h, int lvl, double ***x, double ***b);
static void multigrid(int w, int h, int lvl, double ***x, double ***b);

void reconstruct_from_gradient(int w, int h, int lvl,
                               double ***xx, double ***bb,
                               float **img,
                               float **grad_x, float **grad_y,
                               float **border,
                               float avg)
{
	int i, j;

    set_equation_coef(grad_x, grad_y, w, h, bb[0], border);

    multigrid(w, h, lvl, xx, bb);
    for(i = 0; i < w; i++)
        for(j = 0; j < h; j++) img[i][j] = (float)xx[0][i][j];// + avg;
}

static void set_equation_coef(float **gx, float **gy, int w, int h, double **b, float **bdr)
{
	int i, j;

    // Boundary condition for differenciation of grad(image) should be clamp-to-border-0
    // even when the boundary condition for the image is clamp-to-edge.
    //            image:  a       a   |   a       b       c      d      e
    //      grad(image):          0!  |  b-a     c-a     d-b    e-c
    // laplacian(image):              | c+a-2a  d+a-2b  e+a-2c
    for(i = 0; i < w; i++)
		for(j = 0; j < h; j++)
            //b[i][j] = 1.0 * (gx[i][j] - gx[(i+1)%w][j] + gy[i][j] - gy[i][(j+1)%h]);
			b[i][j] = 1.0 * (gx[(i+w-1)%w][j] - gx[i][j] + gy[i][(j+h-1)%h] - gy[i][j]);
            //b[i][j] = 1.0 * (gx[i<1?w-1:i-1][j] - gx[i][j] + gy[i][j<1?h-1:j-1] - gy[i][j]);

    for(j = 0; j < h; j++)
    {
//        b[0][j] += bdr[0][j];
//        b[w-1][j] += bdr[1][j];
    }
    for(i = 0; i < w; i++)
    {
//        b[i][0] += bdr[2][i];
//        b[i][h-1] += bdr[3][i];
    }
}

static double residual(int w, int h, double **x, double **b)
{
    int i, j;
    double sum = 0, dif;

    for(i = 0; i < w; i++)
    for(j = 0; j < h; j++)
    {
        dif = b[i][j] - 4.0 * x[i][j];
        dif += x[(i+w-1)%w][j];
        dif += x[(i+1)%w][j];
        dif += x[i][(j+h-1)%h];
        dif += x[i][(j+1)%h];
        sum += dif * dif;
    }

    return sqrt(sum) / (double)(w * h);
}

static void one_cycle(int w, int h, double **x, double **b)
{
    int i, j;
    double tmp;

    for(i = 0; i < w; i++) // red-black gauss-seidel
    for(j = (i % 2); j < h; j+=2)
    {
        tmp = b[i][j];
        tmp += x[(i+w-1)%w][j];
        tmp += x[(i+1)%w][j];
        tmp += x[i][(j+h-1)%h];
        tmp += x[i][(j+1)%h];
        x[i][j] = tmp * 0.25;
    }
    for(i = 0; i < w; i++)
    for(j = (1 - (i % 2)); j < h; j+=2)
    {
        tmp = b[i][j];
        tmp += x[(i+w-1)%w][j];
        tmp += x[(i+1)%w][j];
        tmp += x[i][(j+h-1)%h];
        tmp += x[i][(j+1)%h];
        x[i][j] = tmp * 0.25;
    }
}

static void v_cycle(int w, int h, int lvl, double ***x, double ***b)
{
    int i, j;
    double dif;
    
    one_cycle(w, h, x[lvl], b[lvl]);

    if(w <= 1 || h <= 1) return;

    for(i = 0; i < w-1; i+=2) // downsample residual
    for(j = 0; j < h-1; j+=2)
    {
		dif = b[lvl][i][j] + b[lvl][i][j+1] + b[lvl][i+1][j] + b[lvl][i+1][j+1]
			- 2.0 * (x[lvl][i][j] + x[lvl][i][j+1] + x[lvl][i+1][j] + x[lvl][i+1][j+1]);

        dif += x[lvl][(i+w-1)%w][j] + x[lvl][(i+w-1)%w][j+1];
        dif += x[lvl][i][(j+h-1)%h] + x[lvl][i+1][(j+h-1)%h];
        dif += x[lvl][(i+2)%w][j] + x[lvl][(i+2)%w][j+1];
        dif += x[lvl][i][(j+2)%h] + x[lvl][i+1][(j+2)%h];

        b[lvl+1][i/2][j/2] = sqrt(0.5) * dif;
    }

    for(i = 0; i < (w+1)/2; i++)
        for(j = 0; j < (h+1)/2; j++) x[lvl+1][i][j] = 0;
    
    v_cycle(w/2, h/2, lvl+1, x, b);
    
    for(i = 0; i < w; i++)
        for(j = 0; j < h; j++) x[lvl][i][j] += x[lvl+1][i/2][j/2]; // upsampling

    one_cycle(w, h, x[lvl], b[lvl]);
}

static double full_mg(int w, int h, int lvl, int scale, double ***x, double ***b)
{
    int i, j, k;
    double iter;
    double res = 1.0e+10;

    if(w <= 1 || h <= 1) return 0;

    iter = full_mg(w/2, h/2, lvl+1, EPS_SCALE * scale, x, b);
    
    for(i = 0; i < w; i++)
        for(j = 0; j < h; j++) x[lvl][i][j] = x[lvl+1][i/2][j/2]; // upsampling
    
    for(k = 0; res > EPS * (double)scale; k++)
    {
        v_cycle(w, h, lvl, x, b);
        res = residual(w, h, x[lvl], b[lvl]);
        //printf("%5.d: %f\n", k, res);
    }
    printf("multigrid %d: %d\n", lvl, k);

    return iter * 0.25 + (double)k;
}

static void multigrid(int w, int h, int lvl, double ***xx, double ***bb)
{
    int i, j, k, ww, hh;
    double iter;

    ww = w; hh = h;
    for(k = 0; k < lvl-1; k++)
    {
        for(i = 0; i < ww-1; i+=2)
        for(j = 0; j < hh-1; j+=2)
        {
            //xx[k+1][i/2][j/2] = (xx[k][i][j] + xx[k][i+1][j] + xx[k][i][j+1] + xx[k][i+1][j+1]);
            bb[k+1][i/2][j/2] = (bb[k][i][j] + bb[k][i+1][j] + bb[k][i][j+1] + bb[k][i+1][j+1]);
        }
        if(ww % 2)
        {
            for(j = 0; j < hh-1; j+=2)
            {
                //xx[k+1][ww/2][j/2] = xx[k][ww-1][j] + xx[k][ww-1][j+1];
                bb[k+1][ww/2][j/2] = bb[k][ww-1][j] + bb[k][ww-1][j+1];
            }
        }
        if(hh % 2)
        {
            for(i = 0; i < ww-1; i+=2)
            {
                //xx[k+1][i/2][hh/2] = xx[k][i][hh-1] + xx[k][i+1][hh-1];
                bb[k+1][i/2][hh/2] = bb[k][i][hh-1] + bb[k][i+1][hh-1];
            }
        }
        if((ww % 2) && (hh % 2))
        {
            //xx[k+1][ww/2][hh/2] = xx[k][ww-1][hh-1];
            bb[k+1][ww/2][hh/2] = bb[k][ww-1][hh-1];
        }
        ww /= 2; hh /= 2;
    }

    iter = full_mg(w, h, 0, 1, xx, bb);
	//double res = 1.0e+10;
    //for(k = 0; res > EPS; k++)
    //{
    //    v_cycle(w, h, 0, xx, bb);
    //    //res = residual(w, h, xx[0], bb[0]);
    //    //one_cycle_top(w, h, xx[0], bb[0]);
    //    res = residual_top(w, h, xx[0], bb[0]);
    //}
    //iter = (double)k;

    printf("multigrid total: %f\n", iter);
}
