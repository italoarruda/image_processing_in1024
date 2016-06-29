#ifndef __DECONV_UTIL_H__
#define __DECONV_UTIL_H__

#ifdef __cplusplus
extern "C" {
#endif

void cdft2d(int, int, int, double **, double *, int *, double *);

int power2_no_less(int in);

void load_kernel(char *filename, int w, int h, double ***fwd, double ***bwd);

void enlarge(int ow, int oh, float **org, int w, int h, float **img);
void crop(int w, int h, float **img, int org_w, int org_h, float **org);
void simple_edge_taper(int w, int h, int org_w, int org_h, float **img, double **fwd, int taper_width);

void rich_lucy(int w, int h, float **x, float **org, int niter, double **fwd, double **bwd);

#ifdef __cplusplus
}
#endif

#endif
