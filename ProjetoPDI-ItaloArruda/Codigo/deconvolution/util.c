#include <stdio.h>
#include <stdlib.h>
#include "malloc/malloc.h"
#include <math.h>
#include "util.h"


#define MyMalloc(type)                                                       \
type   *malloc_ ## type ## 1d(int n)                                         \
{                                                                            \
    type *buf;                                                               \
                                                                             \
    if((buf = (type *)malloc(sizeof(type) * n)) == NULL)                     \
    {                                                                        \
        fprintf(stderr, "failed to allocate memory\n"); exit(1);             \
    }                                                                        \
    return buf;                                                              \
}                                                                            \
                                                                             \
type  **malloc_ ## type ## 2d(int w, int h)                                  \
{                                                                            \
    type **buf;                                                              \
    int i;                                                                   \
                                                                             \
    if((buf = (type **)malloc(sizeof(type *) * w)) == NULL)                  \
    {                                                                        \
        fprintf(stderr, "failed to allocate memory\n"); exit(1);             \
    }                                                                        \
    if((buf[0] = (type *)malloc(sizeof(type) * w * h)) == NULL)              \
    {                                                                        \
        fprintf(stderr, "failed to allocate memory\n"); exit(1);             \
    }                                                                        \
    for(i = 1; i < w; i++) buf[i] = buf[0] + i * h;                          \
    return buf;                                                              \
}                                                                            \
                                                                             \
type ***malloc_ ## type ## 3d(int w, int h, int d)                           \
{                                                                            \
    type ***buf;                                                             \
    int i;                                                                   \
                                                                             \
    if((buf = (type ***)malloc(sizeof(type **) * w)) == NULL)                \
    {                                                                        \
        fprintf(stderr, "failed to allocate memory\n"); exit(1);             \
    }                                                                        \
    if((buf[0] = (type **)malloc(sizeof(type *) * w * h)) == NULL)           \
    {                                                                        \
        fprintf(stderr, "failed to allocate memory\n"); exit(1);             \
    }                                                                        \
    for(i = 1; i < w; i++) buf[i] = buf[0] + i * h;                          \
    if((buf[0][0] = (type *)malloc(sizeof(type) * w * h * d)) == NULL)       \
    {                                                                        \
        fprintf(stderr, "failed to allocate memory\n"); exit(1);             \
    }                                                                        \
    for(i = 1; i < w * h; i++) buf[0][i] = buf[0][0] + i * d;                \
    return buf;                                                              \
}                                                                            \
                                                                             \
void free_ ## type ## 1d(type   *buf)                                        \
{                                                                            \
    free(buf);                                                               \
}                                                                            \
                                                                             \
void free_ ## type ## 2d(type  **buf)                                        \
{                                                                            \
    free(buf[0]);                                                            \
    free(buf);                                                               \
}                                                                            \
                                                                             \
void free_ ## type ## 3d(type ***buf)                                        \
{                                                                            \
    free(buf[0][0]);                                                         \
    free(buf[0]);                                                            \
    free(buf);                                                               \
}

MyMalloc(uchar)
MyMalloc(int)
MyMalloc(float)
MyMalloc(double)
#undef MyMalloc


static const float gamma_thres = 0.0812;
static const float gamma_exp = 0.45;
static const float gamma_scale = 4.5;
static const float gamma_ofs = 0.099;

float gamma_correct(float color)
{
    if(color < gamma_thres / gamma_scale) return color * gamma_scale;
    else return (1.0 + gamma_ofs) * (float)pow(color, gamma_exp) - gamma_ofs;
}

float inverse_gamma(float color)
{
    if(color < gamma_thres) return color / gamma_scale;
    else return (float)pow((color + gamma_ofs) / (1.0 + gamma_ofs), 1.0 / gamma_exp);
}


