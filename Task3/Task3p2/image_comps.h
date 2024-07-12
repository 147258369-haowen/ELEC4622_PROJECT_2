/*****************************************************************************/
// File: image_comps.h
// Author: David Taubman
// Last Revised: 13 August, 2007
/*****************************************************************************/
// Copyright 2007, David Taubman, The University of New South Wales (UNSW)
/*****************************************************************************/

/*****************************************************************************/
/* STRUCT                        my_image_comp                               */
/*****************************************************************************/
#ifndef IMAGE_COMPS_H
#define IMAGE_COMPS_H
#include <stdio.h>
#include "filter.h"
#define WIDTH 5
#define HIGHT 5
#define PI 3.14159265
class Filter;
//#define DEBUG
#define CLAMP_TO_BYTE(sum) ((sum) = (sum) > 255 ? 255 : ((sum) < 0 ? 0 : (sum)))
#define LAPLACIAN 0
typedef enum {
    CHECK_INPUT,
    LOAD_PICTURE,
    GAUSSIAN_FILTER,
    MOVING_AVERAGE,
    MOVING_AVERAGE_CHECK,
    OUTPUT_PICTURE,
    IMAGE_GRADIENT,
    LAPLACIAN_IMAGE,
    SINC_FILTER,
    LAPLACIAM_PYRAMID,
    INVERT_LAPLACIAN,
}STATE;
typedef struct {
    int width;
    int height;
    int num_comp;
    io_byte* line;
    int MV_Dimension;
    int GaussianDimension;
    bool gradientFlag;
    float alpha;
    int sincWidth;
    int D;
    int H;
    int initwidth;
    int initheight;
    int origion_image_height;
}ImageParam;
#define GAUSSIAN 1
#define MOVINGAVERAGE 0
struct my_image_comp {
    // Data members: (these occupy space in the structure's block of memory)
    int width;
    int height;
    int stride;
    int border; // Extra rows/cols to leave around the boundary
    float* handle; // Points to start of allocated memory buffer
    float* buf; // Points to the first real image sample
    int height_offset;
    // Function members: (these do not occupy any space in memory)
    my_image_comp()
    {
        width = height = stride = border = 0;  handle = buf = NULL;
        height_offset = 0;
    }
    ~my_image_comp()
    {
        if (handle != NULL) delete[] handle;
    }
    void init(int height, int width, int border)
    {
        this->width = width;  this->height = height;  this->border = border;
        stride = width + 2 * border;//扩展后的长宽
        if (handle != NULL)
            delete[] handle; // Delete mem allocated by any previous `init' call
        handle = new float[stride * (height + 2 * border)];
        buf = handle + (border * stride) + border;//the real picture start point
    }
    void perform_boundary_extension();
    void apply_filter_modified_simo(my_image_comp* in, my_image_comp* out, float** inputfilter, int width);
    void vector_filter(my_image_comp* in, int dimension);
    void vector_horizontal_filter(my_image_comp* in, int dimension);
    void GrradientHorizontalFilter(my_image_comp* in, int dimension, int alpha);
    void GrradientverticalFilter(my_image_comp* in, int width, int alpha);
    void GradientFilter(my_image_comp* in, int width, int alpha);
    void SecondGrradientHorizontalFilter(my_image_comp* in, int dimension, ImageParam* imagP, int alpha);
    void SecondGrradientverticalFilter(my_image_comp* in, int dimension, ImageParam* imagP, int alpha);
    // This function is implemented in "filtering_main.cpp".
};
void Image_LPF(my_image_comp** input_comps, my_image_comp** output_comps, Filter** sinc, ImageParam* imageParam);
void apply_filter(my_image_comp* in, my_image_comp* out);
float FilterNormalized(float** input, int dimension);
void apply_filter_modified(my_image_comp* in, my_image_comp* out, float** inputfilter, int width);
void apply_filter_modified_simo(my_image_comp* in, my_image_comp* out, float** inputfilter, int width);
void unsharp_mask_filter(float** inputfilter, int width, float alpha);
void CheckInput(int argc, char* argv[], float* sigma, int* filterChooseFlag, ImageParam* param);
int OutputImage(bmp_out* out, my_image_comp* input_comps, my_image_comp** output_comps, io_byte** line, ImageParam* imageParam, char** argv);
float GaussianFillKernel(int x, int y, float sigma);
int LoadGaussianValue(float** matrix, float sigma, int dimension);
float** allocateMatrix(int dimension);
void FreeMatrix(float** matrix, int width);
int VarianceLoopCheck(float sigma, ImageParam* imageParam);
void LoadImage(bmp_in* in, my_image_comp** input_comps, my_image_comp** output_comps, io_byte** line,
    ImageParam* imageParam, int* filterChoose, char** argv);
void MovingAverageSetValue(float** matrix, int dimension);
int GaussianWindowDimensionChoose(float sigma);
void listshift(float* buffer, float** ptr, int width);
void horizontal(my_image_comp* in, my_image_comp* out, float** inputfilter, int width, int G_MF_flag);
void vertical(my_image_comp* in, my_image_comp* out, float** inputfilter, int width, int G_MF_flag);
void Image_DownSample(my_image_comp** input_comps, my_image_comp** output_comps, ImageParam* imageParam);
void Image_copy(my_image_comp** input_comps, my_image_comp** output_comps, ImageParam* imageParam);
void Image_upsample(my_image_comp** input_comps, my_image_comp** output_comps, ImageParam* imageParam);
void Laplacian_difference(my_image_comp** input_comps, my_image_comp** output_comps, ImageParam* imageParam);
void Image_comps_init(my_image_comp** temp_comps, ImageParam* imageParam, int height, int width, int extention);
int Image_location(my_image_comp** temp_comps, int Height, int width, ImageParam* imageParam);
void Decompoment(my_image_comp* in, my_image_comp** out, int D);
my_image_comp** allocate_laplacian(int D, int height, int width);
my_image_comp* ImageRestore(my_image_comp* image_upsample, my_image_comp* image_laplacian, ImageParam* imageParam);
#endif