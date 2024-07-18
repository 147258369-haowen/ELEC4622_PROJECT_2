// Task_1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <iostream>
#include "io_bmp.h"
#include "filter.h"
#include "image_comps.h"
#include <math.h>
extern int height_offset;
extern float laplacianKernel[3][3];
int main(int argc, char* argv[]) {
    STATE state = CHECK_INPUT;
    int filterChoose = 0;
    bmp_in in;
    bmp_in in2;
    bmp_out out;
    my_image_comp* input_comps = nullptr;
    my_image_comp* input_comps2 = nullptr;
    my_image_comp* output_comps = nullptr;
    io_byte* line = nullptr;
    float** matrix = nullptr;
    my_image_comp* temp_comps = nullptr;
    my_image_comp* temp_comps_out = nullptr;
    my_image_comp* overlaid = nullptr;
    my_image_comp* nextinput = nullptr;
    my_image_comp** lapalacin = nullptr;
    my_image_comp** lapalacin2 = nullptr;
    my_image_comp** lapalacin3 = nullptr;
    my_image_comp* temp_image_2 = nullptr;
    my_image_comp* store_out_image_1 = new my_image_comp[3];
    my_image_comp* store_out_image_2 = new my_image_comp[3];
    ImageParam imageParam;
    float sigma = 0;
    float** lapkernel = allocateMatrix(3);
    Filter* sinc = nullptr;
    int newheight = 0;
    int tempheight = 0;
    bool change_picture = false;
    int first_width,first_height,first_initHeight;
    while (true) {
        switch (state) {
        case CHECK_INPUT:
            CheckInput(argc, argv, &sigma, &filterChoose, &imageParam);
            sinc = new Filter(imageParam.H, imageParam.D);
            sinc->SincKernelGenerate(sinc->sinc_buffer);
            imageParam.sincWidth = ((sinc->length) * 2 + 1);
            printf("[");
            for (int i = 0; i < imageParam.sincWidth; i++) {
                printf("%f,", sinc->sinc_buffer[i]);
            }
            printf("]\n");
            state = LOAD_PICTURE;
            break;
        case LOAD_PICTURE:
            //printf("%f\r\n", FilterInit(matrix, HIGHT, WIDTH));//归一化
            imageParam.GaussianDimension = GaussianWindowDimensionChoose(sigma);
            LoadImage(&in,&in2, &input_comps, &input_comps2 ,&output_comps, &line, &imageParam, &filterChoose, argv);
            state = LAPLACIAM_PYRAMID;
            break;
        case STICH_IMAGE: {
            printf("stich initheight:%d\r\n", imageParam.initheight);
            int D = Image_location(&overlaid, imageParam.initheight, store_out_image_1[0].width, &imageParam);
            imageParam.D = D;
            if (change_picture) {//first time enter
                lapalacin = allocate_laplacian(D, imageParam.origion_image_height, store_out_image_1[0].width, &imageParam);
                lapalacin2 = allocate_laplacian(D, imageParam.origion_image_height, store_out_image_1[0].width, &imageParam);
            }
            if (change_picture) Decompoment(store_out_image_1, lapalacin, D, &imageParam);//first time enter
            else Decompoment(store_out_image_2, lapalacin2, D, &imageParam);
            
            if (!change_picture) {
                int clear_flag = 1;
                Image_comps_init(&output_comps, &imageParam, store_out_image_1[0].height, store_out_image_1[0].width, 0);
                lapalacin3 = allocate_laplacian(D, imageParam.origion_image_height, store_out_image_1[0].width, &imageParam);
                for (int i = 0; i < D + 1; ++i) {
                    for (int n = 0; n < imageParam.num_comp; n++) {
                        lapalacin3[i][n].ImageStich(lapalacin[i] + n, lapalacin2[i] + n, &imageParam);
                    }       
                    Image_copy(&lapalacin3[i], &output_comps,&imageParam, clear_flag);
                    clear_flag = 0;
                }
                free_laplacian(lapalacin, D);
                free_laplacian(lapalacin2, D);
            }
            state = change_picture? STICH_IMAGE: INVERT_LAPLACIAN;
            change_picture = !change_picture;
            
        }break;
        case GAUSSIAN_FILTER:
            printf("Gaussan dimension: %d\r\n", imageParam.GaussianDimension);
            matrix = allocateMatrix(imageParam.GaussianDimension);
            LoadGaussianValue(matrix, sigma, imageParam.GaussianDimension);// load the gaussian value to the matrix
            for (int n = 0; n < imageParam.num_comp; n++)
                input_comps[n].perform_boundary_extension();
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, (imageParam.GaussianDimension - 1) / 2); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                horizontal(input_comps + n, temp_comps + n, matrix, imageParam.GaussianDimension, 0);
            for (int n = 0; n < imageParam.num_comp; n++)
                vertical(temp_comps + n, output_comps + n, matrix, imageParam.GaussianDimension, 0);
            state = (imageParam.gradientFlag) ? IMAGE_GRADIENT : OUTPUT_PICTURE;
            delete[] temp_comps;
            break;
        case INVERT_LAPLACIAN: {
            int D = Image_location(&overlaid, output_comps[0].height, output_comps[0].width, &imageParam);
            imageParam.D = D;
            lapalacin = allocate_laplacian(D, imageParam.origion_image_height, output_comps[0].width, &imageParam);
            Decompoment(output_comps, lapalacin, D, &imageParam);
            my_image_comp* input = lapalacin[D];
            my_image_comp* temp_image_1 = new  my_image_comp[imageParam.num_comp];
            for (int i = D; i >= 1; i--) {
                Image_upsample(&input, &temp_image_1, &imageParam);
                temp_image_2 = ImageRestore(temp_image_1, lapalacin[i - 1], &imageParam);
                printf("lapalacin[D - 1].height:%d\r\n", lapalacin[i - 1][0].height);
                input = temp_image_2;
            }
            imageParam.width = output_comps[0].width;
            printf("INVERT_LAPLACIAN width:%d\r\n", imageParam.width);
            
            /*           imageParam.initheight = lapalacin[D][0].height;
                       imageParam.initwidth = lapalacin[D][0].width;
                       imageParam.D = D;*/
            state = OUTPUT_PICTURE;

        }break;
        case LAPLACIAM_PYRAMID: {
            my_image_comp* temp_diff = new  my_image_comp[imageParam.num_comp];
            my_image_comp* temp_diffback = new  my_image_comp[imageParam.num_comp];
            first_width = imageParam.width;
            first_height = imageParam.height;
            
            Image_comps_init(&temp_diff, &imageParam, imageParam.height, imageParam.width, (imageParam.sincWidth - 1) / 2);
            temp_comps = new  my_image_comp[imageParam.num_comp];
            Image_comps_init(&temp_comps, &imageParam, imageParam.height, imageParam.width, (imageParam.sincWidth - 1) / 2);

            my_image_comp* temp_comps_inital = new  my_image_comp[imageParam.num_comp];
            Image_comps_init(&temp_comps_inital, &imageParam, imageParam.height, imageParam.width, (imageParam.sincWidth - 1) / 2);

            my_image_comp* temp_comps_3 = new  my_image_comp[imageParam.num_comp];
            Image_comps_init(&temp_comps_3, &imageParam, imageParam.height, imageParam.width, (imageParam.sincWidth - 1) / 2);

            my_image_comp* temp_comps_4 = new  my_image_comp[imageParam.num_comp];
            Image_comps_init(&temp_comps_4, &imageParam, imageParam.height, imageParam.width, (imageParam.sincWidth - 1) / 2);

            my_image_comp* temp_comps_5 = new  my_image_comp[imageParam.num_comp];
            Image_comps_init(&temp_comps_5, &imageParam, imageParam.height, imageParam.width, (imageParam.sincWidth - 1) / 2);

            overlaid = new  my_image_comp[imageParam.num_comp];
            first_initHeight = imageParam.initheight;
            tempheight = imageParam.initheight;//extension coloum
            for (int i = 0; i < imageParam.D; i++) {
                imageParam.initheight += (tempheight / 2);
                tempheight = (tempheight / 2);
            }
            my_image_comp* inputptr = nullptr;
            inputptr = change_picture ? input_comps2 : input_comps;
            printf("initheight:%d\r\n", imageParam.initheight);
            Image_comps_init(&overlaid, &imageParam, imageParam.initheight, imageParam.width, (imageParam.sincWidth - 1) / 2);
           
            Image_LPF(&inputptr, &temp_comps_inital, &sinc, &imageParam);//lpf 10
            Image_DownSample(&temp_comps_inital, &temp_comps, &imageParam);// 5
            Image_upsample(&temp_comps, &temp_comps_3, &imageParam);//10
            Laplacian_difference(&inputptr, &temp_comps_3, &imageParam);
            Image_copy(&temp_comps_3, &overlaid, &imageParam,0);
            //Image_DownSample(&temp_comps, &output_comps, &imageParam);//5

            nextinput = temp_comps;

            for (int i = 0; i < imageParam.D - 1; i++) {
                Image_LPF(&nextinput, &temp_comps_4, &sinc, &imageParam);//10
                Image_comps_init(&temp_comps_5, &imageParam, nextinput[0].height, nextinput[0].width, (imageParam.sincWidth - 1) / 2);
                Image_copy_no_offset(&nextinput, &temp_comps_5, &imageParam);
                Image_DownSample(&temp_comps_4, &output_comps, &imageParam);//5
                Image_upsample(&output_comps, &temp_diffback, &imageParam);//10
                Laplacian_difference(&temp_comps_5, &temp_diffback, &imageParam);//10
                Image_copy(&temp_diffback, &overlaid, &imageParam,0);
                nextinput = output_comps;
                //delete[] temp_diffback;
                /* Image_copy(&output_comps, &overlaid, &imageParam);*/
            }
            Image_copy(&output_comps, &overlaid, &imageParam,1);
            height_offset = 0;
            if (!change_picture) {
                Image_comps_init(&store_out_image_1, &imageParam, overlaid[0].height, overlaid[0].width, 0);
                Image_copy_no_offset(&overlaid, &store_out_image_1, &imageParam);
            }
            else {
                Image_comps_init(&store_out_image_2, &imageParam, overlaid[0].height, overlaid[0].width, 0);
                Image_copy_no_offset(&overlaid, &store_out_image_2, &imageParam);
            }
            delete[] overlaid;
            state = change_picture? STICH_IMAGE : LAPLACIAM_PYRAMID;
            if (!change_picture) {
                imageParam.initheight = first_initHeight;//1024
                imageParam.width = first_width;
                imageParam.height = first_height;
            }
            change_picture = true;
            delete[] temp_diff;
            //delete[] temp_diffback;
            delete[] temp_comps_inital;
            delete[] temp_comps_3;
            delete[] temp_comps_4;
            delete[] temp_comps_5;


           
            printf("restored width:%d,height:%d\r\n", imageParam.width, imageParam.height);
        }
                              break;
        case GAUSSIAN_PYRAMID:
            temp_comps = new  my_image_comp[imageParam.num_comp];
            overlaid = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, (imageParam.sincWidth - 1) / 2); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].perform_boundary_extension();
            //*************************************************
            tempheight = imageParam.initheight;//extension coloum
            for (int i = 0; i < imageParam.D; i++) {
                imageParam.initheight += (tempheight / 2);
                tempheight = (tempheight / 2);
            }
            printf("initheight:%d\r\n", imageParam.initheight);
            //*************************************************
            for (int n = 0; n < imageParam.num_comp; n++)
                overlaid[n].init(imageParam.initheight, imageParam.width, (imageParam.sincWidth - 1) / 2); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                overlaid[n].perform_boundary_extension();

            Image_LPF(&input_comps, &temp_comps, &sinc, &imageParam);
            Image_copy(&temp_comps, &overlaid, &imageParam,0);
            Image_DownSample(&temp_comps, &output_comps, &imageParam);
            nextinput = output_comps;
            Image_copy(&output_comps, &overlaid, &imageParam,0);
            for (int i = 0; i < imageParam.D - 1; i++) {
                Image_LPF(&nextinput, &temp_comps, &sinc, &imageParam);
                Image_DownSample(&temp_comps, &output_comps, &imageParam);
                nextinput = output_comps;
                Image_copy(&output_comps, &overlaid, &imageParam,0);
            }
            /* imageParam.width = imageParam.initwidth;
             imageParam.height = imageParam.initheight;*/
             // 48 24 12 6 3
            state = OUTPUT_PICTURE;
            delete[] temp_comps;
            break;
        case MOVING_AVERAGE_CHECK:
            VarianceLoopCheck(sigma, &imageParam);
            state = LOAD_PICTURE;
            break;
        case MOVING_AVERAGE:
            matrix = allocateMatrix(imageParam.MV_Dimension);
            MovingAverageSetValue(matrix, imageParam.MV_Dimension);
            for (int n = 0; n < imageParam.num_comp; n++)
                input_comps[n].perform_boundary_extension();
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, (imageParam.MV_Dimension - 1) / 2); // Don't need a border for output
            for (int n = 0; n < imageParam.num_comp; n++)
                horizontal(input_comps + n, temp_comps + n, matrix, imageParam.MV_Dimension, 1);
            for (int n = 0; n < imageParam.num_comp; n++)
                vertical(temp_comps + n, output_comps + n, matrix, imageParam.MV_Dimension, 1);
            state = (imageParam.gradientFlag) ? IMAGE_GRADIENT : OUTPUT_PICTURE;
            delete[] temp_comps;
            break;
        case IMAGE_GRADIENT:
            printf("Gradient\r\n");
#if LAPLACIAN == 1
            //for (int i = 0; i < 3; i++) {
            //    for (int j = 0; j < 3; j++) {
            //        lapkernel[i][j] = laplacianKernel[i][j];
            //    }
            //}
            //temp_comps = new  my_image_comp[imageParam.num_comp];
            //for (int n = 0; n < imageParam.num_comp; n++)
            //    temp_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            //for (int n = 0; n < imageParam.num_comp; n++)
            //    apply_filter_modified(output_comps+n, temp_comps +n, lapkernel,3);

            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            /************************************************************************************/
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].SecondGrradientHorizontalFilter(output_comps + n, 3, &imageParam, imageParam.alpha);
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].SecondGrradientverticalFilter(temp_comps + n, 3, &imageParam, imageParam.alpha);
            /************************************************************************************/
#else 
            temp_comps = new  my_image_comp[imageParam.num_comp];
            for (int n = 0; n < imageParam.num_comp; n++)
                temp_comps[n].init(imageParam.height, imageParam.width, 1); // Don't need a border for output
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < (temp_comps[i].stride * (temp_comps[i].height + 2 * temp_comps[i].border)); j++) {
                    temp_comps[i].handle[j] = output_comps[i].handle[j];
                }
            }
            for (int n = 0; n < imageParam.num_comp; n++)
                output_comps[n].GradientFilter(temp_comps + n, 3, imageParam.alpha);
#endif
            state = OUTPUT_PICTURE;
            break;
        case OUTPUT_PICTURE:
#if LAPLACIAN == 1
            //OutputImage(&out, input_comps, &temp_comps, &line, &imageParam, argv);
            OutputImage(&out, input_comps, &output_comps, &line, &imageParam, argv);
            bmp_in__close(&in);
            delete[] line;
            delete[] input_comps;
            delete[] temp_comps;

#else
            delete line;
            line = new io_byte[output_comps[0].width * 3];//* num_comps
            OutputImage(&out, input_comps, &temp_image_2, &line, &imageParam, argv);
            bmp_in__close(&in);

            delete[] line;
            delete[] input_comps;
            delete[] output_comps;
#endif


            return 1;
        default:
            break;
        }
    }
}