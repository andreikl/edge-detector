#include <math.h>       // hypot
// stl
#include <ctime>        // clock
#include <cstring>      // memcpy

#include <utils.hpp>

using namespace std;

extern bool is_debug;
extern int debug_x;
extern int debug_y;
extern int width;
extern int height;
extern int frames;

extern string derivatives_x_output;
extern string derivatives_xx_output;
extern string derivatives_y_output;
extern string derivatives_yy_output;
extern string derivatives_xy_output;
extern string suppression_output;
extern string gaussian_output;
extern string output;

void example_harris(double hk, int hr) {
    // converts input parameters to variables
    pixel_t image[width * height];
    pixel_t image_gaussian[width * height];
    pixel_t out[width * height];            

    // initialise output stream
    ostream& os = create_output_stream(output);
   
    for (int k = 0; k < frames; k++) {
        //reads one frame to array
        cin.read((char*)(image), width * height * sizeof(pixel_t));

        if (k == frames - 1) {
            const clock_t calc_begin_time = clock();

            auto gfilter = gaussian_filter<pixel_t, pixel_t>(gaussian_output, 0.8, 5);
            memcpy(image_gaussian, image, width * height * sizeof(pixel_t));
            gfilter.do_process(image, image_gaussian, width, height);

            // borders are -0.5 0.5: -128 to 128, -2, 2: -512 to 512
            const double *dx = new double[3] {-0.5, 0, 0.5};
            derivatives_t image_dx[width * height];
            memset(image_dx, 0, width * height * sizeof(derivatives_t));
            auto dxfilter = convolution_filter(derivatives_x_output, dx, 3, 1);
            dxfilter.do_process(image_gaussian, image_dx, width, height);

            const double *dy = new double[3] {0.5, 0, -0.5};
            derivatives_t image_dy[width * height];
            memset(image_dy, 0, width * height * sizeof(derivatives_t));
            auto dyfilter = convolution_filter(derivatives_y_output, dy, 1, 3);
            dyfilter.do_process(image_gaussian, image_dy, width, height);
            
            derivatives_t image_dxx[width * height];
            derivatives_t image_dyy[width * height];
            derivatives_t image_dxy[width * height];
            //#pragma omp parallel for
            for (int i = 0; i < width * height; i++) {
                image_dxx[i] = image_dx[i] * image_dx[i];
                image_dyy[i] = image_dy[i] * image_dy[i];
                image_dxy[i] = image_dx[i] * image_dy[i];
            }

            auto gfilterdxx = gaussian_filter<derivatives_t, derivatives_t>(derivatives_xx_output, 1.0, 5, is_debug, debug_x, debug_y);
            unique_ptr<derivatives_t> image_gaussian_dxx(new derivatives_t[width * height]);
            derivatives_t* image_gaussian_dxx_ptr = image_gaussian_dxx.get();
            memset(image_gaussian_dxx_ptr, 0, width * height * sizeof(derivatives_t));
            gfilterdxx.do_process(image_dxx, image_gaussian_dxx_ptr, width, height);

            auto gfilterdyy = gaussian_filter<derivatives_t, derivatives_t>(derivatives_yy_output, 1.0, 5, is_debug, debug_x, debug_y);
            unique_ptr<derivatives_t> image_gaussian_dyy(new derivatives_t[width * height]);
            derivatives_t* image_gaussian_dyy_ptr = image_gaussian_dyy.get();
            memset(image_gaussian_dyy_ptr, 0, width * height * sizeof(derivatives_t));
            gfilterdyy.do_process(image_dyy, image_gaussian_dyy_ptr, width, height);

            auto gfilterdxy = gaussian_filter<derivatives_t, derivatives_t>(derivatives_xy_output, 1.0, 5, is_debug, debug_x, debug_y);
            unique_ptr<derivatives_t> image_gaussian_dxy(new derivatives_t[width * height]);
            derivatives_t* image_gaussian_dxy_ptr = image_gaussian_dxy.get();
            memset(image_gaussian_dxy_ptr, 0, width * height * sizeof(derivatives_t));
            gfilterdxy.do_process(image_dxy, image_gaussian_dxy_ptr, width, height);

            //unique_ptr<derivatives_t> R(new derivatives_t[width * height]);
            
            
            
            for (int i = 0; i < width * height; i++) {
                //iy2 * ix2 - pow(ixiy, 2)
                int det = (int)image_gaussian_dxx_ptr[i] * (int)image_gaussian_dyy_ptr[i] - (int)image_gaussian_dxy_ptr[i] * (int)image_gaussian_dxy_ptr[i];
                int trace = image_gaussian_dxx_ptr[i] + image_gaussian_dyy_ptr[i];
                int R = det - hk * trace * trace;
                if (is_debug && i % width == debug_x && i / width == debug_y) {
                    int test = image_gaussian_dxy_ptr[i];
                    cerr << "i: " << i << ", test: " << test << endl;
                    cerr << "image_gaussian_dxx: " << image_gaussian_dxx_ptr[i] << ", image_gaussian_dyy: " << image_gaussian_dyy_ptr[i] << ", image_gaussian_dxy: " << image_gaussian_dxy_ptr[i] << endl;
                    cerr << "R: " << R << ", det: " << det << ", trace: " << hk * trace * trace << endl;
                }
                out[i] = (R > hr)? MAX_BRIGHTNESS: (unsigned char)0;
            }

            //non-maximum suppression
            //derivatives_t nms[width * height];            
            //#pragma omp parallel for
            /*for (int i = 0; i < width * height; i++) {
                const int nn = i - width;
                const int ss = i + width;
                const int ww = i + 1;
                const int ee = i - 1;
                const int nw = nn + 1;
                const int ne = nn - 1;
                const int sw = ss + 1;
                const int se = ss - 1;

                const double dir = (fmod(atan2(image_dy[i], image_dx[i]) + M_PI, M_PI) / M_PI) * 8;
                if (((dir <= 1 || dir > 7) && image_dxy[i] > image_dxy[ee] && image_dxy[i] > image_dxy[ww]) || // 0 deg
                    ((dir > 1 && dir <= 3) && image_dxy[i] > image_dxy[nw] && image_dxy[i] > image_dxy[se]) || // 45 deg
                    ((dir > 3 && dir <= 5) && image_dxy[i] > image_dxy[nn] && image_dxy[i] > image_dxy[ss]) || // 90 deg
                    ((dir > 5 && dir <= 7) && image_dxy[i] > image_dxy[ne] && image_dxy[i] > image_dxy[sw])) {  // 135 deg
                    nms[i] = image_dxy[i];
                }
                else {
                    nms[i] = 0;
                }                
            }*/

            if (is_debug) {
                cerr << "Calculation time: " << float(clock () - calc_begin_time) /  CLOCKS_PER_SEC << endl;
            }
        }
    }

    const clock_t write_begin_time = clock();

    //writes ppm image
    os << "P6\n" << width << " " << height << "\n255\n";
    for (int i = width * height - 1; i >= 0; i--) {
        if (out[i] == MAX_BRIGHTNESS) {
            os << (unsigned char)MAX_BRIGHTNESS;
            os << (unsigned char)MAX_BRIGHTNESS;
            os << (unsigned char)0;
        } else {
            os << image[i];
            os << image[i];
            os << image[i];
        }
    }

    if (is_debug) {
        cerr << "Writing time: " << float(clock () - write_begin_time) /  CLOCKS_PER_SEC << endl;
    }

    delete_output_stream(os);
}
