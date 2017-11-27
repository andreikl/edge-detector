#include <math.h>       // hypot
// stl
#include <ctime>        // clock
#include <cstring>      // memcpy

#include <utils.hpp>

using namespace std;

extern bool is_debug;
extern int width;
extern int height;
extern int frames;

extern string derivatives_x_output;
extern string derivatives_y_output;
extern string derivatives_xy_output;
extern string suppression_output;
extern string output;

void example_harris() {
    // converts input parameters to variables
    pixel_t image[width * height];
    pixel_t out[width * height];            

    // initialise output stream
    ostream& os = create_output_stream(output);
   
    for (int k = 0; k < frames; k++) {
        //reads one frame to array
        cin.read((char*)(image), width * height * sizeof(pixel_t));

        if (k == frames - 1) {
            //exampleContrast(image, width, height);

            const clock_t calc_begin_time = clock();

            // borders are -0.5 0.5: -128 to 128, -2, 2: -512 to 512
            const double dx[] = {-2, 2};
            derivatives_t image_dx[width * height];
            memset(image_dx, 0, width * height * sizeof(derivatives_t));
            auto dxfilter = convolution_filter(derivatives_x_output, dx, 2, 1);
            dxfilter.do_process(image, image_dx, width, height);

            const double dy[] = {2, -2};
            derivatives_t image_dy[width * height];
            memset(image_dy, 0, width * height * sizeof(derivatives_t));
            auto dyfilter = convolution_filter(derivatives_y_output, dx, 1, 2);
            dyfilter.do_process(image, image_dy, width, height);
            
            derivatives_t image_dxy[width * height];
            memset(image_dxy, 0, width * height * sizeof(derivatives_t));

            //#pragma omp parallel for
            for (int i = 0; i < width * height; i++) {
                image_dxy[i] = (derivatives_t)hypot(image_dx[i], image_dy[i]);
            }

            const clock_t dxy_write_begin_time = clock();
            ostream& dxyos = create_output_stream(derivatives_xy_output);
            //writes ppm image
            dxyos << "P6\n" << width << " " << height << "\n255\n";
            for (int i = width * height - 1; i >= 0; i--) {
                dxyos << (unsigned char)image_dxy[i];
                dxyos << (unsigned char)image_dxy[i];
                dxyos << (unsigned char)image_dxy[i];
            }

            if (is_debug) {
                cerr << "Derivative XY filter writing time: " << float(clock () - dxy_write_begin_time) /  CLOCKS_PER_SEC << endl;
            }

            delete_output_stream(dxyos);


            //non-maximum suppression
            derivatives_t nms[width * height];            
            //#pragma omp parallel for
            for (int i = 0; i < width * height; i++) {
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
            }

            const clock_t nms_write_begin_time = clock();
            ostream& nmsos = create_output_stream(suppression_output);
            //writes ppm image
            nmsos << "P6\n" << width << " " << height << "\n255\n";
            for (int i = width * height - 1; i >= 0; i--) {
                nmsos << (unsigned char)nms[i];
                nmsos << (unsigned char)nms[i];
                nmsos << (unsigned char)nms[i];
            }

            if (is_debug) {
                cerr << "Non-maximum suppression writing time: " << float(clock () - nms_write_begin_time) /  CLOCKS_PER_SEC << endl;
            }

            delete_output_stream(nmsos);

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
            os << out[i];
            os << out[i];
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
