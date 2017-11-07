#include <typeinfo>     // typeid
#include <iostream>
#include <cstring>      // memcpy
#include <vector>

#include <string>       // string
#include <array>
#include <ctime>        // clock
#include <map>

#include <math.h>       // exp, pow

// omp
#include <omp.h>

#include <utils.hpp>

using namespace std;

typedef unsigned char pixel_t;
typedef short derivatives_t;

#define MAX_BRIGHTNESS 255

typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
    unsigned char nothing;
} rgb_t;

const string WIDTH("-w");
const int WIDTH_DEF = 640;
const string HEIGHT("-h");
const int HEIGHT_DEF = 480;
const string FRAMES("-f");
const int FRAMES_DEF = 1;
const string DEBUG("-d");
const string EXAMPLE("-e");
const string EXAMPLE_MAIN("main");
const string EXAMPLE_BOOST("boost");
const string HELP("--help");

const string OUTPUT("-o");
const string OUTPUT_DEF = "-";
const string GAUSSIAN_OUTPUT("-g");
const string GAUSSIAN_OUTPUT_DEF = "null";
const string DERIVATIVES_X_OUTPUT("-dx");
const string DERIVATIVES_X_OUTPUT_DEF = "null";
const string DERIVATIVES_Y_OUTPUT("-dy");
const string DERIVATIVES_Y_OUTPUT_DEF = "null";
const string DERIVATIVES_XY_OUTPUT("-dxy");
const string DERIVATIVES_XY_OUTPUT_DEF = "null";
const string SUPPRESSION_OUTPUT("-nms");
const string SUPPRESSION_OUTPUT_DEF = "null";
const string HYSTERESIS_MIN = "-hmi";
const int HYSTERESIS_MIN_DEF = 45;
const string HYSTERESIS_MAX = "-hma";
const int HYSTERESIS_MAX_DEF = 50;

map <string, string> argmap;

bool isDebug = false;

string derivatives_x_output;
string derivatives_y_output;
string derivatives_xy_output;
string suppression_output;
string gaussian_output;
string output;

struct ReturnProxy {
    ReturnProxy (int val) : valInt(val) {}
    ReturnProxy (string val) : valString(val) {}
    ReturnProxy (bool val) : valBool(val) {}

    operator int() const {
        return valInt;
    }

    operator bool() const {
        return valBool;
    }

    operator string() const {
        return valString;
    }

private:
    int valInt;
    bool valBool;
    string valString;
};

ReturnProxy readValue(string name, auto defValue) {
    if (argmap.count(name)) {
        if (typeid(defValue) == typeid(int)) {
            return ReturnProxy(stoi(argmap.at(name)));
        } else if (typeid(defValue) == typeid(string)) {
            return ReturnProxy(argmap.at(name));
        } else if (typeid(defValue) == typeid(bool)) {
            return ReturnProxy(true);
        }
    } else {
        return ReturnProxy(defValue);
    }
};

void printHelp() {
    cout << "edge-detector [options]" << endl;
    cout << "options:" << endl;
    cout << "--help: help" << endl;
    cout << "-d: show debug info" << endl;
    cout << "-w: width, default: " << WIDTH_DEF << endl;
    cout << "-h: height, default: " << HEIGHT_DEF << endl;
    cout << "-f: frame, meaning of this argument depends on example, default: " << FRAMES_DEF << endl;
    cout << GAUSSIAN_OUTPUT << ": gaussian output, default: " << GAUSSIAN_OUTPUT_DEF << ", null: void output" << endl;
    cout << DERIVATIVES_X_OUTPUT << ": derivatives x output, default: " << DERIVATIVES_X_OUTPUT_DEF << ", null: void output" << endl;
    cout << DERIVATIVES_Y_OUTPUT << ": derivatives y output, default: " << DERIVATIVES_Y_OUTPUT_DEF << ", null: void output" << endl;
    cout << DERIVATIVES_XY_OUTPUT << ": derivatives xy output, default: " << DERIVATIVES_XY_OUTPUT_DEF << ", null: void output" << endl;
    cout << SUPPRESSION_OUTPUT << ": Non-maximum supression output, default: " << SUPPRESSION_OUTPUT_DEF << ", null: void output" << endl;
    cout << HYSTERESIS_MIN << ": Hysteresis min, default: " << HYSTERESIS_MIN_DEF << endl;
    cout << HYSTERESIS_MAX << ": Hysteresis max, default: " << HYSTERESIS_MAX_DEF << endl;
    cout << OUTPUT << ": output, default: " << OUTPUT_DEF << ", -: output stream" << endl;
    cout << "-e: testName, default: " << EXAMPLE_MAIN << endl;
    cout << "tests: " << EXAMPLE_MAIN << ", " << EXAMPLE_BOOST << endl;
    exit(0);
}

template<typename T_IN, typename T_OUT> void convolution(const T_IN *in, T_OUT *out, int width, int height, const double *kernel, int kwidth, int kheight) {
    int kwidthhalf = kwidth / 2;
    int kwidthmodulo = kwidth % 2;
    int kheighthalf = kheight / 2;
    int kheightmodulo = kheight % 2;

    /*auto calcConvolution = [=](int row, int column) -> pixel_t {
        double pixel = 0.0;
        size_t c = 0;
        for (int ii = -kheighthalf; ii < kheighthalf + kheightmodulo; ii++) { //rows
            for (int jj = -kwidthhalf; jj < kwidthhalf + kwidthmodulo; jj++) { //columns
                pixel += in[(i - ii) * width + j - jj] * kernel[c];
                c++;
            }
        }
        return static_cast<pixel_t>(pixel);
    };*/

    for (int i = kheighthalf + kheightmodulo; i < height - kheighthalf; i++) { //rows
        for (int j = kwidthhalf + kwidthmodulo; j < width - kwidthhalf; j++) { //columns
            // the perfomance without lambda is 20% faster
            double pixel = 0.0;
            size_t c = 0;
            for (int ii = -kheighthalf; ii < kheighthalf + kheightmodulo; ii++) { //rows
                for (int jj = -kwidthhalf; jj < kwidthhalf + kwidthmodulo; jj++) { //columns
                    pixel += in[(i - ii) * width + j - jj] * kernel[c];
                    /*if (isDebug && i == 20 && j == 91) {
                        int index = (i - ii) * width + j - jj;
                        cerr << "kernel i: " << kernel[c] << endl;
                        cerr << "image [row, column]: [" << index / width << ", " << index % width << "] " << (int)in[(i - ii) * width + j - jj] << endl;
                        cerr << "sum: " << pixel << endl;
                    }

                    if (isDebug && i == 20 && j == 25) {
                        int index = (i - ii) * width + j - jj;
                        cerr << "kernel i: " << kernel[c] << endl;
                        cerr << "image [row, column]: [" << index / width << ", " << index % width << "] " << (int)in[(i - ii) * width + j - jj] << endl;
                        cerr << "sum: " << pixel << endl;
                    }*/
                    c++;
                }
            }
            out[i * width + j] = static_cast<T_OUT>(pixel);
            //out[i * width + j] = calcConvolution(i, j);
        }
    }
}

void gaussian_filter(const pixel_t *in, pixel_t *out, int width, int height, double sigma, int size) {
    const clock_t preparation_begin_time = clock();

    double kernel[size * size];
    double mean = size / 2.0;
    double sum = 0.0;

    if (isDebug) {
        cerr << "Gaussian filter size: " << size << endl;
    }
    size_t c = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            kernel[c] = exp(-0.5 * (pow((i - mean) / sigma, 2.0) +
                                    pow((j - mean) / sigma, 2.0)))
                        / (2 * M_PI * sigma * sigma);
            sum += kernel[c];

            if (isDebug) {
                cerr << i << ", " << j << ": " << kernel[c] << endl;
            }

            c++;
        }
    }

    for (int i = 0; i < c; i++) {
        kernel[i] /= sum;
    }

    if (isDebug) {
        cerr << "Gaussian filter preparation: " << float(clock () - preparation_begin_time) /  CLOCKS_PER_SEC << endl;
    }

    const clock_t gaussian_begin_time = clock();
    convolution<pixel_t, pixel_t>(in, out, width, height, kernel, size, size);

    if (isDebug) {
        cerr << "Gaussian filter calculation: " << float(clock () - gaussian_begin_time) /  CLOCKS_PER_SEC << endl;
    }

    const clock_t write_begin_time = clock();
    ostream& gos = create_output_stream(gaussian_output);

    //writes ppm image
    gos << "P6\n" << width << " " << height << "\n255\n";
    for (int i = width * height - 1; i >= 0; i--) {
        gos << out[i];
        gos << out[i];
        gos << out[i];
    }

    if (isDebug) {
        cerr << "Gaussian filter writing time: " << float(clock () - write_begin_time) /  CLOCKS_PER_SEC << endl;
    }

    delete_output_stream(gos);
}

void contrast_filter(pixel_t *inout, int width, int height) {
    //todo: min and max convert to lambda
    pixel_t min = inout[0];
    pixel_t max = inout[0];

    //#pragma omp parallel for
    for (int i = 0; i < width * height; i++) {
        if (inout[i] < min) {
            min = inout[i];
        } 
        if (inout[i] > max) {
            max = inout[i];
        } 
    }
    //cerr << "min: " << static_cast<int>(min.gr) << endl;
    //cerr << "max: " << static_cast<int>(max.gr) << endl;

    auto adjastColour = [=](pixel_t colour) -> pixel_t {
        // color / min - max = x / MAX_BRIGHTNESS
        return static_cast<pixel_t>((colour - min) * MAX_BRIGHTNESS / (float)(max - min));
    };

    //#pragma omp parallel for
    for (int i = 0; i < width * height; i++) {
        inout[i] = adjastColour(inout[i]);
    }
}

void exampleMain(unsigned char tmin, unsigned char tmax) {
    // converts input parameters to variables
    const int width = readValue(WIDTH, WIDTH_DEF);
    const int height = readValue(HEIGHT, HEIGHT_DEF);
    const int frames = readValue(FRAMES, FRAMES_DEF);

    pixel_t image[width * height];
    pixel_t image_gaussian[width * height];
    pixel_t out[width * height];            

    // initialise output stream
    ostream& os = create_output_stream(output);
    /*ostream& os = ([=]() -> ostream& {
        if (!output.compare("-")) {
            return cerr;
        } else {
            ofstream* file = new ofstream();
            file->open(output, ios::trunc);
            return *file;
        }
    })();*/
    
    for (int k = 0; k < frames; k++) {
        //reads one frame to array
        cin.read((char*)(image), width * height * sizeof(pixel_t));

        if (k == frames - 1) {
            //exampleContrast(image, width, height);

            const clock_t calc_begin_time = clock();

            memcpy(image_gaussian, image, width * height * sizeof(pixel_t));
            gaussian_filter(image, image_gaussian, width, height, 0.8, 5);
            //contrast_filter(image_gaussian, width, height);

            // borders are -0.5 0.5: -128 to 128, -2, 2: -512 to 512
            //const double dx[] = {-2, 2};
            //const double dy[] = {2, -2};
            double dx[] = {-1, 0, 1,
                           -2, 0, 2,
                           -1, 0, 1};
            double dy[] = { 1, 2, 1,
                            0, 0, 0,
                           -1,-2,-1};

            derivatives_t image_dx[width * height];
            memset(image_dx, 0, width * height * sizeof(derivatives_t));
            convolution<pixel_t, derivatives_t>(image_gaussian, image_dx, width, height, dx, 3, 3);
            const clock_t dx_write_begin_time = clock();
            ostream& dxos = create_output_stream(derivatives_x_output);

            //writes ppm image
            dxos << "P6\n" << width << " " << height << "\n255\n";
            for (int i = width * height - 1; i >= 0; i--) {
                dxos << (unsigned char)(image_dx[i] + 127);
                dxos << (unsigned char)(image_dx[i] + 127);
                dxos << (unsigned char)(image_dx[i] + 127);
            }

            if (isDebug) {
                cerr << "Derivative X filter writing time: " << float(clock () - dx_write_begin_time) /  CLOCKS_PER_SEC << endl;
            }

            delete_output_stream(dxos);

            derivatives_t image_dy[width * height];
            memset(image_dy, 0, width * height * sizeof(derivatives_t));
            convolution<pixel_t, derivatives_t>(image_gaussian, image_dy, width, height, dy, 3, 3);
            const clock_t dy_write_begin_time = clock();
            ostream& dyos = create_output_stream(derivatives_y_output);

            //writes ppm image
            dyos << "P6\n" << width << " " << height << "\n255\n";
            for (int i = width * height - 1; i >= 0; i--) {
                dyos << (unsigned char)(image_dy[i] + 127);
                dyos << (unsigned char)(image_dy[i] + 127);
                dyos << (unsigned char)(image_dy[i] + 127);
            }

            if (isDebug) {
                cerr << "Derivative Y filter writing time: " << float(clock () - dy_write_begin_time) /  CLOCKS_PER_SEC << endl;
            }

            delete_output_stream(dyos);
            
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

            if (isDebug) {
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

            if (isDebug) {
                cerr << "Non-maximum suppression writing time: " << float(clock () - nms_write_begin_time) /  CLOCKS_PER_SEC << endl;
            }

            delete_output_stream(nmsos);


            //hysteresis
            //Reuse array, used as a stack of coordinates. width*height/4 elements should be enough
            int *edges = (int*)image_dx;
            memset(out, 0, width * height * sizeof(pixel_t));

            for (int i = 0; i < width * height; i++) {
                if (nms[i] >= tmax && out[i] == 0) { // trace edges
                    out[i] = MAX_BRIGHTNESS;

                    // initialising values
                    int nedges = 1;
                    edges[0] = i;
                    do {
                        // get previous position of line
                        nedges--;
                        const int t = edges[nedges];

                        int nbs[8]; // positions of neighbours
                        nbs[0] = t - width;  // nn
                        nbs[1] = t + width;  // ss
                        nbs[2] = t + 1;      // ww
                        nbs[3] = t - 1;      // ee
                        nbs[4] = nbs[0] + 1; // nw
                        nbs[5] = nbs[0] - 1; // ne
                        nbs[6] = nbs[1] + 1; // sw
                        nbs[7] = nbs[1] - 1; // se

                        for (int k = 0; k < 8; k++) {
                            // add position to handle if it isn't handled yet
                            if (nms[nbs[k]] >= tmin && out[nbs[k]] == 0) {
                                out[nbs[k]] = MAX_BRIGHTNESS;
                                edges[nedges] = nbs[k];
                                nedges++;
                            }
                        }
                    } while (nedges > 0);
                }
            }
            
            if (isDebug) {
                cerr << "Calculation time: " << float(clock () - calc_begin_time) /  CLOCKS_PER_SEC << endl;
            }
        }
    }

    const clock_t write_begin_time = clock();

    //writes ppm image
    os << "P6\n" << width << " " << height << "\n255\n";
    //cerr.write((char*)image.data(), width * height * sizeof(RGB));
    for (int i = width * height - 1; i >= 0; i--) {
        if (out[i] == MAX_BRIGHTNESS) {
        //if (false) {
            os << out[i];
            os << out[i];
            os << (unsigned char)0;
        } else {
            os << image[i];
            os << image[i];
            os << image[i];
        }
    }

    if (isDebug) {
        cerr << "Writing time: " << float(clock () - write_begin_time) /  CLOCKS_PER_SEC << endl;
    }

    delete_output_stream(os);
    /*if (typeid(os) == typeid(ofstream)) {
        ofstream& f1 = static_cast<ofstream&>(os);
        f1.close();
        delete &f1;
    }*/
}

void exampleBoost() {
    int i;
    int threadID = 0;
    #pragma omp parallel for private(i, threadID)
    for (i = 0; i < 16; i++) {
        threadID = omp_get_thread_num();
        # pragma omp critical
        {
            cerr << "Thread: " << threadID << endl;
        }
    }
}

// split to separate files
int main(int argc, char *argv[]) {
    // creates map from input parameters
    for (int i = 0; i < argc; i++) {
        if (argv[i][0] == '-') {
            argmap[argv[i]] = i + 1 < argc? argv[i + 1]: "";
        }
    }

    const string example = readValue(EXAMPLE, EXAMPLE_MAIN);
    isDebug = readValue(DEBUG, false);
    output = (string)readValue(OUTPUT, OUTPUT_DEF);
    gaussian_output = (string)readValue(GAUSSIAN_OUTPUT, GAUSSIAN_OUTPUT_DEF);
    derivatives_x_output = (string)readValue(DERIVATIVES_X_OUTPUT, DERIVATIVES_X_OUTPUT_DEF);
    derivatives_y_output = (string)readValue(DERIVATIVES_Y_OUTPUT, DERIVATIVES_Y_OUTPUT_DEF);
    derivatives_xy_output = (string)readValue(DERIVATIVES_XY_OUTPUT, DERIVATIVES_XY_OUTPUT_DEF);
    suppression_output = (string)readValue(SUPPRESSION_OUTPUT, SUPPRESSION_OUTPUT_DEF);
    int hmin = readValue(HYSTERESIS_MIN, HYSTERESIS_MIN_DEF);
    int hmax = readValue(HYSTERESIS_MAX, HYSTERESIS_MAX_DEF);

    const bool isHelp = readValue(HELP, false);
    if (isHelp) printHelp();
    

    if (!example.compare(EXAMPLE_MAIN)) {
        exampleMain(hmin, hmax);
    } else if (!example.compare(EXAMPLE_BOOST)) {
        exampleBoost();
    } else {
        cerr << "unknown argument: " << example << endl;
        printHelp();
    }
}
