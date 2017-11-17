//stl
#include <iostream>
#include <vector>
#include <array>
#include <map>

// omp
#include <omp.h>

#include <utils.hpp>
#include <canny.hpp>
#include <harris.hpp>

using namespace std;

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
const string EXAMPLE_CANNY("canny");
const string EXAMPLE_HARRIS("harris");
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
const int HYSTERESIS_MIN_DEF = 10;
const string HYSTERESIS_MAX = "-hma";
const int HYSTERESIS_MAX_DEF = 50;

void print_help() {
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
    cout << "-e: example's name, default: " << EXAMPLE_CANNY << endl;
    cout << "examples: " << EXAMPLE_CANNY << "(Canny edge detector), " << EXAMPLE_HARRIS << "(Harris corner detector), " << EXAMPLE_BOOST << endl;
    exit(0);
}

void example_boost() {
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

// global variables
map <string, string> argmap;

bool is_debug = false;
int width = 0;
int height = 0;
int frames = 0;

string derivatives_x_output;
string derivatives_y_output;
string derivatives_xy_output;
string suppression_output;
string gaussian_output;
string output;

int main(int argc, char *argv[]) {
    // creates map from input parameters
    for (int i = 0; i < argc; i++) {
        if (argv[i][0] == '-') {
            argmap[argv[i]] = i + 1 < argc? argv[i + 1]: "";
        }
    }

    width = read_value(argmap, WIDTH, WIDTH_DEF);
    height = read_value(argmap, HEIGHT, HEIGHT_DEF);
    frames = read_value(argmap, FRAMES, FRAMES_DEF);

    is_debug = read_value(argmap, DEBUG, false);
    output = (string)read_value(argmap, OUTPUT, OUTPUT_DEF);
    gaussian_output = (string)read_value(argmap, GAUSSIAN_OUTPUT, GAUSSIAN_OUTPUT_DEF);
    derivatives_x_output = (string)read_value(argmap, DERIVATIVES_X_OUTPUT, DERIVATIVES_X_OUTPUT_DEF);
    derivatives_y_output = (string)read_value(argmap, DERIVATIVES_Y_OUTPUT, DERIVATIVES_Y_OUTPUT_DEF);
    derivatives_xy_output = (string)read_value(argmap, DERIVATIVES_XY_OUTPUT, DERIVATIVES_XY_OUTPUT_DEF);
    suppression_output = (string)read_value(argmap, SUPPRESSION_OUTPUT, SUPPRESSION_OUTPUT_DEF);

    bool is_help = read_value(argmap, HELP, false);
    if (is_help) print_help();

    string example = read_value(argmap, EXAMPLE, EXAMPLE_CANNY);
    if (!example.compare(EXAMPLE_CANNY)) {
        int hmin = read_value(argmap, HYSTERESIS_MIN, HYSTERESIS_MIN_DEF);
        int hmax = read_value(argmap, HYSTERESIS_MAX, HYSTERESIS_MAX_DEF);
        example_canny(hmin, hmax);
    } else if (!example.compare(EXAMPLE_HARRIS)) { 
        example_harris();
    } else if (!example.compare(EXAMPLE_BOOST)) {
        example_boost();
    } else {
        cerr << "unknown argument: " << example << endl;
        print_help();
    }
}
