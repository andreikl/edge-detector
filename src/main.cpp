#include <typeinfo> // typeid
#include <iostream>
#include <fstream>
#include <cstring> // memcpy
#include <vector>
#include <array>
#include <map>

#include <math.h> // exp, pow

// omp
#include <omp.h>

using namespace std;

typedef unsigned char pixel_t;

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
const string OUTPUT("-o");
const string OUTPUT_DEF = "-";
const string DEBUG("-d");
const string EXAMPLE("-e");
const string EXAMPLE_MAIN("main");
const string EXAMPLE_BOOST("boost");
const string HELP("--help");

map <string, string> argmap;
bool isDebug = false;

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
    cout << "-w: width, default: " << WIDTH_DEF << endl;
    cout << "-h: height, default: " << HEIGHT_DEF << endl;
    cout << "-f: frame, meaning of this argument depends on example, default: " << FRAMES_DEF << endl;
    cout << "-o: output, default: " << OUTPUT_DEF << ", -: output stream" << endl;
    cout << "-e: testName, default: " << EXAMPLE_MAIN << endl;
    cout << "tests: " << EXAMPLE_MAIN << ", " << EXAMPLE_BOOST << endl;
    exit(0);
}

void convolution(const pixel_t *in, pixel_t *out, int width, int height, const float *kernel, int ksize) {
    int khalf = ksize / 2;
    int kres = ksize % 2;

    auto calcConvolution = [=](int row, int column) -> pixel_t {
        float pixel = 0.0;
        size_t c = 0;
        for (int i = -khalf; i < khalf + kres; i++) { //rows
            for (int j = -khalf; j < khalf + kres; j++) { //columns
                pixel += in[(row - i) * width + column - j] * kernel[c];
                /*if (isDebug && row == 5 && column == 5) {
                    int index = (row - i) * width + column - j;
                    cout << "kernel i: " << kernel[c] << endl;
                    cout << "image [row, column]: [" << index / width << ", " << index % width << "]";
                }*/
                c++;
            }
        }
        return static_cast<pixel_t>(pixel);
    };

    for (int i = khalf + kres; i < height - khalf; i++) { //rows
        for (int j = khalf + kres; j < width - khalf; j++) { //columns
            out[i * width + j] = calcConvolution(i, j);
        }
    }
}

void gaussian_filter(const pixel_t *in, pixel_t *out, int width, int height, float sigma) {
    const int ksize = 2 * (int)(2 * sigma) + 3;
    const float mean = (float)floor(ksize / 2.0);
    float kernel[ksize * ksize];

    size_t c = 0;
    for (int i = 0; i < ksize; i++) {
        for (int j = 0; j < ksize; j++) {
            kernel[c] = exp(-0.5 * (pow((i - mean) / sigma, 2.0) +
                                    pow((j - mean) / sigma, 2.0)))
                        / (2 * M_PI * sigma * sigma);
            c++;
        }
    }

    convolution(in, out, width, height, kernel, ksize);
}

void contrast_filter(pixel_t *inout, int width, int height) {
    //todo: min and max convert to lambda
    pixel_t min = inout[0];
    pixel_t max = inout[0];

    #pragma omp parallel for
    for (int i = 0; i < width * height; i++) {
        if (inout[i] < min) {
            min = inout[i];
        } 
        if (inout[i] > max) {
            max = inout[i];
        } 
    }
    //cout << "min: " << static_cast<int>(min.gr) << endl;
    //cout << "max: " << static_cast<int>(max.gr) << endl;

    auto adjastColour = [=](pixel_t colour) -> pixel_t {
        // color / min - max = x / 255
        return static_cast<pixel_t>((colour - min) * 255 / (float)(max - min));
    };

    #pragma omp parallel for
    for (int i = 0; i < width * height; i++) {
        inout[i] = adjastColour(inout[i]);
    }
}

void exampleMain() {
    // converts input parameters to variables
    const int width = readValue(WIDTH, WIDTH_DEF);
    const int height = readValue(HEIGHT, HEIGHT_DEF);
    const int frames = readValue(FRAMES, FRAMES_DEF);
    const string output = readValue(OUTPUT, OUTPUT_DEF);

    pixel_t image[width * height];

    // initialise output stream
    ostream& os = ([=]() -> ostream& {
        if (!output.compare("-")) {
            return cout;
        } else {
            ofstream* file = new ofstream();
            file->open(output, ios::trunc);
            return *file;
        }
    })();
    
    for (int k = 0; k < frames; k++) {
        //reads one frame to array
        cin.read((char*)(image), width * height * sizeof(pixel_t));

        if (k == frames - 1) {
            //exampleContrast(image, width, height);

            pixel_t image_gaussian[width * height];
            memcpy(image_gaussian, image, width * height * sizeof(pixel_t));
            gaussian_filter(image, image_gaussian, width, height, 2.0);
            contrast_filter(image_gaussian, width, height);

            const int d_size = 3;
            const float dx[] = {
                -1, 0, 1,
                -2, 0, 2,
                -1, 0, 1,
            };
            const float dy[] = {
                1, 2, 1,
                0, 0, 0,
                -1, -2, -1,
            };

            pixel_t image_dx[width * height];
            memcpy(image_dx, image_gaussian, width * height * sizeof(pixel_t));
            convolution(image_gaussian, image_dx, width, height, dx, d_size);

            pixel_t image_dy[width * height];
            memcpy(image_dy, image_gaussian, width * height * sizeof(pixel_t));
            convolution(image_gaussian, image_dy, width, height, dy, d_size);

            //#pragma omp parallel for
            for (int i = 0; i < width * height; i++) {
                image[i] = (pixel_t)hypot(image_dx[i], image_dy[i]);
            }

            //writes ppm image
            os << "P6\n" << width << " " << height << "\n255\n";
            //std::cout.write((char*)image.data(), width * height * sizeof(RGB));
            for (int i = width * height - 1; i >= 0; i--) {
                os << image[i];
                os << image[i];
                os << image[i];
            }
        }
    }

    if (typeid(os) == typeid(ofstream)) {
        ofstream& f1 = static_cast<ofstream&>(os);
        f1.close();
        delete &f1;
        //cout << "file has been closed";
    }
}

void exampleBoost() {
    int i;
    int threadID = 0;
    #pragma omp parallel for private(i, threadID)
    for (i = 0; i < 16; i++) {
        threadID = omp_get_thread_num();
        # pragma omp critical
        {
            cout << "Thread: " << threadID << endl;
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

    const bool isHelp = readValue(HELP, false);
    if (isHelp) printHelp();
    

    if (!example.compare(EXAMPLE_MAIN)) {
        exampleMain();
    } else if (!example.compare(EXAMPLE_BOOST)) {
        exampleBoost();
    } else {
        cout << "unknown argument: " << example << endl;
        printHelp();
    }
}
