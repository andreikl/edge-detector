#include <math.h>       // exp, pow

// stl
#include <typeinfo>                          // typeid
#include <iostream>                          // ostream
#include <streambuf>                         // streambuf
#include <map>

typedef unsigned char pixel_t;
typedef short derivatives_t;
#define MAX_BRIGHTNESS 255

using namespace std;

ostream& create_output_stream(string name);
void delete_output_stream(ostream& os);
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
};

void contrast_filter(pixel_t *inout, int width, int height);

struct return_proxy {
    return_proxy (int val) : val_int(val) {}
    return_proxy (string val) : val_str(val) {}
    return_proxy (bool val) : val_bool(val) {}

    operator int() const {
        return val_int;
    }

    operator bool() const {
        return val_bool;
    }

    operator string() const {
        return val_str;
    }

private:
    int val_int;
    bool val_bool;
    string val_str;
};

class null_stream: public streambuf, public ostream {
public:
    null_stream(): ostream( this ) {}
    int overflow(int c) { return c; }
};

template<typename T_IN, typename T_OUT> class image_filter {
protected:
    image_filter(string output_name): output_name(output_name), log(create_output_stream("log")) {
    }

public:
    ~image_filter() {
        delete_output_stream(log);
    }

protected:
    void process_before() {
        start_time = clock();
    }

    virtual void process(const T_IN *in, T_OUT *out, int width, int height) = 0;

    void process_after(const T_IN *out, int width, int height) {
        log << typeid(this).name() <<  " filter process time: " << float(clock () - start_time) /  CLOCKS_PER_SEC << endl;
        
        start_time = clock();
        ostream& output = create_output_stream(output_name);

        //writes ppm image
        output << "P6\n" << width << " " << height << "\n255\n";
        for (int i = width * height - 1; i >= 0; i--) {
            output << (unsigned char)out[i];
            output << (unsigned char)out[i];
            output << (unsigned char)out[i];
        }

        log << typeid(this).name() <<  " filter writing time: " << float(clock () - start_time) /  CLOCKS_PER_SEC << endl;
        delete_output_stream(output);
    }

public:
    void do_process(const T_IN *in, T_OUT *out, int width, int height) {
        process_before();
        process(in, out, width, height);
        process_after(out, width, height);
    }

protected:
    ostream& log;
    string output_name;
    clock_t start_time;
};

//void gaussian_filter(const pixel_t *in, pixel_t *out, int width, int height, double sigma, int size);
class gaussian_filter: public image_filter<pixel_t, pixel_t> {
public:
    gaussian_filter(string output_name, double sigma, int size): image_filter(output_name), sigma(sigma), size(size) {
        log << "Gaussian filter: sigma, size: " << sigma << ", " << size << endl;
    }

protected:
    void process(const pixel_t *in, pixel_t *out, int width, int height) {
        double kernel[size * size];
        double mean = size / 2.0;
        double sum = 0.0;

        size_t c = 0;
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                kernel[c] = exp(-0.5 * (pow((i - mean) / sigma, 2.0) +
                                        pow((j - mean) / sigma, 2.0)))
                            / (2 * M_PI * sigma * sigma);
                sum += kernel[c];

                log << i << ", " << j << ": " << kernel[c] << endl;

                c++;
            }
        }

        for (int i = 0; i < c; i++) {
            kernel[i] /= sum;
        }

        log << "Gaussian filter preparation: " << float(clock() - start_time) /  CLOCKS_PER_SEC << endl;

        // reset start time
        start_time = clock();

        convolution<pixel_t, pixel_t>(in, out, width, height, kernel, size, size);
    }

private:
    double sigma;
    int size;
};

return_proxy read_value(map<string, string> argmap, string name, auto def_val) {
    if (argmap.count(name)) {
        if (typeid(def_val) == typeid(int)) {
            return return_proxy(stoi(argmap.at(name)));
        } else if (typeid(def_val) == typeid(string)) {
            return return_proxy(argmap.at(name));
        } else if (typeid(def_val) == typeid(bool)) {
            return return_proxy(true);
        }
    } else {
        return return_proxy(def_val);
    }
};
