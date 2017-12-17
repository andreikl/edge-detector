#include <math.h>       // exp, pow

// stl
#include <typeinfo>                          // typeid
#include <iostream>                          // ostream
#include <streambuf>                         // streambuf
#include <memory>                            // auto_ptr
#include <map>

typedef unsigned char pixel_t;
typedef short derivatives_t;
#define MAX_BRIGHTNESS 255

using namespace std;

ostream& create_output_stream(string name);
void delete_output_stream(ostream& os);

template<typename T_IN, typename T_OUT>
void convolution(const T_IN *in, T_OUT *out, int width, int height, const double *kernel, int kwidth, int kheight) {
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
    return_proxy (double val) : val_double(val) {}
    return_proxy (string val) : val_str(val) {}
    return_proxy (bool val) : val_bool(val) {}
    return_proxy (int val) : val_int(val) {}

    operator int() const {
        return val_int;
    }

    operator bool() const {
        return val_bool;
    }

    operator string() const {
        return val_str;
    }

    operator double() const {
        return val_double;
    }

private:
    int val_int;
    bool val_bool;
    string val_str;
    double val_double;
};

class null_stream: public streambuf, public ostream {
public:
    null_stream(): ostream( this ) {}
    int overflow(int c) { return c; }
};

template<typename T_IN, typename T_OUT>
class image_filter {
protected:
    image_filter(string output_name, bool is_debug, int debug_x, int debug_y):
        output_name(output_name),
        is_debug(is_debug),
        debug_x(debug_x),
        debug_y(debug_y),
        log(create_output_stream("log")) {
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

    void process_after(const T_OUT *out, int width, int height) {
        log << typeid(this).name() <<  " filter process time: " << float(clock () - start_time) /  CLOCKS_PER_SEC << endl;

        T_OUT min = out[0];
        T_OUT max = out[0];

        //#pragma omp parallel for
        for (int i = 0; i < width * height; i++) {
            if (out[i] < min) {
                min = out[i];
            } 
            if (out[i] > max) {
                max = out[i];
            } 
        }

        start_time = clock();
        ostream& output = create_output_stream(output_name);

        //writes ppm image
        output << "P6\n" << width << " " << height << "\n255\n";
        for (int i = width * height - 1; i >= 0; i--) {
            unsigned char colour = static_cast<unsigned char>((out[i] - min) * MAX_BRIGHTNESS / (double)(max - min));

            if (this->is_debug && i % width == this->debug_x && i / width == this->debug_y) {
                log << typeid(this).name() <<  " point output[" << this->debug_x << ", " << this->debug_y << "]: "<< out[i] << ", c: " <<  (int)colour << endl;
                output << (unsigned char)MAX_BRIGHTNESS;
                output << (unsigned char)0;
                output << (unsigned char)0;
            } else {
                output << colour;
                output << colour;
                output << colour;
            }
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

protected:
    bool is_debug;
    int debug_x;
    int debug_y;
};

template<typename T_IN, typename T_OUT>
class gaussian_filter: public image_filter<T_IN, T_OUT> {
public:
    gaussian_filter(string output_name, double sigma, int size): image_filter<T_IN, T_OUT>::image_filter(output_name, false, -1, -1), sigma(sigma), size(size) {
        image_filter<T_IN, T_OUT>::log << "Gaussian filter: sigma, size: " << sigma << ", " << size << endl;
    }

    gaussian_filter(string output_name, double sigma, int size, bool is_debug, int debug_x, int debug_y): image_filter<T_IN, T_OUT>::image_filter(output_name, is_debug, debug_x, debug_y), sigma(sigma), size(size) {
        image_filter<T_IN, T_OUT>::log << "Gaussian filter: sigma, size: " << sigma << ", " << size << endl;
    }

protected:
    void process(const T_IN *in, T_OUT *out, int width, int height) {
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

                image_filter<T_IN, T_OUT>::log << i << ", " << j << ": " << kernel[c] << endl;

                c++;
            }
        }

        for (int i = 0; i < c; i++) {
            kernel[i] /= sum;
        }

        image_filter<T_IN, T_OUT>::log << "Gaussian filter preparation: " << float(clock() - image_filter<T_IN, T_OUT>::start_time) /  CLOCKS_PER_SEC << endl;

        // reset start time
        image_filter<T_IN, T_OUT>::start_time = clock();

        convolution<T_IN, T_OUT>(in, out, width, height, kernel, size, size);
    }

private:
    double sigma;
    int size;
};

class convolution_filter: public image_filter<pixel_t, derivatives_t> {
public:
    convolution_filter(string output_name, const double conv[], int size_x, int size_y): image_filter(output_name, false, -1, -1), conv(conv), size_x(size_x), size_y(size_y) {
        log << "Convolution filter: convolution, size: " << conv << ", " << size_x << ", " << size_y << endl;
    }

    convolution_filter(string output_name, const double conv[], int size_x, int size_y, bool is_debug, int debug_x, int debug_y): image_filter(output_name, is_debug, debug_x, debug_y), conv(conv), size_x(size_x), size_y(size_y) {
        log << "Convolution filter: convolution, size: " << conv << ", " << size_x << ", " << size_y << endl;
    }

protected:
    void process(const pixel_t *in, derivatives_t *out, int width, int height) {
        convolution<pixel_t, derivatives_t>(in, out, width, height, conv.get(), size_x, size_y);
    }

private:
    unique_ptr<const double> conv;
    int size_x;
    int size_y;
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
