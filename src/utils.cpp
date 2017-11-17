// stl
#include <typeinfo>     // typeid
#include <fstream>      // ofstream
#include <string>       // string
#include <map>

#include "utils.hpp"

extern bool is_debug;

ostream& create_output_stream(string name) {
    if (!name.compare("-")) {
        return cout;
    } else if (!name.compare("log")) {
        if (is_debug) {
            return cerr;
        } else {
            ostream *nstream = new null_stream();
            return *nstream;
        }
    } else if (!name.compare("null")) {
        ostream *nstream = new null_stream();
        return *nstream;
    } else {
        ofstream *file = new ofstream();
        file->open(name, ios::trunc);
        return *file;
    }
}

void delete_output_stream(ostream& os) {
    if (typeid(os) == typeid(ofstream)) {
        ofstream& f1 = static_cast<ofstream&>(os);
        f1.close();
        delete &f1;
    } else if (typeid(os) == typeid(ostream)) {
        //TODO: check if null stream
        //delete &os;
        if (is_debug) {
            cerr << "The null output stream pointer has been deleted." << endl;
        }
    }
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

// boost
//#include <algorithm>                       //copy
//#include <ios>                             // ios_base::beg
//#include <iosfwd>                          // streamsize
//#include <boost/iostreams/stream.hpp>
//#include <boost/iostreams/categories.hpp>  // source_tag, sink_tag

//namespace io = boost::iostreams;

//template<typename Container>
//class DebugOutput {
//public:
//    // boost standart value_type, size_type and category in boost/iostreams/categories.hpp
//    typedef typename Container::value_type char_type;
//    typedef io::sink_tag category;
//
//    //constructor
//    DebugOutput(Container& container): container_(container) {}
//
//    streamsize write(const char_type* s, streamsize n);
//private:
//    Container container_;
//};

//template<typename Container>
//streamsize DebugOutput<Container>::write(const char_type* s, streamsize n) {
//    container_.insert(container_.end(), s, s + n);
//}
