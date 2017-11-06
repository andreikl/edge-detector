#include <string>       // string
#include <fstream>      // ofstream
#include <iostream>     // ostream
#include <typeinfo>     // typeid

#include "utils.hpp"

extern bool isDebug;

ostream& create_output_stream(string name) {
    if (!name.compare("-")) {
        return cout;
    } else if (!name.compare("null")) {
        ostream *null_stream = new NullStream();
        return *null_stream;
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
        //if (isDebug) {
        //    cerr << "The output file has been closed and pointer has been deleted." << endl;
        //}
    } else if (typeid(os) == typeid(ostream)) {
        delete &os;
        if (isDebug) {
            cerr << "The null output stream pointer has been deleted." << endl;
        }
    }
}

//template<typename Container>
//streamsize DebugOutput<Container>::write(const char_type* s, streamsize n) {
//    container_.insert(container_.end(), s, s + n);
//}
