#include <algorithm>                       //copy
// boost
//#include <ios>                             // ios_base::beg
#include <iosfwd>                          // streamsize
//#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/categories.hpp>  // source_tag, sink_tag

using namespace std;
namespace io = boost::iostreams;

template<typename Container>
class DebugOutput {
public:
    // boost standart value_type, size_type and category in boost/iostreams/categories.hpp
    typedef typename Container::value_type char_type;
    typedef io::sink_tag category;

    //constructor
    DebugOutput(Container& container): container_(container) {}

    streamsize write(const char_type* s, streamsize n);
private:
    Container container_;
};