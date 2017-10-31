#include "utils.hpp"

extern bool isDebug;

template<typename Container>
streamsize DebugOutput<Container>::write(const char_type* s, streamsize n) {
    container_.insert(container_.end(), s, s + n);
}
