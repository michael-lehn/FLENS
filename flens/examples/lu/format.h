#ifndef LU_FORMAT_H
#define LU_FORMAT_H 1

#include <memory>
#include <iostream>
#include <string>
#include <cstdio>

namespace flens {

template<typename ... Args>
std::string
format(std::string fmt, Args ... args)
{
    size_t size = 1 + snprintf(nullptr, 0, fmt.c_str(), args ...);
    std::unique_ptr<char[]> buf(new char[size]);
    snprintf(buf.get(), size, fmt.c_str(), args ...);

    return std::string(buf.get(), buf.get() + size-1);
}

} // namespace flens

#endif
