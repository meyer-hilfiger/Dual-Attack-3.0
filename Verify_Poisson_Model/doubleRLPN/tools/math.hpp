#ifndef HEADER_MATH
#define HEADER_MATH
constexpr size_t power(size_t x, size_t p) {
    size_t result = 1;
    while (p) {
        if (p & 0x1) {
            result *= x;
        }
        x *= x;
        p >>= 1;
    }
    return result;
}
#endif