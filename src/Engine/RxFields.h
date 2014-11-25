#ifndef RX_FIELDS_H
#define RX_FIELDS_H

#include <unordered_map>
#include <vector>
#include "Ray.h"
#include "Complex.h"

// Rx field
struct RxField
{
    ComplexVector field;
    double offset; // offset from the ray to the center of the rx sphere

    RxField(const ComplexVector &field, double offset) : field(field), offset(offset) {}
};

class RxFields
{
private:
    std::unordered_map<RayPath, std::vector<RxField> > mapping;

public:
    void AddField(const ComplexVector &field, const RayPath &path, double offset);
    ComplexVector Sum();
};

#endif
