#include "RxFields.h"

void RxFields::AddField(const ComplexVector &field, const RayPath &path, double offset)
{
    if (mapping.count(path) == 0)
    {
        mapping[path] = std::vector<RxField>();
    }
    mapping[path].push_back(RxField(field, offset));
}

ComplexVector RxFields::Sum()
{
    ComplexVector sum(
        ComplexNumber(0, 0),
        ComplexNumber(0, 0),
        ComplexNumber(0, 0));

    if (mapping.size() == 0)
    {
        return sum;
    }
    else
    {
        std::unordered_map<RayPath, std::vector<RxField> >::iterator it;
        for (it = mapping.begin(); it != mapping.end(); ++it) // for each path
        {
            // use the field with the min distance
            double minOffset = DBL_MAX;
            ComplexVector minField(
                ComplexNumber(0, 0),
                ComplexNumber(0, 0),
                ComplexNumber(0, 0));

            for (unsigned int i = 0; i < it->second.size(); i++)
            {
                if (it->second[i].offset < minOffset)
                {
                    minOffset = it->second[i].offset;
                    minField = it->second[i].field;
                }
            }

            sum = sum + minField;
        }
        return sum;
    }
}
