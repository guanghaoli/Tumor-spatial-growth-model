#ifndef sort_function_hpp
#define sort_function_hpp

#include <stdio.h>
#include <algorithm>

using namespace std;

struct element
{
    double data;
    int index;
};
int compare(const void *a, const void *b)
{
    return (*(const element*)a).data > (*(const element*)b).data?1:-1;
}

#endif
