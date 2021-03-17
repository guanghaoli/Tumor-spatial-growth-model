#ifndef initiation_visual_range_hpp
#define initiation_visual_range_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>

using namespace blitz;
void initiation_visual_range(Array<double,3> &Visual_range,int Visual_range_x,int Visual_range_y,int &N0,int R0){
    int r0 = R0/2;
    int center_x = Visual_range_x/2;
    int center_y = Visual_range_y/2;
    for(int x=1;x<=Visual_range_x;x++){
        for(int y =1;y<=Visual_range_y;y++){
            if(pow((x-center_x),2)+pow((y-center_y),2) <= pow(r0,2)){
                N0 += 1;
                Visual_range(x,y,1) = 1;
            }
        }
    }
}

#endif
