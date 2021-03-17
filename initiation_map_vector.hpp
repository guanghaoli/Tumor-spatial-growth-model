#ifndef initiation_map_vector_hpp
#define initiation_map_vector_hpp

#include <stdio.h>
#include <map>
#include <iostream>
#include <blitz/blitz.h>
#include <blitz/array.h>

using namespace blitz;
using namespace std;
void initiation_map_vector(Array<double,2> &cell_array,map<int,vector<int> > &clone_cell_number,map<int,vector<int> > &clone_max_driver,map<int,vector<double> > &clone_ave_driver){
    int ROW = cell_array.rows();
    cout << "initiate the map vector"<<endl;
    for(int ii=1; ii<=ROW; ii++){
        int cell_index = (int)cell_array(ii,5);
        clone_cell_number[cell_index].push_back(1);
        clone_max_driver[cell_index].push_back(0);
        clone_ave_driver[cell_index].push_back(0.0);
    }

}

#endif