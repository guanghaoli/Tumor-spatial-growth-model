#ifndef stat_map_vector_hpp
#define stat_map_vector_hpp

#include <stdio.h>
#include <map>
#include <blitz/blitz.h>
#include <blitz/array.h>

using namespace blitz;

void stat_map_vector(Array<double,2> &cell_array,map<int,vector<int> > &clone_cell_number,map<int,vector<int> > &clone_max_driver,map<int,vector<double> > &clone_ave_driver){
    
 
    int vec_size;
    map<int,vector<int> >::iterator map_vec_it;
    for(map_vec_it=clone_cell_number.begin(); map_vec_it!=clone_cell_number.end(); map_vec_it++)
    {
        clone_cell_number[map_vec_it->first].push_back(0);
        vec_size = clone_cell_number[map_vec_it->first].size();
    }
    
    for(map_vec_it=clone_max_driver.begin(); map_vec_it!=clone_max_driver.end(); map_vec_it++)
    {
        clone_max_driver[map_vec_it->first].push_back(0);
    }

    map<int,vector<double> >::iterator map_vec;
    for(map_vec=clone_ave_driver.begin(); map_vec!=clone_ave_driver.end(); map_vec++)
    {
        clone_ave_driver[map_vec->first].push_back(0.0);
    }
    
    int ROW = cell_array.rows();
    for(int ii=1; ii<=ROW; ii++)
    {
        int cell_index = (int)cell_array(ii,5);
        int driver_num = (int)cell_array(ii,6);
        clone_cell_number[cell_index][vec_size-1]++;
        clone_ave_driver[cell_index][vec_size-1] += driver_num;
        if(driver_num > clone_max_driver[cell_index][vec_size-1])
        {
            clone_max_driver[cell_index][vec_size-1] = driver_num;
        }
    }
    for(map_vec=clone_ave_driver.begin(); map_vec!=clone_ave_driver.end(); map_vec++)
    {
        int cell_index = map_vec->first;
        if(clone_cell_number[cell_index][vec_size-1] > 0)
        {
            clone_ave_driver[cell_index][vec_size-1] /= clone_cell_number[cell_index][vec_size-1];
        }
    }

}

#endif
