//
//  judge_death.hpp
//  Cancer_simu_20181121
//
//  Created by Guanghao Li on 2018/11/30.
//  Copyright © 2018 Guanghao Li. All rights reserved.
//

#ifndef judge_death_hpp
#define judge_death_hpp

#include <stdio.h>
#include <vector>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <algorithm>

using namespace std;
using namespace blitz;
void judge_death(int Visual_range_x, int Visual_range_y,Array<double, 2> &cell_array,Array<double,3> &Visual_range,Array<double,2> cell_array_temp,int cell_column)
{
    Range all = Range::all();
    int C0= cell_array.rows();
    int alive_cell = 0;
    for(int i=1;i<=C0;i++)
    {
        if((int)cell_array(i,9)!=0){
            alive_cell += 1;
        }else{
            Visual_range((int)cell_array(i,1),(int)cell_array(i,2),all) = 0;
        }
    }

    cell_array_temp.resize(alive_cell,cell_column);
    cell_array_temp = 0;
    int j = 1; 
    for (int i =1;i<=C0;i++){
        if(cell_array(i,9)!=0){
            cell_array_temp(j,all) = cell_array(i,all);
            j++;
        }
    }
    cell_array.resize(alive_cell,cell_column);
    cell_array=0;
    cell_array(all,all)=cell_array_temp(all,all);
    cell_array_temp.resize(1,cell_column);
}

#endif /* judge_death_hpp */
