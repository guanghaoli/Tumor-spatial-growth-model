//
//  save_cell_array.hpp
//  Cancer_simu_No_R
//
//  Created by Guanghao Li on 2019/1/11.
//  Copyright © 2019 Guanghao Li. All rights reserved.
//

#ifndef save_cell_array_hpp
#define save_cell_array_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <map>

using namespace blitz;
void save_cell_array(double initial_b,int h, Array<double,2> &cell_array)
{
    char filedir_cell [100] = {'\0'};
    sprintf(filedir_cell, "./growth_%.2f/Cell_array_growth_%.2f_h_%.1d.txt",initial_b,initial_b,h);
    FILE * fid_cell;
    fid_cell=fopen (filedir_cell,"w+");

    int RO= cell_array.rows();
    int CO = cell_array.columns();
    for (int ro=1;ro<=RO;ro++)
    {
        for(int co=1;co<=CO;co++)
        {
            if(co<CO)
            {
                fprintf(fid_cell,"%lf\t",cell_array(ro,co));
            }
            else
            {
                fprintf(fid_cell,"%lf\n",cell_array(ro,co));
            }
        }
    }
    fclose(fid_cell);

}

#endif /* save_cell_array_hpp */
