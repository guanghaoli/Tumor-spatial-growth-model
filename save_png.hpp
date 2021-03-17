//
//  save_data_2019.hpp
//  Cancer_simu_No_R
//
//  Created by Guanghao Li on 2019/1/11.
//  Copyright © 2019 Guanghao Li. All rights reserved.
//

#ifndef save_png_hpp
#define save_png_hpp

#include <stdio.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <pngwriter.h>
#include <map>

using namespace blitz;
void save_png(int N0,int h,double initial_b, int Visual_range_x, int Visual_range_y, Array<double,3> &Visual_range, Array<double,2> &cell_array, Array<double,2> &colorspace,Array<unsigned int,2> &Mut_trace,Array<double,2> colorClone,map<int,int> map_color,Array<double,2> &color_growth)
{

    int color_growth_number = 1000;
    int cell_index_i;
    double max_growth = 1.0;
    double min_growth = initial_b;
    double interstep = (max_growth-min_growth)/color_growth_number;
    double growth_i;
    int growth_color;

    char filedir_clone [100] = {'\0'};
    sprintf(filedir_clone, "./growth_%.2f_clonepics/%.1d.png",initial_b,h);
    char filedir_clone2 [100] = {'\0'};
    sprintf(filedir_clone2, "%.03d h",h);
    char filedir_clone3 [100] = {'\0'};
    //sprintf(filedir_clone3, "/Library/Fonts/Times New Roman.ttf");
    sprintf(filedir_clone3, "/usr/share/fonts/google-crosextra-caladea/Caladea-Regular.ttf");
    //sprintf(filedir_clone3,"/usr/share/fonts/liberation/LiberationSans-Regular.ttf");
    FILE *fid_clone;
    fid_clone=fopen(filedir_clone,"wb");
    pngwriter image_clone(Visual_range_x, Visual_range_y, 0.0, filedir_clone);
   
    char filedir_growth[100] = {'\0'};
    sprintf(filedir_growth, "./growth_%.2f_growth/%.1d.png",initial_b,h);
    char filedir_growth2[100] = {'\0'};
    sprintf(filedir_growth2, "%.03d h",h);
    char filedir_growth3[100] = {'\0'};
    //sprintf(filedir_growth3, "/Library/Fonts/Times New Roman.ttf");
    sprintf(filedir_growth3, "/usr/share/fonts/google-crosextra-caladea/Caladea-Regular.ttf");
    //sprintf(filedir_growth3,"/usr/share/fonts/liberation/LiberationSans-Regular.ttf");
    FILE *fid_growth;
    fid_growth=fopen(filedir_growth,"wb");
    pngwriter image_growth(Visual_range_x, Visual_range_y, 0.0, filedir_growth);

    int C0 = cell_array.rows();
    for(int ii=1;ii<=C0;ii++)
    {
        int x = cell_array(ii,1);
        int y = cell_array(ii,2);

        /*clone*/
        cell_index_i = (int)cell_array(ii,5);
        image_clone.plot(x, y, colorClone(map_color[cell_index_i],1), colorClone(map_color[cell_index_i],2), colorClone(map_color[cell_index_i],3));

        /*growth*/
        growth_i = cell_array(ii,10);
        growth_color = (int)((growth_i-min_growth)/interstep);
        growth_color = growth_color>color_growth_number?color_growth_number:growth_color;
        growth_color = growth_color<1?1:growth_color;
        image_growth.plot(x, y, color_growth(growth_color,2),color_growth(growth_color,3), 0.0);

    }
    int text_xpos = Visual_range_x-200;
    int text_ypos = Visual_range_y-100;
    int font_size = 30;

    image_clone.plot_text(filedir_clone3, font_size, text_xpos, text_ypos, 0.0, filedir_clone2, 1.0, 1.0, 1.0);
    image_clone.close();
    fclose(fid_clone);

    image_growth.plot_text(filedir_growth3, font_size, text_xpos, text_ypos, 0.0, filedir_growth2, 1.0, 1.0, 1.0);
    image_growth.close();
    fclose(fid_growth);


}

#endif /* save_png_hpp */
