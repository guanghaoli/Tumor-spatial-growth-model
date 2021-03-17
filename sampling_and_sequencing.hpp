//
// sampling_and_sequencing.hpp
//
//  Created by Guanghao Li on 2019/2/01.
//  Copyright © 2019 Guanghao Li. All rights reserved.
//
#ifndef sampling_and_sequencing_hpp
#define sampling_and_sequencing_hpp
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <pngwriter.h>
#include <blitz/blitz.h>
#include <blitz/array.h>

using namespace blitz;
using namespace std;

int sign(int x,int y){
    if(x>y){
        return 1;
    }else{
        return -1;
    }
}

void sampling_and_sequencing(int Visual_range_x,int Visual_range_y, Array<double,3> &Visual_range,Array<unsigned int,2> &Mut_trace, Array<double,2> &cell_array,double initial_b,int h,Array<double,2> colorClone,map<int,int> map_color,int split_number){
    
    Range all = Range::all();
    vector<std::uint32_t> random_data(624);
    random_device r;
    generate(random_data.begin(),random_data.end(),ref(r));
    seed_seq seed(random_data.begin(),random_data.end());
    mt19937 RNG(seed);
        
    int total_range = 15;
    int sample_range = (total_range-1)/2;
    int RO = cell_array.rows();
    int sample_number = (split_number+1)*(split_number+1);
    int time_h = h;

    map<unsigned int,int> cell_id_index;
    for(int ii=1;ii<=RO;ii++){
        cell_id_index[(unsigned int)cell_array(ii,7)] = ii;
    }

    vector<int> x_cor;
    vector<int> y_cor;
    for(int ii=1;ii<=RO;++ii){
        x_cor.push_back((int)cell_array(ii,1));
        y_cor.push_back((int)cell_array(ii,2));
    }
    vector<int>::iterator min_x = min_element(x_cor.begin(),x_cor.end());
    vector<int>::iterator max_x = max_element(x_cor.begin(),x_cor.end());
    vector<int>::iterator min_y = min_element(y_cor.begin(),y_cor.end());
    vector<int>::iterator max_y = max_element(y_cor.begin(),y_cor.end());

    int x_len = *max_x - *min_x;
    int y_len = *max_y - *min_y;
    vector<int> x_pos,y_pos;
    int center_x = (*min_x + *max_x)/2;
    int center_y = (*min_y + *max_y)/2;

    // define values
    double border = 0.05;

    if(split_number==2){

        int x_small = *min_x + int(x_len*border);
        int y_small = *min_y + int(y_len*border);
        int x_mid = (*min_x+*max_x)/2;
        int y_mid = (*min_y+*max_y)/2;
        int x_max = *max_x - int(x_len*border);
        int y_max = *max_y - int(y_len*border);

        x_pos.push_back((x_small+x_mid)/2-x_len*0.08);
        y_pos.push_back((y_mid+y_small)/2-y_len*0.08);
        x_pos.push_back(x_small);
        y_pos.push_back(y_mid);
        x_pos.push_back((x_small+x_mid)/2-x_len*0.08);
        y_pos.push_back((y_mid+y_max)/2+y_len*0.08);

        x_pos.push_back(x_mid);
        y_pos.push_back(y_small);
        x_pos.push_back(x_mid);
        y_pos.push_back(y_mid);
        x_pos.push_back(x_mid);
        y_pos.push_back(y_max);

        x_pos.push_back((x_mid+x_max)/2+x_len*0.08);
        y_pos.push_back((y_small+y_mid)/2-y_len*0.08);
        x_pos.push_back(x_max);
        y_pos.push_back(y_mid);
        x_pos.push_back((x_mid+x_max)/2+x_len*0.08);
        y_pos.push_back((y_max+y_mid)/2+y_len*0.08);
    }
    
    if(split_number==3){
 
        double x_y_inter = 0.9/3;
        int x0 = *min_x + x_len*border;
        int y0 = *min_y + y_len*border;
        int x1 = x0 + x_len*x_y_inter;
        int y1 = y0 + y_len*x_y_inter;
        int x2 = x0 + 2*x_len*x_y_inter;
        int y2 = y0 + 2*y_len*x_y_inter;
        int x3 = x0 + 3*x_len*x_y_inter;
        int y3 = y0 + 3*y_len*x_y_inter;

        x_pos.push_back((x0+x1)/2-x_len*0.04);
        y_pos.push_back((y1+y0)/2-y_len*0.04);
        x_pos.push_back(x0);
        y_pos.push_back(y1);
        x_pos.push_back(x0);
        y_pos.push_back(y2);
        x_pos.push_back((x0+x1)/2-x_len*0.04);
        y_pos.push_back((y2+y3)/2+y_len*0.04);

        x_pos.push_back(x1);
        y_pos.push_back(y0);
        x_pos.push_back(x1);
        y_pos.push_back(y1);
        x_pos.push_back(x1);
        y_pos.push_back(y2);
        x_pos.push_back(x1);
        y_pos.push_back(y3);

        x_pos.push_back(x2);
        y_pos.push_back(y0);
        x_pos.push_back(x2);
        y_pos.push_back(y1);
        x_pos.push_back(x2);
        y_pos.push_back(y2);
        x_pos.push_back(x2);
        y_pos.push_back(y3);

        x_pos.push_back((x2+x3)/2+x_len*0.04);
        y_pos.push_back((y0+y1)/2-y_len*0.04);
        x_pos.push_back(x3);
        y_pos.push_back(y1);
        x_pos.push_back(x3);
        y_pos.push_back(y2);
        x_pos.push_back((x2+x3)/2+x_len*0.04);
        y_pos.push_back((y2+y3)/2+y_len*0.04);
        
    }
    
    vector<unsigned int> region_cell_id;
    vector<unsigned int> out_cell_id,out_cell_depth;
    Array<double,3> sub_visual(total_range,total_range,4,FortranArray<3>());
    sub_visual = 0;
    int curr_x,curr_y;
    map<unsigned int,int> id_mutation; 
    map<unsigned int,int> id_count; 
    map<unsigned int,int> id_count_depth; 
    double sum_driver = 0; 
    double sum_generation = 0; 
    int cell_array_index;
    map<unsigned int,int>::iterator map_iter;
    int region_cell_number;
    unsigned int curr_id,curr_parent;
    int mean_depth = 100;
    int size_par = 2; 
    double p_par = size_par*1.0/(mean_depth + size_par); 
    double purity = 1.0; 
   
    negative_binomial_distribution<int> nega_dis(15.8,0.155);
    int mut_index,odd_or_even;
    int satisfied_number = 0; 
    
    char region_sampling_name [100] = {'\0'};
    sprintf(region_sampling_name, "./growth_%.2f_sampling/Mutation_id_h_%d_split_%d.txt",initial_b,time_h,split_number);
    FILE *region_mutation_file;
    region_mutation_file = fopen(region_sampling_name,"w+");

    char mutation_depth_name [100] = {'\0'};
    sprintf(mutation_depth_name, "./growth_%.2f_sampling/Mutation_fre_h_%d_split_%d.txt",initial_b,time_h,split_number);
    FILE *mutation_depth_file;
    mutation_depth_file = fopen(mutation_depth_name,"w+");
        
    for(int ii=0;ii<sample_number;++ii)
    {
        curr_x = x_pos[ii];
        curr_y = y_pos[ii];
        sub_visual(all,all,all) = Visual_range(Range(curr_x-sample_range,curr_x+sample_range),Range(curr_y-sample_range,curr_y+sample_range),all);
        for(int jj=1;jj<=total_range;++jj)
        {
            for(int kk=1;kk<=total_range;++kk)
            {
                if((int)sub_visual(jj,kk,1)!=0)
                {
                    region_cell_id.push_back((unsigned int)sub_visual(jj,kk,4));
                }
            }
        }

        region_cell_number = region_cell_id.size();
        if(region_cell_number>=100)
        {
            satisfied_number += 1;
            for(int ll=0;ll<region_cell_number;ll++)
            {
                cell_array_index = cell_id_index[region_cell_id[ll]];
                sum_generation += cell_array(cell_array_index,11); 
                sum_driver += cell_array(cell_array_index,6); 
                curr_parent = region_cell_id[ll];
                while(curr_parent!=0)
                {
                    mut_index = (int)((curr_parent+1)/2);
                    if(curr_parent%2==1)
                    {
                        odd_or_even = 0;
                    }else
                    {
                        odd_or_even = 3;
                    }
                    ++id_count[curr_parent]; 
                    id_mutation[curr_parent] = Mut_trace(mut_index,odd_or_even+3);
                    curr_parent = Mut_trace(mut_index,odd_or_even+2); 
                }
            }

            for(map_iter=id_count.begin();map_iter!=id_count.end();++map_iter)
            {
                int id_count_value = map_iter->second;
                if(0.5*id_count_value/region_cell_number>=0.1)
                {
                    fprintf(region_mutation_file,"%d\t",map_iter->first);
                    out_cell_id.push_back(map_iter->first);
                }
            }

            fprintf(region_mutation_file,"\n");
            for(auto &jj:out_cell_id){
                fprintf(region_mutation_file,"%d\t",id_mutation[jj]);
            }
            fprintf(region_mutation_file,"\n"); 

            char depth_afs[100] = {'\0'};
            sprintf(depth_afs,"./growth_%.2f_sampling/Sample_MAFA_%d_h_%d_split_%d.txt",initial_b,ii+1,time_h,split_number);
            FILE *depth_file;
            depth_file = fopen(depth_afs,"w+");

            for(map_iter = id_count.begin(); map_iter != id_count.end(); map_iter++)
            {
                double curr_fre = 0.5*purity*map_iter->second/region_cell_number;
                if(curr_fre > 0.001)
                { 
                    for(int mut=0; mut<id_mutation[map_iter->first]; mut++)
                    {
                        int site_depth = nega_dis(RNG); 
                        binomial_distribution<int> bino_dis(site_depth,curr_fre);
                        int var_depth = bino_dis(RNG); 
                        double now_fre = var_depth*1.0/site_depth; 
                        if(site_depth >= 10 && var_depth >= 6)
                        { 
                            fprintf(depth_file,"%d\t%d\t%f\t%f\n",map_iter->first,site_depth,curr_fre,now_fre);
                            if(now_fre >= 0.1)
                            {
                                ++id_count_depth[map_iter->first];
                            }
                        }
                    }
                }
            }
            fclose(depth_file);
            
            for(map_iter=id_count_depth.begin();map_iter!=id_count_depth.end();++map_iter)
            {
                fprintf(mutation_depth_file,"%d\t",map_iter->first);
                out_cell_depth.push_back(map_iter->first);
            }
            fprintf(mutation_depth_file,"\n");
            for(auto &jj:out_cell_depth)
            {
                fprintf(mutation_depth_file,"%d\t",id_count_depth[jj]);
            }
            fprintf(mutation_depth_file,"\n"); 
        }else
        {
            x_pos[ii] = curr_x + sign(center_x,curr_x)*x_len*0.04;
            y_pos[ii] = curr_y + sign(center_y,curr_y)*y_len*0.04;
            ii--;
            
        }

        vector<unsigned int>().swap(region_cell_id);
        vector<unsigned int>().swap(out_cell_id);
        vector<unsigned int>().swap(out_cell_depth);
        id_count.clear();
        id_mutation.clear();
        id_count_depth.clear();

    }

    fclose(region_mutation_file);
    fclose(mutation_depth_file);

    char spatial_position[100] = {'\0'};
    sprintf(spatial_position,"./growth_%.2f_sampling/Sampling_Image_%d_%d.png",initial_b,time_h,split_number);
    FILE *spatial_file;
    spatial_file = fopen(spatial_position,"wb");
    pngwriter image_spatial(Visual_range_x,Visual_range_y,0.0,spatial_position);
    char spatial_font[100] = {'\0'};
    // sprintf(spatial_font, "/Library/Fonts/Times New Roman.ttf");
    sprintf(spatial_font,"/usr/share/fonts/google-crosextra-caladea/Caladea-Regular.ttf");
    // sprintf(spatial_font,"/usr/share/fonts/liberation/LiberationSans-Regular.ttf");
    int font_size = 30;

    int cell_index_i;
    for(int ii=1;ii<=RO;ii++)
    {
        int x = (int)cell_array(ii,1);
        int y = (int)cell_array(ii,2);

        cell_index_i = (int)cell_array(ii,5);
        image_spatial.plot(x, y, colorClone(map_color[cell_index_i],1), colorClone(map_color[cell_index_i],2), colorClone(map_color[cell_index_i],3));
    }

    for(int jj=0;jj<sample_number;++jj)
    {
        char sample_label[100] = {'\0'};
        sprintf(sample_label,"X%d",jj+1);
        curr_x = x_pos[jj];
        curr_y = y_pos[jj];
        for(int xx=curr_x-sample_range;xx<=curr_x+sample_range;++xx){
            for(int yy=curr_y-sample_range;yy<=curr_y+sample_range;++yy){
                // image_spatial.plot(xx,yy,0.93,0.42,0.31);//coral2
                image_spatial.plot(xx,yy,1.0,0.0,0.0);//coral2
            }
        }
        image_spatial.plot_text(spatial_font,font_size,curr_x,curr_y,0.0,sample_label,1.0,0.0,0.0);
    }
    image_spatial.close();
    fclose(spatial_file);
    
}


#endif