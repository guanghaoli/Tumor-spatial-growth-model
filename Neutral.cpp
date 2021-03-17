//
//  main_new.cpp
//  Cancer_20190124_speed_no_float
//
//  Created by Guanghao Li on 2019/1/24.
//  Copyright © 2019 Guanghao Li. All rights reserved.
//

#include <iostream>
#include <time.h>
#include <memory>
#include <stdio.h>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <map>
#include <set>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <random>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_block.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_matrix.h>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include "initiation_visual_range.hpp"
#include "initiation_cell_array.hpp"
#include "color_clone.hpp"
#include "save_png.hpp"
#include "save_mut_trace.hpp"
#include "save_cell_array.hpp"
#include "judge_death.hpp"
#include "sampling_and_sequencing.hpp"
#include "division_invasion_slow.hpp"
#include "initiation_map_vector.hpp"
#include "stat_map_vector.hpp"
#include "save_clone_related.hpp"

using namespace blitz;
using namespace std;


int main(int argc,char *argv[]){


    int compe_cutoff = atoi(argv[1]);
    int mutation_rate = atoi(argv[2]);
    double in_birth = atof(argv[3]);
    double quiescent = atof(argv[4]);
    double initial_b = atof(argv[5]);
    double inside_pro = atof(argv[6]);


    int Visual_range_x = 2000;
    int Visual_range_y = 2000;
    int R0 = 1;
    int N0 = 0;

    Range all = Range::all();
 
    int cell_column = 13;
    int mut_column = 6;
    Array<double,3> Visual_range(Visual_range_x,Visual_range_y,4,FortranArray<3>());
    Visual_range(all,all,all) = 0;

    Array<double,2> cell_array(1,cell_column,FortranArray<2>());
    cell_array = 0;

    Array<unsigned int,2> Mut_trace(1,mut_column,FortranArray<2>());
    Mut_trace = 0;
    Array<unsigned int,2> Mut_temp(1,mut_column,FortranArray<2>());
    Mut_temp = 0;
    int Mut_index = 1;

    unsigned int cell_label = 1;

    //double genome_rate = 1.0*pow(10,-8); // mutation rate for whole genome
    //double Driver_prob = pow(10,-6); // driver mutation rate
    double Driver_prob = 0.0; // driver mutation rate

    double Delete_prob = 5.0*pow(10, -3); // 5*10^6 deleterious sites in genome
    int pos_lambda = mutation_rate;

    double initial_d = 1-initial_b;
    double s_coef = 0.1;

    int Time_generation = 360000;

    int GUARD = 0;

    
    const gsl_rng_type *T00;
    gsl_rng *r00;
    gsl_rng_env_setup();
    T00 = gsl_rng_ranlxs0;
    gsl_rng_default_seed = ((unsigned long)(time(NULL)));
    r00 = gsl_rng_alloc(T00);
    vector<std::uint32_t> random_data(624);
    random_device r;
    generate(random_data.begin(),random_data.end(),ref(r));
    seed_seq sed(random_data.begin(),random_data.end());
    mt19937 RNG(sed);
    uniform_real_distribution<> dis(0.05,0.95);


    int color_number = 4 * R0 * R0;
    Array<double,2> colorspace(color_number,4,FortranArray<2>());
    colorspace = 0;
    for (int i=1;i<=color_number;i++){
        for (int j=1;j<=4;j++){
            if(j==1){
                colorspace(i,j) = 1;
            }
            else{
                colorspace(i,j) = dis(RNG);
            }
        }
    }

    double interstep=0.001;
    int growth_number = (int)(1.0/interstep);
    Array<double,2> color_growth(growth_number,3,FortranArray<2>());
    color_growth = 0;
    for(int i=1;i<=growth_number;i++){
        color_growth(i,1) = 1;
        color_growth(i,2) = i*interstep;
        color_growth(i,3) = 1-i*interstep;
    }

    initiation_visual_range(Visual_range,Visual_range_x,Visual_range_y,N0,R0);

    cell_array.resize(N0,cell_column);
    cell_array = 0;
    initiation_cell_array(cell_array,Visual_range,N0,Visual_range_x,Visual_range_y,cell_label,initial_b,in_birth,quiescent);

    Mut_trace.resize(800000000,mut_column);
    Mut_trace = 0;
    for(int i=1;i<=N0;++i)
    {
        if(i%2==1)
        {
            int j = (i+1)/2;
            Mut_trace(j,1) = (unsigned int)cell_array(i,7);
            Mut_trace(j,2) = 0; 
            Mut_trace(j,3) = 0; 
        }else
        {
            int j = i/2;
            Mut_trace(j,4) = (unsigned int)cell_array(i,7);
            Mut_trace(j,5) = 0; 
            Mut_trace(j,6) = 0; 
        }
    }

    Array<double,2> colorClone(1,3,FortranArray<2>());
    colorClone = 0;
    map<int,int> map_color;
    for(int i=1;i<=N0;++i)
    {
        colorClone(i,1) = dis(RNG);
        colorClone(i,2) = dis(RNG);
        colorClone(i,3) = dis(RNG);
        map_color[i] = i;
    }

    /*Parameter Output*/
    /**************************/
    char dirname [100] = {'\0'};
    sprintf(dirname, "mkdir ./growth_%.2f",initial_b);
    system(dirname);
    char dirname1 [100] = {'\0'};
    sprintf(dirname1, "mkdir ./growth_%.2f_growth",initial_b);
    system(dirname1);
    char dirname2 [100] = {'\0'};
    sprintf(dirname2, "mkdir ./growth_%.2f_clonepics",initial_b);
    system(dirname2);
    char dirname3[100] = {'\0'};
    sprintf(dirname3, "mkdir ./growth_%.2f_sampling",initial_b);
    system(dirname3);

    char filedir [100] = {'\0'};
    sprintf(filedir, "./Parameters.txt");
    FILE * fid1;
    fid1=fopen (filedir,"w+");
    fprintf(fid1, "%s %s %d\n" ,"R0", "=", R0);
    fprintf(fid1, "%s %s %d\n" ,"Visual_range_x", "=", Visual_range_x);
    fprintf(fid1, "%s %s %d\n" ,"Visual_range_y", "=", Visual_range_y);
    fprintf(fid1, "%s %s %f\n","initial_birth_rate","=",initial_b);
    fprintf(fid1, "%s %s %f\n","initial_inner_birth_rate","=",in_birth);
    fprintf(fid1, "%s %s %d\n","Mutation rate","=",mutation_rate);
    fprintf(fid1, "%s %s %f\n","Driver_prob","=",Driver_prob);
    fprintf(fid1, "%s %s %f\n","Delete_prob","=",Delete_prob);
    fprintf(fid1, "%s %s %d\n","pos_lambda","=",pos_lambda);
    fprintf(fid1, "%s %s %f\n","s_coef","=",s_coef);
    fprintf(fid1, "%s %s %d\n","compe_cutoff","=",compe_cutoff);
    fclose(fid1);

    int h = 0;

    Array<double, 2> cell_array_temp(1,cell_column,FortranArray<2>());
    cell_array_temp = 0;


    int h_out=0;
    // save_png(N0,h,initial_b,Visual_range_x,Visual_range_y,Visual_range,cell_array,colorspace,Mut_trace,colorClone,map_color,color_growth);
    save_cell_array(initial_b,h,cell_array);
    h_out += 1;

    vector<int> h_record,cell_number,clone_number;
    set<int> clone_name;

    map<int,vector<int> > clone_cell_number,clone_max_driver;
    map<int,vector<double> > clone_ave_driver;

    int out_index=0;
    int split_number;
    while(h < Time_generation)
    {

        judge_death(Visual_range_x,Visual_range_y,cell_array,Visual_range,cell_array_temp,cell_column);
        if(h>60)
        {
            stat_map_vector(cell_array,clone_cell_number,clone_max_driver,clone_ave_driver);
        }
        time_t rawtime;
        struct tm *timeinfo;
        time(&rawtime);
        timeinfo = localtime(&rawtime);
        cout << "h = " << h << "  time: " << asctime(timeinfo) << endl;
        int C1 = cell_array.rows();
        cout <<"C1: "<< C1<<endl;
        if(C1==0){
            cout << "Another random number !"<<endl;
            exit(1);
        }

        if(h==h_out){
           
            if(h_out==60)
            {
                for(int i=1;i<=C1;++i)
                {
                    cell_array(i,5) = cell_array(i,7); 
                }
                color_clone(cell_array,colorClone,map_color,initial_b);
                initiation_map_vector(cell_array,clone_cell_number,clone_max_driver,clone_ave_driver);
            }
            // save_png(N0,h,initial_b,Visual_range_x,Visual_range_y,Visual_range,cell_array,colorspace,Mut_trace,colorClone,map_color,color_growth);
      
            h_out += 1;
        }

        if(out_index<1){
            if(C1>500000){
                out_index += 1;
                cout<<"The size is Enough!"<<endl;
                cout<<"Cell_label:"<<cell_label<<endl;
                split_number = 3;
                // save_png(N0,h,initial_b,Visual_range_x,Visual_range_y,Visual_range,cell_array,colorspace,Mut_trace,colorClone,map_color,color_growth);
                save_cell_array(initial_b,h,cell_array);
                sampling_and_sequencing(Visual_range_x,Visual_range_y,Visual_range,Mut_trace,cell_array,initial_b,h,colorClone,map_color,split_number);
                save_clone_related(h_record,cell_number,clone_number,cell_label,h,mut_column,initial_b,clone_cell_number,clone_max_driver,clone_ave_driver,split_number);
                // save_mut_trace(Mut_trace,cell_array,h_record,cell_number,clone_number,Mut_temp,cell_label,h,mut_column,initial_b);
                // return 0;
            }
        }

        if(out_index<2){
            if(C1>1000000){
                out_index += 1;
                cout<<"The size is Enough!"<<endl;
                cout<<"Cell_label:"<<cell_label<<endl;
                split_number = 3;
                // save_png(N0,h,initial_b,Visual_range_x,Visual_range_y,Visual_range,cell_array,colorspace,Mut_trace,colorClone,map_color,color_growth);
                save_cell_array(initial_b,h,cell_array);
                sampling_and_sequencing(Visual_range_x,Visual_range_y,Visual_range,Mut_trace,cell_array,initial_b,h,colorClone,map_color,split_number);
                save_clone_related(h_record,cell_number,clone_number,cell_label,h,mut_column,initial_b,clone_cell_number,clone_max_driver,clone_ave_driver,split_number);
                // save_mut_trace(Mut_trace,cell_array,h_record,cell_number,clone_number,Mut_temp,cell_label,h,mut_column,initial_b);
                // return 0;
            }
        }

        if(out_index<3){
            if(C1>1500000){
                out_index += 1;
                cout<<"The size is Enough!"<<endl;
                cout<<"Cell_label:"<<cell_label<<endl;
                split_number = 3;
                // save_png(N0,h,initial_b,Visual_range_x,Visual_range_y,Visual_range,cell_array,colorspace,Mut_trace,colorClone,map_color,color_growth);
                save_cell_array(initial_b,h,cell_array);
                sampling_and_sequencing(Visual_range_x,Visual_range_y,Visual_range,Mut_trace,cell_array,initial_b,h,colorClone,map_color,split_number);
                save_clone_related(h_record,cell_number,clone_number,cell_label,h,mut_column,initial_b,clone_cell_number,clone_max_driver,clone_ave_driver,split_number);
                // save_mut_trace(Mut_trace,cell_array,h_record,cell_number,clone_number,Mut_temp,cell_label,h,mut_column,initial_b);
                // return 0;
            }
        }

        if(out_index<4){
            if(C1>2000000){
                out_index += 1;
                cout<<"The size is Enough!"<<endl;
                cout<<"Cell_label:"<<cell_label<<endl;
                split_number = 3;
                // save_png(N0,h,initial_b,Visual_range_x,Visual_range_y,Visual_range,cell_array,colorspace,Mut_trace,colorClone,map_color,color_growth);
                save_cell_array(initial_b,h,cell_array);
                sampling_and_sequencing(Visual_range_x,Visual_range_y,Visual_range,Mut_trace,cell_array,initial_b,h,colorClone,map_color,split_number);
                save_clone_related(h_record,cell_number,clone_number,cell_label,h,mut_column,initial_b,clone_cell_number,clone_max_driver,clone_ave_driver,split_number);
                // save_mut_trace(Mut_trace,cell_array,h_record,cell_number,clone_number,Mut_temp,cell_label,h,mut_column,initial_b);
                return 0;
            }
        }


        cell_array_temp.resize(C1, cell_column);
        int division_time = 0;
        division_invasion_slow(h,cell_array,Visual_range,cell_label,Mut_trace,cell_array_temp,division_time,clone_name,Driver_prob,Delete_prob,pos_lambda,initial_b,s_coef,compe_cutoff,quiescent,in_birth,GUARD,inside_pro);
        h_record.push_back(h);
        cell_number.push_back(C1);
        clone_number.push_back(clone_name.size());
        clone_name.clear();

        h+=1;

    }

    return 0;

}

