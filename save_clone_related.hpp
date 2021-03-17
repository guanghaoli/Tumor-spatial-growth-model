//
//  save_last.hpp
//
//  Created by Guanghao Li on 2019/2/01.
//  Copyright © 2019 Guanghao Li. All rights reserved.
//
#ifndef save_clone_related_hpp
#define save_clone_related_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <map>
#include <blitz/blitz.h>
#include <blitz/array.h>

using namespace blitz;
using namespace std;

void save_clone_related(vector<int> &h_record,vector<int> &cell_number,vector<int> &clone_number,unsigned int cell_label,int h,double mut_column,double initial_b,map<int,vector<int> > &clone_cell_number,map<int,vector<int> > &clone_max_driver,map<int,vector<double> > &clone_ave_driver,int split_number){

    string time_h = to_string(h);
    string sigma_g = to_string(initial_b);
    sigma_g = sigma_g.substr(0,sigma_g.size()-4);
    string split_str = to_string(split_number);

    string filename_cell_clone = "./growth_"+ sigma_g +"_sampling/" + "Clone_cell_number_"+ time_h + "_split_"+ split_str + ".txt";
    cout << "clone_stat: " << filename_cell_clone<<endl;
    ofstream ofs_cell_clone(filename_cell_clone);
    if(ofs_cell_clone.bad()){
        cerr << "Unable to write to file: "<< filename_cell_clone <<endl;
        exit(1);
    }
    int h0 = 60;
    ofs_cell_clone << "time";
    map<int,vector<int> >::iterator map_vec_it=clone_cell_number.begin();
    int cell_index_begin = map_vec_it->first;
    int cell_index_size = clone_cell_number[cell_index_begin].size();
    for(map_vec_it=clone_cell_number.begin(); map_vec_it!=clone_cell_number.end(); map_vec_it++){
        int cell_index = map_vec_it->first;
        ofs_cell_clone << "\t" << cell_index;
    }
    ofs_cell_clone << "\n";
    for(int i=0;i<cell_index_size;i++){
        ofs_cell_clone << h0;
        h0++;
        for(map_vec_it=clone_cell_number.begin(); map_vec_it!=clone_cell_number.end(); map_vec_it++){
            int cell_index = map_vec_it->first;
            ofs_cell_clone << "\t" <<  clone_cell_number[cell_index][i];
        }
        ofs_cell_clone << "\n";
    }

    string filename_max_driver = "./growth_"+ sigma_g +"_sampling/" + "Clone_max_driver_"+ time_h + "_split_"+ split_str +".txt";
    cout << "max_driver_stat: " << filename_max_driver<<endl;
    ofstream ofs_max_driver(filename_max_driver);
    if(ofs_max_driver.bad()){
        cerr << "Unable to write to file: "<< filename_max_driver <<endl;
        exit(1);
    }
    h0 = 60;
    ofs_max_driver << "time";
    map_vec_it=clone_max_driver.begin();
    cell_index_begin = map_vec_it->first;
    cell_index_size = clone_max_driver[cell_index_begin].size();
    for(map_vec_it=clone_max_driver.begin(); map_vec_it!=clone_max_driver.end(); map_vec_it++){
        int cell_index = map_vec_it->first;
        ofs_max_driver << "\t" << cell_index;
    }
    ofs_max_driver << "\n";
    for(int i=0;i<cell_index_size;i++){
        ofs_max_driver << h0;
        h0++;
        for(map_vec_it=clone_max_driver.begin(); map_vec_it!=clone_max_driver.end(); map_vec_it++){
            int cell_index = map_vec_it->first;
            ofs_max_driver << "\t" <<  clone_max_driver[cell_index][i];
        }
        ofs_max_driver << "\n";
    }

    string filename_ave_driver = "./growth_"+ sigma_g +"_sampling/" + "Clone_ave_driver_"+ time_h + "_split_"+ split_str +".txt";
    cout << "max_driver_stat: " << filename_ave_driver<<endl;
    ofstream ofs_ave_driver(filename_ave_driver);
    if(ofs_ave_driver.bad()){
        cerr << "Unable to write to file: "<< filename_ave_driver <<endl;
        exit(1);
    }
    h0 = 60;
    ofs_ave_driver << "time";
    map<int,vector<double> >::iterator map_vec_it_d=clone_ave_driver.begin();
    cell_index_begin = map_vec_it_d->first;
    cell_index_size = clone_ave_driver[cell_index_begin].size();
    for(map_vec_it_d=clone_ave_driver.begin(); map_vec_it_d!=clone_ave_driver.end(); map_vec_it_d++){
        int cell_index = map_vec_it_d->first;
        ofs_ave_driver << "\t" << cell_index;
    }
    ofs_ave_driver << "\n";
    for(int i=0;i<cell_index_size;i++){
        ofs_ave_driver << h0;
        h0++;
        for(map_vec_it_d=clone_ave_driver.begin(); map_vec_it_d!=clone_ave_driver.end(); map_vec_it_d++){
            int cell_index = map_vec_it_d->first;
            ofs_ave_driver << "\t" <<  clone_ave_driver[cell_index][i];
        }
        ofs_ave_driver << "\n";
    }

    string filename_clone = "./growth_"+ sigma_g +"_sampling/" + "Cell_number_clone_"+ time_h + "_split_"+ split_str +".txt";
    cout << "clone: " << filename_clone <<endl;
    ofstream ofs_clone(filename_clone);
    if(ofs_clone.bad()){
        cerr << "Unable to write to file:" << filename_clone <<endl;
        exit(1);
    }

    int h_number = h_record.size();
    for(int ii=0;ii<h_number;++ii){
        ofs_clone << h_record[ii] <<"\t"<< cell_number[ii] <<"\t"<< clone_number[ii] <<"\n";
    }
    ofs_clone.close();

}


#endif
