//
//  save_last.hpp
//
//  Created by Guanghao Li on 2019/2/01.
//  Copyright © 2019 Guanghao Li. All rights reserved.
//
#ifndef save_mut_trace_hpp
#define save_mut_trace_hpp

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <blitz/blitz.h>
#include <blitz/array.h>

using namespace blitz;
using namespace std;

void save_mut_trace(Array<unsigned int,2> &Mut_trace, Array<double,2> &cell_array,vector<int> &h_record,vector<int> &cell_number,vector<int> &clone_number,Array<unsigned int,2> Mut_temp,unsigned int cell_label,double h,double mut_column,double initial_b){

    Range all = Range::all();
    int mut_temp_size = (int)((cell_label)/2);
    Mut_temp.resize(mut_temp_size,mut_column);
    Mut_temp(all,all) = Mut_trace(Range(1,mut_temp_size),all);
    string time_h = to_string((int)h);
    string sigma_g = to_string(initial_b);
    sigma_g = sigma_g.substr(0,sigma_g.size()-4);
 
    string filename = "./growth_"+ sigma_g +"_sampling/" + "Mut_trace_" + time_h +".txt";
    cout << "Mut_trace: " << filename <<endl;
    ofstream ofs(filename.c_str());
    if(ofs.bad()){
        cerr << "Unable to write to file: "<<filename<<endl;
        exit(1);
    }
    ofs << Mut_temp<<endl;
    ofs.close();

    string filename_cell = "./growth_"+ sigma_g +"_sampling/" +"cell_array_" + time_h +".txt";
    cout << "cell_array: " << filename_cell <<endl;
    ofstream ofs_cell(filename_cell.c_str());
    if(ofs_cell.bad()){
        cerr << "Unable to write to file: "<<filename_cell<<endl;
        exit(1);
    }

    ofs_cell << fixed;
    ofs_cell << cell_array<<endl;
    ofs_cell.close();

}


#endif