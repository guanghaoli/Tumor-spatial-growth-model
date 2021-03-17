#ifndef color_clone_hpp
#define color_clone_hpp

#include <blitz/blitz.h>
#include <blitz/array.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
using namespace blitz;
using namespace std;

void color_clone(Array<double,2> &cell_array,Array<double,2> &colorClone,map<int,int> &map_color,double initial_b){
    
    vector<std::uint32_t> random_data(624);
    random_device r;
    generate(random_data.begin(),random_data.end(),ref(r));
    seed_seq sed(random_data.begin(),random_data.end());
    mt19937 RNG(sed);
    uniform_real_distribution<> dis(0.05,0.95);
    
    string sigma_g = to_string(initial_b);
    sigma_g = sigma_g.substr(0,sigma_g.size()-4);// keep two number 

    int RO = cell_array.rows();
    int cell_index;
    //map_color using map to give the color number for cell id in this time, cell_id: color_number
    for(int i=1;i<=RO;++i){
        cell_index = (int)cell_array(i,5); // The cell index generate the clone colors
        map_color[cell_index] = i;
    }
    colorClone.resize(RO,3);
    for(int j=1;j<=RO;++j){
        colorClone(j,1) = dis(RNG);
        colorClone(j,2) = dis(RNG);
        colorClone(j,3) = dis(RNG);
    }

    /*Output colorclone*/
    string filename_color = "./growth_"+ sigma_g +"_sampling/" + "colorClone.txt";
    ofstream ofs_color(filename_color);
    if(ofs_color.bad()){
        cerr << "Unable to write to file: "<<filename_color<<endl;
        exit(1);
    }
    ofs_color << colorClone<<endl;
    ofs_color.close();

    /*output colorclone_new.txt*/
    string filename_color_new = "./growth_"+ sigma_g +"_sampling/" + "colorClone_new.txt";
    ofstream ofs_color_new(filename_color_new);
    if(ofs_color_new.bad()){
        cerr << "Unable to write to file: "<<filename_color_new<<endl;
        exit(1);
    }
    for(int ii=1;ii<=RO;ii++){
        ofs_color_new << colorClone(ii,1) <<"\t"<< colorClone(ii,2) <<"\t"<< colorClone(ii,3) <<"\n";
    }
    ofs_color_new.close();


    /*Output map_color*/
    string filename_mapcolor = "./growth_"+ sigma_g +"_sampling/" + "mapcolor.txt";
    ofstream ofs_mapcolor(filename_mapcolor);
    if(ofs_mapcolor.bad()){
        cerr << "Unable to write to file:" << filename_mapcolor <<endl;
        exit(1);
    }
    map<int,int>::iterator map_iter;
    for(map_iter=map_color.begin();map_iter!=map_color.end();map_iter++){
        ofs_mapcolor << map_iter->first << "\t" << map_iter->second << "\n";
    }
    ofs_mapcolor.close();

}

#endif
