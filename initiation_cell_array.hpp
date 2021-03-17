#ifndef initiation_cell_array_hpp
#define initiation_cell_array_hpp

#include <stdio.h>
#include <random>
#include <cmath>
#include <algorithm>
#include <functional>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cstring>
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
using namespace blitz;
using namespace std;


void initiation_cell_array(Array<double,2> &cell_array,Array<double,3> &Visual_range,int N0,int Visual_range_x,int Visual_range_y,unsigned int &cell_label,double initial_b,double in_birth,double quiescent)
{
    
    /*The coordinate*/
    Array<int,2> Visual_coor(2,N0,FortranArray<2>());
    int number_cor=1;
    for(int x=1;x<=Visual_range_x;x++){
        for(int y=1;y<=Visual_range_y;y++){
            if(Visual_range(x,y,1)==1){
                Visual_coor(1,number_cor) = x;
                Visual_coor(2,number_cor) = y;
                number_cor += 1;
            }
        }
        
    }
    
    for (int x=1;x<=N0;x++){
        int current_x = Visual_coor(1,x);
        int current_y = Visual_coor(2,x);
        cell_array(x,1) = current_x; 
        cell_array(x,2) = current_y; 

        cell_array(x,5) = x; 

        cell_array(x,7) = x; 
        cell_array(x,9) = 1; 
        cell_array(x,10) = initial_b;
        cell_array(x,11) = 0; 
        cell_array(x,12) = in_birth;
        cell_array(x,13) = quiescent;

        Visual_range(current_x,current_y,2) = x; 
        Visual_range(current_x,current_y,3) = 0; 
        Visual_range(current_x,current_y,4) = x; 

        cell_label += 1;

    }

}

#endif
