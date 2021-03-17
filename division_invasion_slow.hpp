//
//  division_invasion_slow.hpp
//  Cancer_20190124_speed_no_float
//
//  Created by Guanghao Li on 2019/1/24.
//  Copyright © 2019 Guanghao Li. All rights reserved.
//

#ifndef division_invasion_slow_hpp
#define division_invasion_slow_hpp

#include <stdio.h>
#include <random>
#include <algorithm>
#include <cmath>
#include <functional>
#include <vector>
#include <set>
#include <blitz/blitz.h>
#include <blitz/array.h>
#include <gsl/gsl_sf_lambert.h>
using namespace blitz;
using namespace std;

void division_invasion_slow(int h,Array<double, 2> &cell_array, Array<double, 3> &Visual_range,unsigned int &cell_label,Array<unsigned int, 2> &Mut_trace,Array<double,2> &cell_array_temp,int &division_time,set<int> &clone_name,double Driver_prob,double Delete_prob,int pos_lambda,double initial_b,double s_coef,int compe_cutoff,double quiescent,double in_birth,int &GUARD,double inside_pro)
{
    Range all = Range::all();
    vector<std::uint32_t> random_data(624);
    random_device r;
    generate(random_data.begin(),random_data.end(),ref(r));
    seed_seq seed(random_data.begin(),random_data.end());
    mt19937 RNG(seed);
    
    poisson_distribution<int> mut_pos_distribution(pos_lambda);
    poisson_distribution<int> mut_Driver_distribution(Driver_prob);
    poisson_distribution<int> mut_Delete_distribution(Delete_prob);
    uniform_real_distribution<double> unif_real_dis(0.0,1.0);
    int mut_number1;
    int mut_number_Driver1;
    int mut_number2;
    int mut_number_Driver2;
    
    int cell_column = 13;
    int mut_column_temp = 3;
    Array<double, 2> cell_array_temp2(1,cell_column,FortranArray<2>());
    cell_array_temp2 = 0;
    Array<unsigned int,2> Mut_temp1(1,mut_column_temp,FortranArray<2>());
    Array<unsigned int,2> Mut_temp2(1,mut_column_temp,FortranArray<2>());
    Array<double,2> cell_temp(1,cell_column,FortranArray<2>());
    Array<double,2> cell_copy(1,cell_column,FortranArray<2>());
    
    int division_range=3;
    int half_range = (division_range-1)/2;
    Array<double,3> sub_visual(division_range,division_range,4,FortranArray<3>());
    sub_visual = 0;
    

    int x,y,left_most_x,left_most_y;
    int number,loci_number;
    int cell_index;
    vector<int> cell_count;
 
    int C1 = cell_array.rows();
    map<unsigned int,int> cell_id_index;
    vector<int> x_cor; 
    vector<int> y_cor; 
    for(int ii=1;ii<=C1;ii++)
    {
        cell_id_index[(unsigned int)cell_array(ii,7)] = ii;
        x_cor.push_back((int)cell_array(ii,1));
        y_cor.push_back((int)cell_array(ii,2));
    }

    map<unsigned int,int>::iterator map_cell;
    map_cell = cell_id_index.end();
    map_cell--;
    unsigned int max_cell_id = map_cell->first;


    vector<int>::iterator min_x = min_element(x_cor.begin(),x_cor.end());
    vector<int>::iterator max_x = max_element(x_cor.begin(),x_cor.end());
    vector<int>::iterator min_y = min_element(y_cor.begin(),y_cor.end());
    vector<int>::iterator max_y = max_element(y_cor.begin(),y_cor.end());
    int x_len = *max_x - *min_x;
    int y_len = *max_y - *min_y;
    int center_x = (*min_x + *max_x)/2;
    int center_y = (*min_y + *max_y)/2;
    int judge_len;
    if(x_len>y_len)
    {
        judge_len = y_len/2;
    }else{
        judge_len = x_len/2;
    }

    unsigned int curr_cell_id; 
    int mut_index,odd_or_even;

    for(int i=1;i<=C1;++i)
    {
        clone_name.insert((int)cell_array(i,5));
        Mut_temp1 = 0;
        Mut_temp2 = 0;
        cell_temp = 0;
        cell_copy = 0;
        curr_cell_id = (unsigned int)cell_array(i,7);

        x=(int)cell_array(i,1);
        y=(int)cell_array(i,2);
        left_most_x = x - half_range;
        left_most_y = y - half_range;
        
        sub_visual(all,all,all) = Visual_range(Range(x-half_range,x+half_range),Range(y-half_range,y+half_range),all);
        
        cell_copy(1,all) = cell_array(i,all);
        
        if(h<60)
        {
            double divi_prob = unif_real_dis(RNG);
            if(C1<10)
            {
                while(divi_prob > cell_array(i,10))
                {
                    divi_prob = unif_real_dis(RNG);
                }
            }
            if(divi_prob <= cell_array(i,10))
            {
                number = 0;
                for(int i=1;i<=division_range;i++)
                {
                    for(int j=1;j<=division_range;j++)
                    {
                        number += 1;
                        if((int)sub_visual(i,j,1)==0)
                        {
                            cell_count.push_back(number);
                        }
                    }
                }
                loci_number = (int)cell_count.size();
                if (loci_number!=0)
                {
                    
                    mut_number1 = mut_pos_distribution(RNG);
                    mut_number_Driver1 = mut_Driver_distribution(RNG);

                    mut_number2 = mut_pos_distribution(RNG);
                    mut_number_Driver2 = mut_Driver_distribution(RNG);
             
                    cell_array(i,4) = mut_number1 + mut_number_Driver1;
                    cell_array(i,7) = cell_label;
                    cell_array(i,8) = cell_copy(1,8) + cell_array(i,4); 
                    cell_array(i,11) = cell_copy(1,11) + 1;

                    if(mut_number_Driver1 > 0)
                    {
                        cell_array(i,6) += 1;
                        cell_array(i,10) += cell_array(i,10)*s_coef;
                        if(cell_array(i,10)>1)
                        {
                            cell_array(i,10) = 1;
                        }
                    }
                    Visual_range(x,y,4) = cell_label;
                    Visual_range(x,y,3) = cell_array(i,6);
                    
                    Mut_temp1(1,1) = cell_label; 
                    Mut_temp1(1,2) = cell_copy(1,7);
                    Mut_temp1(1,3) = mut_number1 + mut_number_Driver1;
                    cell_id_index[cell_label] = i; 


                    cell_label=cell_label+1; 
                    cell_index=cell_copy(1,5);
                    shuffle(cell_count.begin(), cell_count.end(),RNG);
                    
                    int loci=cell_count[0]; 
                    int loci_x = left_most_x + (loci-half_range)/division_range;
                    int loci_y = left_most_y + (loci-half_range)%division_range;
                    cell_temp(1,1)= loci_x; 
                    cell_temp(1,2)= loci_y; 
                    cell_temp(1,4) = mut_number2 + mut_number_Driver2;
                    cell_temp(1,5) = cell_index; 
                    cell_temp(1,6) = cell_copy(1,6); 
                    cell_temp(1,7) = cell_label; 
                    cell_temp(1,8) = cell_copy(1,8) + cell_temp(1,4); 
                    cell_temp(1,9) = cell_copy(1,9); 
                    cell_temp(1,10) = cell_copy(1,10); 
                    cell_temp(1,11) = cell_copy(1,11) + 1;
                    cell_temp(1,12) = cell_copy(1,12);
                    cell_temp(1,13) = cell_copy(1,13);

                    if(mut_number_Driver2 > 0)
                    {
                        cell_temp(1,3) += 1;
                    }
                    
                
                    Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),1)=1;
                    Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),2)=cell_index;
                    Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),3)=cell_temp(1,6); 
                    Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),4)=cell_label;
                    
                    Mut_temp2(1,1) = cell_label;
                    Mut_temp2(1,2) = cell_copy(1,7); 
                    Mut_temp2(1,3) = mut_number2 + mut_number_Driver2;

               
                    division_time += 1;
                   

                    cell_array_temp(division_time,all) = cell_temp(1,all);
                    cell_id_index[cell_label] = C1 + division_time;
     
                    cell_label += 1;
                    mut_index = (int)((cell_label-2)/2);

                    Mut_trace(mut_index,Range(4,6)) = Mut_temp1(1,all);
                    Mut_trace(mut_index+1,Range(1,3)) = Mut_temp2(1,all);
          
                }
            }
            else 
            {
                cell_array(i,9) = 0;        
            }

        }
        else 
        {
           
            if(cell_array(i,6)<compe_cutoff)
            {
                double c_distance = sqrt(pow(x-center_x,2)+pow(y-center_y,2));

                if(c_distance >= inside_pro*judge_len)
                {
                    double divi_prob_out = unif_real_dis(RNG);
                    if(divi_prob_out < cell_array(i,10))
                    {
      
                        number = 0;
                        for(int i=1;i<=division_range;i++)
                        {
                            for(int j=1;j<=division_range;j++)
                            {
                                number += 1;
                                if((int)sub_visual(i,j,1)==0)
                                {
                                    cell_count.push_back(number);
                                }
                            }
                        }
                        loci_number = (int)cell_count.size();
        
                        if (loci_number!=0)
                        {
                    
                            mut_number1 = mut_pos_distribution(RNG);
                            mut_number_Driver1 = mut_Driver_distribution(RNG);
             

                            mut_number2 = mut_pos_distribution(RNG);
                            mut_number_Driver2 = mut_Driver_distribution(RNG);
                  
                            cell_array(i,4) = mut_number1 + mut_number_Driver1;
                            cell_array(i,7) = cell_label;
                            cell_array(i,8) = cell_copy(1,8) + cell_array(i,4); 
                            cell_array(i,11) = cell_copy(1,11) + 1;
          
                            if(mut_number_Driver1 > 0)
                            {
                               
                                cell_array(i,3) += 1;
        
                            }
                           
                            Visual_range(x,y,4) = cell_label;
                            
                            Mut_temp1(1,1) = cell_label; 
                            Mut_temp1(1,2) = cell_copy(1,7); 
                            Mut_temp1(1,3) = mut_number1 + mut_number_Driver1;
                            cell_id_index[cell_label] = i; 

                 
                            cell_label=cell_label+1; 
                            cell_index=cell_copy(1,5);
                            shuffle(cell_count.begin(), cell_count.end(),RNG);
                            
                            int loci=cell_count[0]; 
                            int loci_x = left_most_x + (loci-half_range)/division_range;
                            int loci_y = left_most_y + (loci-half_range)%division_range;
                            cell_temp(1,1)= loci_x; 
                            cell_temp(1,2)= loci_y; 
                            cell_temp(1,4) = mut_number2 + mut_number_Driver2;
                            cell_temp(1,5) = cell_index; 
                            cell_temp(1,6) = cell_copy(1,6); 
                            cell_temp(1,7) = cell_label; 
                            cell_temp(1,8) = cell_copy(1,8) + cell_temp(1,4); 
                            cell_temp(1,9) = cell_copy(1,9); 
                            cell_temp(1,10) = cell_copy(1,10); 
                            cell_temp(1,11) = cell_copy(1,11) + 1;
                            cell_temp(1,12) = cell_copy(1,12);
                            cell_temp(1,13) = cell_copy(1,13);

                            if(mut_number_Driver2 > 0)
                            {
                         
                                cell_temp(1,3) += 1;

                            }
                            
          
                            Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),1)=1;
                            Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),2)=cell_index;
    
                            Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),4)=cell_label;
                            
                            Mut_temp2(1,1) = cell_label;
                            Mut_temp2(1,2) = cell_copy(1,7); 
                            Mut_temp2(1,3) = mut_number2 + mut_number_Driver2;

            
                            division_time += 1;
                       
                            cell_array_temp(division_time,all) = cell_temp(1,all);
                            cell_id_index[cell_label] = C1 + division_time;
                        
                            cell_label += 1;
                            mut_index = (int)((cell_label-2)/2);

                            Mut_trace(mut_index,Range(4,6)) = Mut_temp1(1,all);
                            Mut_trace(mut_index+1,Range(1,3)) = Mut_temp2(1,all);
                         
                        }
                        else
                        { 
                            cell_array(i,9) = 0;
                        }
                    }
                    else
                    { 
                        cell_array(i,9) = 0;
                    }
                }
                else 
                { 
                    double divi_prob_inn = unif_real_dis(RNG);
    
                    if(divi_prob_inn <= cell_array(i,12) + cell_array(i,13) && divi_prob_inn > cell_array(i,13))
                    {
                        
                        number = 0;
                        for(int i=1;i<=division_range;i++)
                        {
                            for(int j=1;j<=division_range;j++)
                            {
                                number += 1;
                                if((int)sub_visual(i,j,1)==0)
                                {
                                    cell_count.push_back(number);
                                }
                            }
                        }
                        loci_number = (int)cell_count.size();
                      
                        if (loci_number!=0)
                        {
                     
                            mut_number1 = mut_pos_distribution(RNG);
                            mut_number_Driver1 = mut_Driver_distribution(RNG);
                  

                            mut_number2 = mut_pos_distribution(RNG);
                            mut_number_Driver2 = mut_Driver_distribution(RNG);
                   
                            cell_array(i,4) = mut_number1 + mut_number_Driver1;
                            cell_array(i,7) = cell_label;
                            cell_array(i,8) = cell_copy(1,8) + cell_array(i,4); 
                            cell_array(i,11) = cell_copy(1,11) + 1;
                      
                            if(mut_number_Driver1 > 0 && GUARD!=1)
                            {
                                cell_array(i,6) += 1;
                                cell_array(i,12) += cell_array(i,12)*s_coef;
                                if(cell_array(i,12)>1)
                                {
                                    cell_array(i,12) = 1;
                                }
                            }
                            if(cell_array(i,6)>=compe_cutoff)
                            {
                                GUARD=1;
                            }

                  
                            Visual_range(x,y,4) = cell_label;
                            Visual_range(x,y,3) = cell_array(i,6); 
                            
                            Mut_temp1(1,1) = cell_label; 
                            Mut_temp1(1,2) = cell_copy(1,7); 
                            Mut_temp1(1,3) = mut_number1 + mut_number_Driver1;
                            cell_id_index[cell_label] = i; 

         
                            cell_label=cell_label+1; 
                            cell_index=cell_copy(1,5);
                            shuffle(cell_count.begin(), cell_count.end(),RNG);
                            
                            int loci=cell_count[0]; 
                            int loci_x = left_most_x + (loci-half_range)/division_range;
                            int loci_y = left_most_y + (loci-half_range)%division_range;
                            cell_temp(1,1)= loci_x; 
                            cell_temp(1,2)= loci_y; 
                            cell_temp(1,4) = mut_number2 + mut_number_Driver2;
                            cell_temp(1,5) = cell_index; 
                            cell_temp(1,6) = cell_copy(1,6); 
                            cell_temp(1,7) = cell_label; 
                            cell_temp(1,8) = cell_copy(1,8) + cell_temp(1,4); 
                            cell_temp(1,9) = cell_copy(1,9); 
                            cell_temp(1,10) = cell_copy(1,10); 
                            cell_temp(1,11) = cell_copy(1,11) + 1;
                            cell_temp(1,12) = cell_copy(1,12);
                            cell_temp(1,13) = cell_copy(1,13);

                            if(mut_number_Driver2 > 0 && GUARD!=1)
                            {
                                cell_temp(1,6) += 1;
                                cell_temp(1,12) += cell_temp(1,12)*s_coef;
                                if(cell_temp(1,12)>1)
                                {
                                    cell_temp(1,12) = 1;
                                }
                            }
                            
                            if(cell_temp(1,6)>=compe_cutoff)
                            {
                                GUARD=1;
                            }
                            
                            Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),1)=1;
                            Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),2)=cell_index;
                            Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),3)=cell_temp(1,6);
                            Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),4)=cell_label;
                            
                            Mut_temp2(1,1) = cell_label;
                            Mut_temp2(1,2) = cell_copy(1,7); 
                            Mut_temp2(1,3) = mut_number2 + mut_number_Driver2;

                   
                            division_time += 1;

                            cell_array_temp(division_time,all) = cell_temp(1,all);
                            cell_id_index[cell_label] = C1 + division_time;
                       
                            cell_label += 1;
                            mut_index = (int)((cell_label-2)/2);

                            Mut_trace(mut_index,Range(4,6)) = Mut_temp1(1,all);
                            Mut_trace(mut_index+1,Range(1,3)) = Mut_temp2(1,all);
                      
                        }
                   
                    }
                   
                    if (divi_prob_inn > cell_array(i,12) + cell_array(i,13))
                    {
                        cell_array(i,9) = 0;
                    }         

                }

            }
            else 
            {
              
                number = 0;
                for(int i=1;i<=division_range;i++)
                {
                    for(int j=1;j<=division_range;j++)
                    {
                        number += 1;
                        if((int)sub_visual(i,j,1)==0)
                        {
                            cell_count.push_back(number);
                        }
                    }
                }
                loci_number = (int)cell_count.size();
                
                if(loci_number!=0)
                {
                    
                    mut_number1 = mut_pos_distribution(RNG);
                    mut_number_Driver1 = mut_Driver_distribution(RNG);
              

                    mut_number2 = mut_pos_distribution(RNG);
                    mut_number_Driver2 = mut_Driver_distribution(RNG);
                    
                    cell_array(i,4) = mut_number1 + mut_number_Driver1;
                    cell_array(i,7) = cell_label;
                    cell_array(i,8) = cell_copy(1,8) + cell_array(i,4); 
                    cell_array(i,11) = cell_copy(1,11) + 1;
                
                    if(mut_number_Driver1 > 0 && GUARD!=1)
                    {
                        cell_array(i,6) += 1;
                        cell_array(i,12) += cell_array(i,12)*s_coef;
                        if(cell_array(i,12)>1)
                        {
                            cell_array(i,12) = 1;
                        }
                    }
                  
                    Visual_range(x,y,4) = cell_label;
                    Visual_range(x,y,3) = cell_array(i,6); 
                    
                    Mut_temp1(1,1) = cell_label; 
                    Mut_temp1(1,2) = cell_copy(1,7); 
                    Mut_temp1(1,3) = mut_number1 + mut_number_Driver1;
                    cell_id_index[cell_label] = i; 

               
                    cell_label=cell_label+1; 
                    cell_index=cell_copy(1,5);
                    shuffle(cell_count.begin(), cell_count.end(),RNG);
                    
                    int loci=cell_count[0]; 
                    int loci_x = left_most_x + (loci-half_range)/division_range;
                    int loci_y = left_most_y + (loci-half_range)%division_range;
                    cell_temp(1,1)= loci_x; 
                    cell_temp(1,2)= loci_y; 
                    cell_temp(1,4) = mut_number2 + mut_number_Driver2;
                    cell_temp(1,5) = cell_index; 
                    cell_temp(1,6) = cell_copy(1,6); 
                    cell_temp(1,7) = cell_label; 
                    cell_temp(1,8) = cell_copy(1,8) + cell_temp(1,4); 
                    cell_temp(1,9) = cell_copy(1,9); 
                    cell_temp(1,10) = cell_copy(1,10); 
                    cell_temp(1,11) = cell_copy(1,11) + 1;
                    cell_temp(1,12) = cell_copy(1,12);
                    cell_temp(1,13) = cell_copy(1,13);

                    if(mut_number_Driver2 > 0 && GUARD!=1)
                    {
                        cell_temp(1,6) += 1;
                        cell_temp(1,12) += cell_temp(1,12)*s_coef;
                        if(cell_temp(1,12)>1)
                        {
                            cell_temp(1,12) = 1;
                        }
                    }
                    
                    Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),1)=1;
                    Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),2)=cell_index;
                    Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),3)=cell_temp(1,6); 
                    Visual_range((int)cell_temp(1,1),(int)cell_temp(1,2),4)=cell_label;
                    
                    Mut_temp2(1,1) = cell_label;
                    Mut_temp2(1,2) = cell_copy(1,7); 
                    Mut_temp2(1,3) = mut_number2 + mut_number_Driver2;

                    division_time += 1;

                    cell_array_temp(division_time,all) = cell_temp(1,all);
                    cell_id_index[cell_label] = C1 + division_time;
                  
                    cell_label += 1;
                    mut_index = (int)((cell_label-2)/2);

                    Mut_trace(mut_index,Range(4,6)) = Mut_temp1(1,all);
                    Mut_trace(mut_index+1,Range(1,3)) = Mut_temp2(1,all);
                
                }
                else 
                {
                    unsigned int cell_id;
                    int dri_number;
                    vector<int> not_empty;
                    vector<int> adv_number;
                    for(int i=1;i<=division_range;i++)
                    {
                        for(int j=1;j<=division_range;j++)
                        {
                            cell_id = sub_visual(i,j,4); 
                            dri_number = sub_visual(i,j,3);
                            if((int)sub_visual(i,j,1) != 0 && sub_visual(i,j,3) < compe_cutoff)
                            {
                                not_empty.push_back(cell_id);
                                adv_number.push_back(dri_number);
                            }
                        }
                    }
                    int not_empty_number = (int)not_empty.size();
              
                    if(not_empty_number > 0)
                    {
              
                            mut_number1 = mut_pos_distribution(RNG);
                            mut_number_Driver1 = mut_Driver_distribution(RNG);
                           
                            mut_number2 = mut_pos_distribution(RNG);
                            mut_number_Driver2 = mut_Driver_distribution(RNG);
                                    
                            cell_array(i,4) = mut_number1 + mut_number_Driver1;
                            cell_array(i,7) = cell_label;
                            cell_array(i,8) = cell_copy(1,8) + cell_array(i,4); 
                            cell_array(i,11) = cell_copy(1,11) + 1;
                            
                            if(mut_number_Driver1 > 0)
                            {
                                cell_array(i,6) += 1;
                                cell_array(i,12) += cell_array(i,12)*s_coef;
                                if(cell_array(i,12)>1)
                                {
                                    cell_array(i,12) = 1;
                                }
                            }
                            
                            Visual_range(x,y,4) = cell_label;
                            Visual_range(x,y,3) = cell_array(i,6); 
                            
                            Mut_temp1(1,1) = cell_label; 
                            Mut_temp1(1,2) = cell_copy(1,7); 
                            Mut_temp1(1,3) = mut_number1 + mut_number_Driver1;
                            cell_id_index[cell_label] = i; 

                            cell_label=cell_label+1;
                            cell_index=cell_copy(1,5);
                            shuffle(not_empty.begin(),not_empty.end(),RNG);
                            int cell_array_id = not_empty[0];
                       
                            int cell_array_index = cell_id_index[cell_array_id];
                           
                            int curr_x,curr_y;
                            
                            if(cell_array_index <= C1)
                            { 
                                
                                cell_array(cell_array_index,4) = mut_number2 + mut_number_Driver2;
                                cell_array(cell_array_index,5) = cell_copy(1,5);
                                cell_array(cell_array_index,6) = cell_copy(1,6);
                                cell_array(cell_array_index,7) = cell_label;
                                cell_array(cell_array_index,8) = cell_copy(1,8) + cell_array(cell_array_index,4);
                                cell_array(cell_array_index,9) = cell_copy(1,9);
                                cell_array(cell_array_index,10) = cell_copy(1,10);
                                cell_array(cell_array_index,11) = cell_copy(1,11) + 1;
                                cell_array(cell_array_index,12) = cell_copy(1,12);
                                cell_array(cell_array_index,13) = cell_copy(1,13);

                               
                                if(mut_number_Driver2 > 0 && GUARD!=1)
                                {
                                    cell_array(cell_array_index,6) += 1;
                                    cell_array(cell_array_index,12) += cell_array(cell_array_index,12)*s_coef;
                                    if(cell_array(cell_array_index,12)>1)
                                    {
                                        cell_array(cell_array_index,12) = 1;
                                    }
                                }
                                curr_x = cell_array(cell_array_index,1);
                                curr_y = cell_array(cell_array_index,2);
                                Visual_range(curr_x,curr_y,3) = cell_array(cell_array_index,6);
                
                                cell_id_index[cell_label] = cell_array_index;
                            }else
                            { 
                                int cell_temp_index = cell_array_index - C1;
                                cell_array_temp(cell_temp_index,4) = mut_number2 + mut_number_Driver2;
                                cell_array_temp(cell_temp_index,5) = cell_copy(1,5);
                                cell_array_temp(cell_temp_index,6) = cell_copy(1,6);
                                cell_array_temp(cell_temp_index,7) = cell_label;
                                cell_array_temp(cell_temp_index,8) = cell_copy(1,8) + cell_array_temp(cell_temp_index,4);
                                cell_array_temp(cell_temp_index,9) = cell_copy(1,9);
                                cell_array_temp(cell_temp_index,10) = cell_copy(1,10);
                                cell_array_temp(cell_temp_index,11) = cell_copy(1,11) + 1;
                                cell_array_temp(cell_temp_index,12) = cell_copy(1,12);
                                cell_array_temp(cell_temp_index,13) = cell_copy(1,13);

                                if(mut_number_Driver2 > 0 && GUARD!=1)
                                {
                                    cell_array_temp(cell_temp_index,6) += 1;
                                    cell_array_temp(cell_temp_index,12) += cell_array_temp(cell_temp_index,12)*s_coef;
                                    if(cell_array_temp(cell_temp_index,12) > 1)
                                    {
                                        cell_array_temp(cell_temp_index,12) = 1;
                                    }
                                }
                                curr_x = cell_array_temp(cell_temp_index,1);
                                curr_y = cell_array_temp(cell_temp_index,2);
                                Visual_range(curr_x,curr_y,3) = cell_array_temp(cell_temp_index,6); 
                                cell_id_index[cell_label] = cell_array_index;
                            }

                            
                            Visual_range(curr_x,curr_y,1) = 1;
                            Visual_range(curr_x,curr_y,2) = cell_index;
                            Visual_range(curr_x,curr_y,4) = cell_label;
                            
                            Mut_temp2(1,1) = cell_label;
                            Mut_temp2(1,2) = cell_copy(1,7); 
                            Mut_temp2(1,3) = mut_number2 + mut_number_Driver2;

        
                            cell_label += 1;
                            mut_index = (int)((cell_label-2)/2);

                            Mut_trace(mut_index,Range(4,6)) = Mut_temp1(1,all);
                            Mut_trace(mut_index+1,Range(1,3)) = Mut_temp2(1,all);
        

                    }
                    vector<int>().swap(not_empty); 
                }
            }
        }
        vector<int>().swap(cell_count);
    }
 
    if(division_time>=1)
    {
        cell_array_temp2.resize(C1+division_time,cell_column);
        for(int i=1;i<=division_time+C1;i++)
        {
            if(i<=C1)
            {
                cell_array_temp2(i,all) = cell_array(i,all);
            }else
            {
                cell_array_temp2(i,all) = cell_array_temp(i-C1,all);
            }
        }
        cell_array.resize(C1+division_time,cell_column);
        cell_array(all,all) = cell_array_temp2(all,all);
        cell_array_temp2.resize(1,cell_column);
        cell_array_temp.resize(1,cell_column);
    }
    
}

#endif /* division_invasion_slow_hpp */
