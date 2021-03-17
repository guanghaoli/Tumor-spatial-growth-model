# Tumor-spatial-growth-model

The codes are used to model tumor growth in spatial structure (2D). With different evolutionary driven forces, different main programs are applied. File "Neutral.cpp" is used to simulate tumor growth under neutral evolution. "Selection.cpp" is used to simulate tumor growth under natural selection.

Several C++ libraries are needed to complile the codes, including blitz++, PNG Writer and gsl.

Usage: ./Tumor_growth driver_number mutation_rate inner_cell_growth_rate quiescent_rate outer_cell_growth_rate distance_proportion_of_inner_cell
For example: ./Tumor_growth 4 15 0.3 0.5 0.53 0.9

Need Help? Contact me via e-mail: liguanghaoli@163.com
