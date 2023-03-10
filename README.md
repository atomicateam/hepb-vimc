# hepb-vimc

The hepatitis B model is run using Python (≥v3.7). Make sure you have the packages numpy, pandas, matplotlib and Atomica ( https://github.com/atomicateam/atomica) installed before you run the model.

Notebook to run:
-	`Run_Model_v3.ipynb` demonstrates how to run the model in order to collect the Parameter Set, Central Results, and Stochastic Results. All result sets may be saved as a csv afterwards
-	`Plot_All_Results.ipynb` runs the Hepatitis B model separately and allows you to plot various outputs in the model. In addition, there are functions that can be used to view plots of the central and stochastic results, as well as a boxplot of each input parameter in the parameter 


The notebook creates Pandas DataFrame results relating to the VIMC outputs
- Input results
- Central Output Results
- Stochastic Output Results

The central estimate of our model is generated by using all mean input parameter values. Each column in the central estimate sheet is constructed following the VIMC guidelines. Specifically, each burden outcome is extracted as follows:
- Cohort Size
  - Cohort size is the total number of people in the age group in the region during a given year. 
  - For each model age group, we used the VIMC Interpolated (1-year time and age) population dataset to distribute the model’s age groups into 1-year age bins 
- Total Cases
  - Total cases are taken from the variable tot_inc, which takes the total incidence in the age-group in a given year
  - Similar to cohort size, since this variable is taken for each age group and distributed using the VIMC Interpolated (1-year time and age) population datasets
- DALYs
  - DALYs are taken from the variable, dalys, which calculate the total DALYs accrued per age group in a given year
  - Similar to cohort size, since this variable is taken for each age group and distributed using the VIMC Interpolated (1-year time and age) population dataset 
- Deaths
  - Deaths are taken as a sum of deaths associated with acute hepatitis B (cl_acu), decompensated cirrhosis (cl_dc), and Hepatocellular carcinoma (cl_hcc) in an age group at a given year
  - Similar to cohort size, since this variable is taken for each age group and distributed using the VIMC Interpolated (1-year time and age) population dataset

Stochastic Estimates

Each stochastic estimate uses one of the 30 parameter sets as part of the input parameter set to generate a unique model outcome. Besides the run­_id­ column, all the columns of the stochastic result spreadsheet are defined in the same way as the central estimate sheet.

All stochastic estimate results of the same scenario are found in the same spreadsheet, labelled either stochastic_burden_est_hepb_burnet_novax.csv (No vaccination scenario), and stochastic_burden_est_hepb_burnet_default.csv (Default scenario), with each run within each spreadsheet characterised by their run_id­.
