import numpy as np
import pandas as pd
import atomica as at

def load_cen_project(fw = 'hbv_v14_gamma_2.xlsx', db = "AFR_db_v1_2_1.xlsx", cl = 'AFR_calib_v1_2.xlsx'):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Loads an Atomica Project with the results using the baseline input parameters

    Inputs:
    - fw:   Excel sheet of the Atomica Framework of the model (str)
    - db:   Excel sheet of the Atomica Databook of the model (str)
    - cl:   Excel sheet of the Model Calibration of the model - specific to a region estimate (str)

    Outputs:
    - P:    Atomica Project attributed to the inputted framework and databook (Project)
    - res:  Project results from the model (Result)
    '''
    ## Load projects and input parameters
    P =at.Project(framework=fw, databook=db, sim_start=1990, sim_end=2101, sim_dt=0.25, do_run=False)

    cal = P.make_parset()
    cal.load_calibration(cl)

    res = P.run_sim(parset=cal, result_name = 'Status Quo')

    return P, res

def load_sto_project(n_samples, fw = 'hbv_v14_gamma_2.xlsx', db = "AFR_db_v1_2_1.xlsx", cl = 'AFR_calib_v1_2.xlsx', seed = 310123):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Loads a set of Atomica Projects with the results using the baseline input parameters with variances in select input parameters

    Inputs:
    - n_samples:    Number of stochastic input samples (int)
    - fw:           Excel sheet of the Atomica Framework of the model (str)
    - db:           Excel sheet of the Atomica Databook of the model (str)
    - cl:           Excel sheet of the Model Calibration of the model - specific to a region estimate (str)
    - seed:         Random seed (int)

    Outputs:
    - P:        Atomica Project attributed to the inputted framework and databook (Project)
    - parsets:  Dictionary of input parameter sets (dict)
    - results:  Dictionary of project results from the model (dict)
    '''
    ## Load projects and input parameters
    P =at.Project(framework=fw, databook=db, sim_start=1990, sim_end=2101, sim_dt=0.25, do_run=False)

    cal = P.make_parset()
    cal.load_calibration(cl)

    np.random.seed(seed) # Set a random seed for consistent results
    afr_ua=P.parsets[0]

    parsets = {}
    results = {}

    ## Sample a set of parameters
    for i in range(1, n_samples+1):
        parsets[i] = afr_ua.sample()

    ## Copy the male treatment/vaccine parameters onto the female parameters
    ages = ['0-4', '5-14', '15-49', '50-69', '70+'] # Age groups
    nathis = ['ci_p', 'm_acu', 'm_dc', 'm_hcc']     # Natural History-related paramters

    trtinv = ['te_dc_cc', 'te_icl_ict', 'te_icl_cc', 'te_cc_dc', 'te_cc_hcc', 'te_ie_cc', 'te_m_dc', 'te_m_hcc', \
              'te_ict_hcc', 'te_ie_hcc', 'te_icl_hcc', 'te_dc_hcc'] # Treatment-related efficacy parameters
    vacinv = ['eag_ve', 'sag_ve', 'hb3_ve', 'mav_ve'] # Vaccination-related efficacy parameters
    tot_pars = trtinv + vacinv
    tot_pars.remove('te_icl_ict')

    # Ensure that the randomly sampled input parameters are reasonable (nathis vals >=0, others between 0 and 1)
    for i in range(1, n_samples+1):
        for age in P.data.pops:
            for par in nathis + ['te_icl_ict']:
                if parsets[i].get_par(par).ts[age].assumption < 0:
                    parsets[i].get_par(par).ts[age].assumption = 0

            for par in tot_pars:
                if parsets[i].get_par(par).ts[age].assumption < 0:
                    parsets[i].get_par(par).ts[age].assumption = 0
                elif parsets[i].get_par(par).ts[age].assumption > 1:
                    parsets[i].get_par(par).ts[age].assumption = 1

    # Duplicate treatment and vaccine efficacy parameters for male onto females
    for i in range(1, n_samples+1):
        for age in ages:
            for par in tot_pars:
                parsets[i].get_par(par).ts[f'{age}F'].assumption = parsets[i].get_par(par).ts[f'{age}M'].assumption

    ## Generate Results
    for i in range(1, n_samples+1):
        results[i] = P.run_sim(parsets[i], result_name = f'sample result {i}')

    return P, parsets, results

def input_results(n_samples, fw = 'hbv_v14_gamma_2.xlsx', db = "AFR_db_v1_2_1.xlsx", cl = 'AFR_calib_v1_2.xlsx', seed = 310123):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Generates a DataFrame listing all input results.

    Inputs:
    - n_samples:    Number of stochastic input samples (int)
    - fw:           Excel sheet of the Atomica Framework of the model (str)
    - db:           Excel sheet of the Atomica Databook of the model (str)
    - cl:           Excel sheet of the Model Calibration of the model - specific to a region estimate (str)
    - seed:         Random seed (int)

    Outputs:
    - in_df:        DataFrame listing all stochastic parameter sets (DataFrame)
    '''

    ## Load projects and input parameters
    P, parsets, results = load_sto_project(n_samples, fw= fw, db=db, cl=cl, seed=seed)

    sexes = ['M', 'F']
    ages = ['0-4', '5-14', '15-49', '50-69', '70+']
    nathis = ['ci_p', 'm_acu', 'm_dc', 'm_hcc']
    trtinv = ['te_dc_cc', 'te_icl_ict', 'te_icl_cc', 'te_cc_dc', 'te_cc_hcc', 'te_ie_cc', 'te_m_dc', 'te_m_hcc', 'te_ict_hcc', 'te_ie_hcc', 'te_icl_hcc', 'te_dc_hcc']
    vacinv = ['eag_ve', 'sag_ve', 'hb3_ve', 'mav_ve']
    all_pars = trtinv + vacinv + nathis

    # Create a list of input parameters (which vary my age and sex)
    in_df_pars = []
    for par in all_pars:
        for age in ages:
            for sex in sexes:
                in_df_pars.append(f'{par}_{age}{sex}')

    ## Create a dataframe for the input parameters
    in_df = pd.DataFrame(columns = ['run_id']+in_df_pars)

    # Input all parameters
    for sim in range(1,n_samples+1):
        in_df.loc[sim,'run_id'] = sim
        for par in all_pars:
            for age in ages:
                for sex in sexes:
                    in_df.loc[sim,f'{par}_{age}{sex}'] = parsets[sim].get_par(par).ts[f'{age}{sex}'].assumption #TODO: get different input parameters per population.loc[sim,col] = parsets[sim].get_par(col).ts[0].assumption #TODO: get different input parameters per population

    # Drop 0-value columns
    drop_cols = []
    for col in in_df.columns:
        if in_df.loc[:,col].sum() == 0:
            drop_cols.append(col)

    in_df.drop(columns = drop_cols, inplace = True)

    return in_df

def central_results(fw = 'hbv_v14_gamma_2.xlsx', db = "AFR_db_v1_2_1.xlsx", cl = 'AFR_calib_v1_2.xlsx'):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Generates a DataFrame listing all central output results.

    Inputs:
    - fw:       Excel sheet of the Atomica Framework of the model (str)
    - db:       Excel sheet of the Atomica Databook of the model (str)
    - cl:       Excel sheet of the Model Calibration of the model - specific to a region estimate (str)

    Outputs:
    - cen_df:   DataFrame listing the central estimates of 1-year age groups between 2000-2100 (DataFrame)
    '''
    P, res = load_cen_project(fw= fw, db=db, cl=cl)
    ## Locations for the important dataframes
    loc = 'Data/Templates/' # The location for the output templates may be different
    loc2 = 'Data/Demographics/' # The location for the input data may be different (depending on where you store it)

    ## Load and process the input data
    df1 = pd.read_csv(loc2+'202212rfp-1_dds-202208_int_pop_both.csv')

    df1 = df1[(df1.year >= 1990) & (df1.year <=2100)]
    age_groups = ['0-4', '5-14', '15-49', '50-69', '70+']

    ## Assort data by age group
    for idx,row in df1.iterrows():
        if row.age_from < 5:
            df1.loc[idx,'age_group'] = '0-4'
        elif row.age_from < 15:
            df1.loc[idx,'age_group'] = '5-14'
        elif row.age_from < 50:
            df1.loc[idx,'age_group'] = '15-49'
        elif row.age_from < 70:
            df1.loc[idx,'age_group'] = '50-69'
        else:
            df1.loc[idx,'age_group'] = '70+'


    ## Start generating output results
    output_dict = {'cohort_size': 'alive', 'cases': 'tot_inc', 'dalys':'dalys'}
    cen_df = pd.read_csv(loc+'/central-burden-template.202212rfp-1.RFP_standard template_HepB.csv')

    for opt in output_dict:
        tot_val = res.get_variable(output_dict[opt])[0].vals + res.get_variable(output_dict[opt])[1].vals
        for age in range(0,5):
            for year in range(2000,2101):
                ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('0-4',year),'value'] # Ratio of age relative to age group
                res_val = tot_val[(year-1990)*4] # The '4' should be timestep (if future runs are different)

                idx = cen_df[(cen_df.year == year) & (cen_df.age == age)].index[0] #Index of age group year sim entry
                cen_df.loc[idx, opt] = ratio * res_val

        tot_val = res.get_variable(output_dict[opt])[2].vals + res.get_variable(output_dict[opt])[3].vals
        for age in range(5,15):
            for year in range(2000,2101):
                ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('5-14',year),'value']
                res_val = tot_val[(year-1990)*4]

                idx = cen_df[(cen_df.year == year) & (cen_df.age == age)].index[0]
                cen_df.loc[idx, opt] = ratio * res_val

        tot_val = res.get_variable(output_dict[opt])[4].vals + res.get_variable(output_dict[opt])[5].vals
        for age in range(15,50):
            for year in range(2000,2101):
                ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('15-49',year),'value']
                res_val = tot_val[(year-1990)*4]

                idx = cen_df[(cen_df.year == year) & (cen_df.age == age)].index[0]
                cen_df.loc[idx, opt] = ratio * res_val

        tot_val = res.get_variable(output_dict[opt])[6].vals + res.get_variable(output_dict[opt])[7].vals
        for age in range(50,70):
            for year in range(2000,2101):
                ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('50-69',year),'value']
                res_val = tot_val[(year-1990)*4]

                idx = cen_df[(cen_df.year == year) & (cen_df.age == age)].index[0]
                cen_df.loc[idx, opt] = ratio * res_val

        tot_val = res.get_variable(output_dict[opt])[8].vals + res.get_variable(output_dict[opt])[9].vals
        for age in range(70,101):
            for year in range(2000,2101):
                ratio = float(df1[(df1.age_from == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('70+',year),'value']
                res_val = tot_val[(year-1990)*4]

                idx = cen_df[(cen_df.year == year) & (cen_df.age == age)].index[0]
                cen_df.loc[idx, opt] = ratio * res_val

    ## Deaths are calculated separately
    tot_val = 0
    for var in ['cl_acu', 'cl_cir', 'cl_hcc']:
        tot_val += res.get_variable(var)[0].vals + res.get_variable(var)[1].vals # Total deaths in age group

    for age in range(0,5):
        for year in range(2000,2101):
            ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('0-4',year),'value']
            res_val = tot_val[(year-1990)*4]

            idx = cen_df[(cen_df.year == year) & (cen_df.age == age)].index[0]
            cen_df.loc[idx, 'deaths'] = ratio * res_val

    tot_val = 0
    for var in ['cl_acu', 'cl_cir', 'cl_hcc']:
        tot_val += res.get_variable(var)[2].vals + res.get_variable(var)[3].vals
    for age in range(5,15):
        for year in range(2000,2101):
            ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('5-14',year),'value']
            res_val = tot_val[(year-1990)*4]

            idx = cen_df[(cen_df.year == year) & (cen_df.age == age)].index[0]
            cen_df.loc[idx, 'deaths'] = ratio * res_val

    tot_val = 0
    for var in ['cl_acu', 'cl_cir', 'cl_hcc']:
        tot_val += res.get_variable(var)[4].vals + res.get_variable(var)[5].vals
    for age in range(15,50):
        for year in range(2000,2101):
            ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('15-49',year),'value']
            res_val = tot_val[(year-1990)*4]

            idx = cen_df[(cen_df.year == year) & (cen_df.age == age)].index[0]
            cen_df.loc[idx, 'deaths'] = ratio * res_val

    tot_val = 0
    for var in ['cl_acu', 'cl_cir', 'cl_hcc']:
        tot_val += res.get_variable(var)[6].vals + res.get_variable(var)[7].vals
    for age in range(50,70):
        for year in range(2000,2101):
            ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('50-69',year),'value']
            res_val = tot_val[(year-1990)*4]

            idx = cen_df[(cen_df.year == year) & (cen_df.age == age)].index[0]
            cen_df.loc[idx, 'deaths'] = ratio * res_val

    tot_val = 0
    for var in ['cl_acu', 'cl_cir', 'cl_hcc']:
        tot_val += res.get_variable(var)[8].vals + res.get_variable(var)[9].vals
    for age in range(70,101):
        for year in range(2000,2101):
            ratio = float(df1[(df1.age_from == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('70+',year),'value']
            res_val = tot_val[(year-1990)*4]

            idx = cen_df[(cen_df.year == year) & (cen_df.age == age)].index[0]
            cen_df.loc[idx, 'deaths'] = ratio * res_val

    return cen_df

def stochastic_results(n_samples, fw = 'hbv_v14_gamma_2.xlsx', db = "AFR_db_v1_2_1.xlsx", cl = 'AFR_calib_v1_2.xlsx', seed = 310123):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Generates a DataFrame listing all stochastic output results. Each stochastic output is indexed by a 'run_id' which can also refer to the stochastic input parameter set.

    Inputs:
    - n_samples:    Number of stochastic input samples (int)
    - fw:           Excel sheet of the Atomica Framework of the model (str)
    - db:           Excel sheet of the Atomica Databook of the model (str)
    - cl:           Excel sheet of the Model Calibration of the model - specific to a region estimate (str)
    - seed:         Random seed (int)

    Outputs:
    - final_df:     DataFrame listing all stochastic estimates of 1-year age groups between 2000-2100 (DataFrame)
    '''
    P, parsets, results = load_sto_project(n_samples, fw= fw, db=db, cl=cl, seed=seed)
    ## Locations for the important dataframes
    loc = 'Data/Templates/' # The location for the output templates may be different
    loc2 = 'Data/Demographics/' # The location for the input data may be different (depending on where you store it)

    df = pd.read_csv(loc+'/stochastic-burden-template.202212rfp-1.RFP_standard template_HepB.csv') #template

     ## Load and process the input data
    df1 = pd.read_csv(loc2+'202212rfp-1_dds-202208_int_pop_both.csv')

    df1 = df1[(df1.year >= 1990) & (df1.year <=2100)]
    age_groups = ['0-4', '5-14', '15-49', '50-69', '70+']

    for idx,row in df1.iterrows():
        if row.age_from < 5:
            df1.loc[idx,'age_group'] = '0-4'
        elif row.age_from < 15:
            df1.loc[idx,'age_group'] = '5-14'
        elif row.age_from < 50:
            df1.loc[idx,'age_group'] = '15-49'
        elif row.age_from < 70:
            df1.loc[idx,'age_group'] = '50-69'
        else:
            df1.loc[idx,'age_group'] = '70+'

    dfs = {}
    output_dict = {'cohort_size': 'alive', 'cases': 'tot_inc', 'dalys':'dalys'}
    for sim in range(1,n_samples+1):
        stores = results[sim]
        print(f'Sim no {sim}')

        dfs[sim] = pd.read_csv(loc+'/stochastic-burden-template.202212rfp-1.RFP_standard template_HepB.csv')
        dfs[sim].run_id = sim
        ## Other outputs
        for opt in output_dict:
            print(opt)
            tot_val = stores.get_variable(output_dict[opt])[0].vals + stores.get_variable(output_dict[opt])[1].vals
            for age in range(0,5):
                for year in range(2000,2101):
                    ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('0-4',year),'value'] # Ratio of age relative to age group
                    res_val = tot_val[(year-1990)*4] # '4' refers to the timestep, 1990 refers to the beginning of the simulation

                    idx = df[(dfs[sim].year == year) & (dfs[sim].age == age)].index[0]
                    dfs[sim].loc[idx, opt] = ratio * res_val

            tot_val = stores.get_variable(output_dict[opt])[2].vals + stores.get_variable(output_dict[opt])[3].vals
            for age in range(5,15):
                for year in range(2000,2101):
                    ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('5-14',year),'value']
                    res_val = tot_val[(year-1990)*4]

                    idx = df[(dfs[sim].year == year) & (dfs[sim].age == age)].index[0]
                    dfs[sim].loc[idx, opt] = ratio * res_val

            tot_val = stores.get_variable(output_dict[opt])[4].vals + stores.get_variable(output_dict[opt])[5].vals
            for age in range(15,50):
                for year in range(2000,2101):
                    ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('15-49',year),'value']
                    res_val = tot_val[(year-1990)*4]

                    idx = df[(dfs[sim].year == year) & (dfs[sim].age == age)].index[0]
                    dfs[sim].loc[idx, opt] = ratio * res_val

            tot_val = stores.get_variable(output_dict[opt])[6].vals + stores.get_variable(output_dict[opt])[7].vals
            for age in range(50,70):
                for year in range(2000,2101):
                    ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('50-69',year),'value']
                    res_val = tot_val[(year-1990)*4]

                    idx = df[(dfs[sim].year == year) & (dfs[sim].age == age)].index[0]
                    dfs[sim].loc[idx, opt] = ratio * res_val

            tot_val = stores.get_variable(output_dict[opt])[8].vals + stores.get_variable(output_dict[opt])[9].vals
            for age in range(70,101):
                for year in range(2000,2101):
                    ratio = float(df1[(df1.age_from == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('70+',year),'value']
                    res_val = tot_val[(year-1990)*4]

                    idx = df[(dfs[sim].year == year) & (dfs[sim].age == age)].index[0]
                    dfs[sim].loc[idx, opt] = ratio * res_val

        ## Deaths
        tot_val = 0
        for var in ['cl_acu', 'cl_cir', 'cl_hcc']:
            tot_val += stores.get_variable(var)[0].vals + stores.get_variable(var)[1].vals

        for age in range(0,5):
            for year in range(2000,2101):
                ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('0-4',year),'value']
                res_val = tot_val[(year-1990)*4]

                idx = df[(dfs[sim].year == year) & (dfs[sim].age == age)].index[0]
                dfs[sim].loc[idx, 'deaths'] = ratio * res_val

        tot_val = 0
        for var in ['cl_acu', 'cl_cir', 'cl_hcc']:
            tot_val += stores.get_variable(var)[2].vals + stores.get_variable(var)[3].vals
        for age in range(5,15):
            for year in range(2000,2101):
                ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('5-14',year),'value']
                res_val = tot_val[(year-1990)*4]

                idx = df[(dfs[sim].year == year) & (dfs[sim].age == age)].index[0]
                dfs[sim].loc[idx, 'deaths'] = ratio * res_val

        tot_val = 0
        for var in ['cl_acu', 'cl_cir', 'cl_hcc']:
            tot_val += stores.get_variable(var)[4].vals + stores.get_variable(var)[5].vals
        for age in range(15,50):
            for year in range(2000,2101):
                ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('15-49',year),'value']
                res_val = tot_val[(year-1990)*4]

                idx = df[(dfs[sim].year == year) & (dfs[sim].age == age)].index[0]
                dfs[sim].loc[idx, 'deaths'] = ratio * res_val

        tot_val = 0
        for var in ['cl_acu', 'cl_cir', 'cl_hcc']:
            tot_val += stores.get_variable(var)[6].vals + stores.get_variable(var)[7].vals
        for age in range(50,70):
            for year in range(2000,2101):
                ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('50-69',year),'value']
                res_val = tot_val[(year-1990)*4]

                idx = df[(dfs[sim].year == year) & (dfs[sim].age == age)].index[0]
                dfs[sim].loc[idx, 'deaths'] = ratio * res_val

        tot_val = 0
        for var in ['cl_acu', 'cl_cir', 'cl_hcc']:
            tot_val += stores.get_variable(var)[8].vals + stores.get_variable(var)[9].vals
        for age in range(70,101):
            for year in range(2000,2101):
                ratio = float(df1[(df1.age_from == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('70+',year),'value']
                res_val = tot_val[(year-1990)*4]

                idx = df[(dfs[sim].year == year) & (dfs[sim].age == age)].index[0]
                dfs[sim].loc[idx, 'deaths'] = ratio * res_val

    final_df = dfs[1]
    for i in range(2,n_samples+1):
        final_df = pd.concat([final_df,dfs[i]])

    return final_df

def plot_output(df, age, run_id, output):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Plots a time series of a specific output within an age group of a specific run_id among the stochastic outputs

    Inputs:
    - df:       Stochastic Result Dataframe (DataFrame)
    - age:      Age of result to plot (int)
    - run_id:   Run_id of result to plot (int)
    - output:   Output column of result to plot - either 'cohort_size', 'cases', 'dalys', 'deaths' (str)
    '''
    plt.plot(np.arange(2000,2101), df[(df.age == age) & (df.run_id == run_id)][output])

def plot_cen_output(df, age, output):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Plots a time series of a specific output within an age group of among the central outputs

    Inputs:
    - df:       Central Result Dataframe (DataFrame)
    - age:      Age of result to plot (int)
    - output:   Output column of result to plot - either 'cohort_size', 'cases', 'dalys', 'deaths' (str)
    '''
    plt.plot(np.arange(2000,2101), df[(df.age == age)][output],color = 'black',linewidth=2)

def plot_all(cen_df, sto_df, age):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Plots time series of all outputs from the central run and all stochastic runs. The central run is specified by a thicker black line, while all stochastic outputs are characterised by thinner time series. Ideally up to 100 stochastic outputs will generate the best idea of the distribution from the uncertainty analysis.

    Inputs:
    - cen_df:   Central Result Dataframe (DataFrame)
    - sto_df:   Stochastic Result Dataframe (DataFrame)
    - age:      Age of result to plot (int)
    '''
    plt.figure()
    for i in range(1,31):
        plot_output(sto_df, age, i, 'cohort_size')
        plot_cen_output(cen_df, age, 'cohort_size')
        plt.title(f'Cohort Size, age {age}')

    plt.figure()
    for i in range(1,31):
        plot_output(sto_df, age, i, 'cases')
        plot_cen_output(cen_df, age, 'cases')
        plt.title(f'Cases, age {age}')

    plt.figure()
    for i in range(1,31):
        plot_output(sto_df, age, i, 'dalys')
        plot_cen_output(cen_df, age, 'dalys')
        plt.title(f'DALYs, age {age}')


    plt.figure()
    for i in range(1,31):
        plot_output(sto_df, age, i, 'deaths')
        plot_cen_output(cen_df, age, 'deaths')
        plt.title(f'Total Deaths, age {age}')
