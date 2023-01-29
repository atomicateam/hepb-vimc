import numpy as np
import pandas as pd
import atomica as at

def load_cen_project(fw = 'hbv_v14_gamma_2.xlsx', db = "AFR_db_v1_2_1.xlsx", cl = 'AFR_calib_v1_2.xlsx'):
    ## Load projects and input parameters
    P =at.Project(framework=fw, databook=db, sim_start=1990, sim_end=2101, sim_dt=0.25, do_run=False) #test PSA on some parameters (vax eff, mtct_infx, mortality rates)

    cal = P.make_parset()
    cal.load_calibration(cl)

    res = P.run_sim(parset=cal, result_name = 'Status Quo')

    return P, res

def load_sto_project(n_samples, fw = 'hbv_v14_gamma_2.xlsx', db = "AFR_db_v1_2_1.xlsx", cl = 'AFR_calib_v1_2.xlsx', seed = 310123):
    ## Load projects and input parameters
    P =at.Project(framework=fw, databook=db, sim_start=1990, sim_end=2101, sim_dt=0.25, do_run=False) #test PSA on some parameters (vax eff, mtct_infx, mortality rates)

    cal = P.make_parset()
    cal.load_calibration(cl)

    np.random.seed(seed)
    afr_ua=P.parsets[0]

    parsets = {}
    results = {}

    ## Sample a set of parameters
    for i in range(1, n_samples+1):
        parsets[i] = afr_ua.sample()

    ## Copy the male treatment/vaccine parameters onto the female parameters
    ages = ['0-4', '5-14', '15-49', '50-69', '70+']
    trtinv = ['te_dc_cc', 'te_icl_ict', 'te_icl_cc', 'te_cc_dc', 'te_cc_hcc', 'te_ie_cc', 'te_m_dc', 'te_m_hcc', 'te_ict_hcc', 'te_ie_hcc', 'te_icl_hcc', 'te_dc_hcc']
    vacinv = ['eag_ve', 'sag_ve', 'hb3_ve', 'mav_ve']
    tot_pars = trtinv + vacinv

    for i in range(1, n_samples+1):
        for age in ages:
            for par in tot_pars:
                parsets[i].get_par(par).ts[f'{age}F'].assumption = parsets[i].get_par(par).ts[f'{age}M'].assumption

    ## Generate Results
    for i in range(1, n_samples+1):
        results[i] = P.run_sim(parsets[i], result_name = f'sample result {i}')

    return P, parsets, results

def input_results(n_samples, fw = 'hbv_v14_gamma_2.xlsx', db = "AFR_db_v1_2_1.xlsx", cl = 'AFR_calib_v1_2.xlsx', seed = 310123):

    ## Load projects and input parameters
    P, parsets, results = load_sto_project(n_samples, fw= fw, db=db, cl=cl, seed=seed)

    sexes = ['M', 'F']
    ages = ['0-4', '5-14', '15-49', '50-69', '70+']
    nathis = ['ci_p', 'm_acu', 'm_dc', 'm_hcc']
    trtinv = ['te_dc_cc', 'te_icl_ict', 'te_icl_cc', 'te_cc_dc', 'te_cc_hcc', 'te_ie_cc', 'te_m_dc', 'te_m_hcc', 'te_ict_hcc', 'te_ie_hcc', 'te_icl_hcc', 'te_dc_hcc']
    vacinv = ['eag_ve', 'sag_ve', 'hb3_ve', 'mav_ve']
    tot_pars = trtinv + vacinv

    in_df_pars = []
    for par in tot_pars:
        for age in ages:
            for sex in sexes:
                in_df_pars.append(f'{par}_{age}{sex}')

    ## Create a dataframe for the input parameters
    in_df = pd.DataFrame(columns = ['run_id']+in_df_pars)

    for sim in range(1,n_samples+1):
        in_df.loc[sim,'run_id'] = sim
        for par in tot_pars:
            for age in ages:
                for sex in sexes:
                    in_df.loc[sim,f'{par}_{age}{sex}'] = parsets[sim].get_par(par).ts[f'{age}{sex}'].assumption #TODO: get different input parameters per population.loc[sim,col] = parsets[sim].get_par(col).ts[0].assumption #TODO: get different input parameters per population

    return in_df

def central_results(fw = 'hbv_v14_gamma_2.xlsx', db = "AFR_db_v1_2_1.xlsx", cl = 'AFR_calib_v1_2.xlsx'):
    P, res = load_cen_project(fw= fw, db=db, cl=cl)
    ## Locations for the important dataframes
    loc = 'Data/Templates/' # The location for the output templates may be different
    loc2 = 'Data/Demographics/' # The location for the input data may be different (depending on where you store it)

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


    ## Start generating output results
    output_dict = {'cohort_size': 'alive', 'cases': 'tot_inc', 'dalys':'dalys'}
    cen_df = pd.read_csv(loc+'/central-burden-template.202212rfp-1.RFP_standard template_HepB.csv')

    for opt in output_dict:
        tot_val = res.get_variable(output_dict[opt])[0].vals + res.get_variable(output_dict[opt])[1].vals
        for age in range(0,5):
            for year in range(2000,2101):
                ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('0-4',year),'value']
                res_val = tot_val[(year-1990)*4]

                idx = cen_df[(cen_df.year == year) & (cen_df.age == age)].index[0]
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
        tot_val += res.get_variable(var)[0].vals + res.get_variable(var)[1].vals

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

    P, parsets, results = load_sto_project(n_samples, fw= fw, db=db, cl=cl, seed=seed)
    ## TODO: load project for stochastic results

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
                    ratio = float(df1[(df1.age_to == age) & (df1.year == year)].value)/df1.groupby(by=['age_group','year']).sum().loc[('0-4',year),'value']
                    res_val = tot_val[(year-1990)*4]

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