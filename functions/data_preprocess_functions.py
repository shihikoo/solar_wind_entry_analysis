#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  9 20:08:54 2021

@author: jliao
"""

#import pandas as pd
from scipy import stats
import math
import datetime as datetime
import pandas as pd

def extract_date(input_datetime_obj):
    date = input_datetime_obj.strftime("%Y-%m-%d")
    return(date)

def extract_time(input_datetime_obj):
    time = input_datetime_obj.strftime("%H-%M-%S")
    return(time)
    
def extract_dispersion_list(mydata, direction_name = 'PARA'):
    estimated_distance_name = 'ESTIMATED_DISTANCE_' + direction_name
    energy_name = 'EN_' + direction_name
    chisq_name = 'DIS_FITTING_CHISQ_' + direction_name
    dof_name = 'DIS_FITTING_DOF_' + direction_name
    rsquare_name = 'DIS_FITTING_RSQUARE_' + direction_name
    n_dispersion_name = 'N_DISPERSION_' + direction_name
    model_field_length_name = 'FLLEN_' + direction_name
    index = mydata.loc[:,estimated_distance_name].notna()
    mydata = mydata.loc[index,:]
       
    dispersion = mydata.groupby([estimated_distance_name,'date',  n_dispersion_name]).agg({'GSE_X':'count'
                               , chisq_name:'mean', rsquare_name:'mean', dof_name:'mean',  energy_name:'mean', model_field_length_name:'mean'
                               , 'TIME':'mean', 'GSM_X':'mean', 'GSM_Y':'mean', 'GSM_Z':'mean', 'MLT':'median', 'L':'mean',  'STORM_PHASE':'max', 'BX_GSM':'mean'
                               , 'DIST':'mean', 'BETA':'mean', 'datetime_str':'min', 'KP':'mean', 'SW_P':'mean', 'DST':'mean', 'IMF_BY':'mean', 'IMF_BZ':'mean'
                               }).reset_index()
    
    dispersion = dispersion.rename(columns={estimated_distance_name:'estimated_distance', n_dispersion_name:'n_dispersion', 'GSE_X':'dispersion_length',  chisq_name:'chisq', rsquare_name:'rsquare', dof_name:'dof', energy_name:'energy',model_field_length_name:'model_field_line_length_idl'})

    dispersion["direction"] = direction_name

    return(dispersion)

def calculate_cdf(dispersion):
    p = 1-stats.chi2.cdf(dispersion['chisq'], dispersion['dof'])
    return(p)
    
def identify_region(onedata):
    if (onedata['MLT'] >= 8.) & (onedata['MLT'] < 16.):
        region = 'Dayside'
    else:
        if onedata.BETA < 0.05:
            region = 'Lobe'
        elif onedata.BETA < 1.:
            region = 'BL'
        else:
            region = 'PS'
            
    return(region)
    
    
def calculate_velocity(energy, ion_mass = 15.89):
    Avogadro_constant = 6.02214086e23 # moe-1
    electron_charge = 1.60217662e-19  #coulombs
    velocity = math.sqrt(2.*energy*electron_charge/(ion_mass/Avogadro_constant/1e3))/1e3 # 4.577e7*sqrt(data_energy/1e3/AMU) 
    return(velocity)
    
    
def preprocess_data(data):
    cooked_data = data
    cooked_data['datetime_str'] = data.loc[:,'TIME'].apply(datetime.datetime.utcfromtimestamp)
    cooked_data['date'] = data.loc[:,'datetime_str'].apply(extract_date)
    return(cooked_data)

def preprocess_dispersion_list(dispersion_list):
    dispersion_list['p_value'] = dispersion_list.apply(calculate_cdf,axis = 1)
    dispersion_list['region'] = dispersion_list.apply(identify_region, axis=1)

    dispersion_list['index'] = dispersion_list.index

    dispersion_list.loc[ (dispersion_list['GSM_Z'] < 0), 'location'] = 'south'
    dispersion_list.loc[ (dispersion_list['GSM_Z'] > 0), 'location'] = 'north'

    dispersion_list.loc[((dispersion_list['GSM_X'] < -1) & (((dispersion_list['direction'] == 'PARA') & (dispersion_list['BX_GSM'] > 0)) | ((dispersion_list['direction'] == 'ANTI') & (dispersion_list['BX_GSM'] < 0)))) | ((dispersion_list['GSM_X'] > -1) & (((dispersion_list['direction'] == 'ANTI') & (dispersion_list['GSM_Z'] < 0)) | ((dispersion_list['direction'] == 'PARA') & (dispersion_list['GSM_Z'] > 0)))), 'direction_et'] = 'earthward'

    dispersion_list.loc[((dispersion_list['GSM_X'] < -1) & (((dispersion_list['direction'] == 'ANTI') & (dispersion_list['BX_GSM'] > 0)) | ((dispersion_list['direction'] == 'PARA') & (dispersion_list['BX_GSM'] < 0)))) | ((dispersion_list['GSM_X'] > -1) & (((dispersion_list['direction'] == 'PARA') & (dispersion_list['GSM_Z'] < 0)) | ((dispersion_list['direction'] == 'ANTI') & (dispersion_list['GSM_Z'] > 0)))) , 'direction_et'] = 'outward'

    dispersion_list['dispersion_time'] = 2. * (dispersion_list['dof']+2)

    return(dispersion_list)

def extract_dispersions(data, save_to_filename = 'output/dispersion_list.csv'):
    dispersion_para = extract_dispersion_list(data, direction_name = 'PARA')
    dispersion_anti = extract_dispersion_list(data, direction_name = 'ANTI')
    dispersion_list = pd.concat([dispersion_para,dispersion_anti],ignore_index = True)  
    
    dispersion_list = preprocess_dispersion_list(dispersion_list)
    
    dispersion_list.to_csv(save_to_filename)
    
    return(dispersion_list)

