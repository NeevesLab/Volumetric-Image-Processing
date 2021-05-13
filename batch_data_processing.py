import pandas as pd
import surface_coverage_processing as s
import fluorescence_processing as f
import volumetric_processing as v
import numpy as np
import os
import bioformats
import javabridge
import batch_analysis_functions as b
import re

# read in data from spreadsheet
sheet_path="G:\My Drive\Grad School\Lab\Segmented Injury and Pressure Drop Paper\Homogenous Vs Heterogenous - Peptides & TF - No Peptides in TF Strip\Flow Assay Log.xlsx"
sheet_path=sheet_path.replace('\\','/')
master_info_df=pd.read_excel(sheet_path,sheet_name='Kinetic Data')
dynamic_master=pd.DataFrame()
metric_master=pd.DataFrame()
old_dynamic_master=pd.read_csv('time_series_data.csv')
old_metric_master=pd.read_csv('metrics_data.csv')

# -------info to copy from into df to results of analysis
info=['Date','Donor','Donor No','dP (kPa)',
      'Surface','Channel Number','Assay ID']

current_path=os.getcwd()
javabridge.start_vm(class_path=bioformats.JARS)
# A for loop to loop through the cells in the excel sheet where I store all of the relevant information for the assays
for index, row in master_info_df.iterrows():
    kinetic_path = row['File Path']
    print('---------Row: ', str(index+2),' ',kinetic_path)
    # Only do analysis for cells where I specify a file location since some cells don't have this
    if not pd.isna(kinetic_path):
        # since a relative path is passed in from the excel sheet we need to add the absolute path 
        # to that 
        # from Cell Sense they dump the channels of all photos into the same folder but distinguish them by file type
        kinetic_path=kinetic_path.replace('\\','/')
        path= kinetic_path
        # We have an option in the excel sheet to pull the data from the old master csvs for a given assay if we don't need to analyze it
        if row['Analyze Assay']=='N':
           print('\t','Read in From Previous Analysis')
           dynamic_df=old_dynamic_master[old_dynamic_master['Assay ID']==row['Assay ID']]
           met_df=old_metric_master[old_metric_master['Assay ID']==row['Assay ID']]
        else:
            # Generate a dictionary that allows us to specify the options we want to pass into our function to evaluate mean fluorescence and metrics
            options = dict(
                             stats=True,  # argument telling function to return summary statistics
                             show_linear=True, # argument telling the function to show show the summary stats plotted over the raw data so we can verify they are
                             # being evaluated properly
                             show=False, # show the analyzed images
                             t_lag_level=200E3,
                             cycle_vm=False,
                             image_channel=1,
                             meta_stage_loop=False,
                             t_sample=1,
                             rebin=False
                             )
            # Optional things that are specified for the functions
            if not pd.isna(row['Zero Index']):
                options['zero_index']=row['Zero Index']
            # if a background is specified we feed that into the fluorescence processing 
            if not pd.isna(row['Background']):
                options['background']=row['Background']
            if row['Stage Loop']=='N':
                options['meta_stage_loop']=False
            if not pd.isna(row['Meta Number']):
                options['meta_number']=int(row['Meta Number'])
            if not pd.isna(row['Zero Index']):
                options['zero_index']=row['Zero Index']
            
            dio_vol,dio_m=v.single_volume_time_series(path,**options)
            dio_m.columns = ['DIO'+ str(col) for col in dio_m.columns]
            for col in dio_vol.columns:
                if col !='time (s)':
                    dio_vol=dio_vol.rename(columns={str(col):'DIO'+str(col)})
            options['image_channel']=0
            options['edge']=False
            options['manual_thresh']=True
            options['m_thresh_low']=110
            # options['hist_cutoff']=.001
            fib_vol,fib_m=v.single_volume_time_series(path,**options)
            fib_m.columns = ['FIB'+ str(col) for col in fib_m.columns]
            fib_vol=fib_vol.drop(columns=['time (s)'])
            fib_vol.columns=['FIB'+ str(col) for col in fib_vol.columns]
            
            digits=re.findall(r'\d+', path)
            last_digit=int(digits[-1])
            path=path.replace('0'+str(last_digit),'0'+str(last_digit+1))
            ps_vol,ps_m=v.single_volume_time_series(path,**options)
            ps_m.columns = ['PS'+ str(col) for col in ps_m.columns]
            ps_vol=ps_vol.drop(columns=['time (s)'])
            ps_vol.columns=['PS'+ str(col) for col in ps_vol.columns]
            dynamic_df=pd.concat([dio_vol,fib_vol,ps_vol],axis=1)
            
            met_df=pd.concat([dio_m,fib_m,ps_m],axis=1)
            # met_df=fl_m
            # store info about assay conditions into analsysi
            b.store_info(dynamic_df,row,info)
            b.store_info(met_df,row,info)
        # Store all of this data in dataframes so we can easily plot all of it together
        dynamic_master=dynamic_master.append(dynamic_df,ignore_index=True)
        metric_master=metric_master.append(met_df,ignore_index=True)
    dynamic_master.to_csv('time_series_data.csv',index=False)
    metric_master.to_csv('metrics_data.csv',index=False)
        # metric_platelet_df=metric_platelet_df.append(m,ignore_index=True)        
javabridge.kill_vm()
# ----------------Normalize data to max of that day
# Function to do normalization with

# Variables I want to normalize by
# normalization_variables=['Donor No']
# for i in range(len(normalization_variables)):
#     # Perform normalization on all of the storage dfs
#     dynamic_df=normalize_by(normalization_variables[i],'Surface Coverage',dynamic_df,master_info_df)
#     metric_platelet_df=b.normalize_by(normalization_variables[i],'Max',metric_platelet_df,master_info_df)
#     metric_platelet_df=b.normalize_by(normalization_variables[i],'Slope',metric_platelet_df,master_info_df)



# metric_platelet_df.to_csv('metrics_data.csv',index=False)