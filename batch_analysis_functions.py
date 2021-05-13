import numpy as np
import pandas as pd

def store_info(df1,df2,info):
    # df1 = df you are copying info to
    # df2 = df you are copying from
    # info = list of strings that you are storing
    for i in range(len(info)):
        df1[info[i]]=np.array([df2[info[i]]]*len(df1))
def print_info(df,info):
    string='-----'
    for i in range(len(info)):
        string+=info[i]+': '+str(row[info[i]])+'__'
    print(string)
    
def normalize_by(norm_var,column_string,target_df,info_df):
    # sort target df by date
    # if not target_df[norm_var].dtypes=='string'
#    target_df=target_df.sort_values(by=norm_var)
    # get all the dates to normalize by an drop duplicates
    norm_column=target_df[norm_var]
    norm_column=norm_column.drop_duplicates()
    # make new column label for normalized variable
    norm_column_label='Normalized '+column_string+' by '+norm_var
    # series to store normalized variables in
    normalized_series=np.array([])
    for x in norm_column:
        df=target_df[target_df[norm_var]==x] # get sub data frame of normalization variable
        variable_column=df[column_string] # get the target variable from that sub data frame
        normalized_variable_column=variable_column/np.max(variable_column) # normalize by dividing by max
        normalized_variable_column=normalized_variable_column.values # take just the values of this normalization
        # concat into array
        normalized_series=np.concatenate((normalized_series,normalized_variable_column))
    # store array as new normalized column in df and return said df
    target_df[norm_column_label]=normalized_series
    return target_df