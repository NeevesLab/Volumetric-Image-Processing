import pandas as pd
import numpy as np
import os
import glob
from skimage import io 
import matplotlib.pyplot as plt
import javabridge
import bioformats
import vsi_metadata as v


def fluorescence_time_series (filepath,interval=18,threshold=100,
                               csv_path='',stats_path='',store_csv=False,
                               zero_index_time=0,stats=False,show_linear=False,
                               t_lag_level=250,rescale='None',background='None',
                               zero_index='None',threshold_filter=True,vsi=False,
                               cycle_vm=True,meta_number=None,image_channel=1,meta_stage_loop=True,
                               t_sample=1
                               ):
    
    # if a vsi file is specified then read in through vsi means rather than manually reading in tifs
    if vsi:
        # start javabridge
        if cycle_vm:
            javabridge.start_vm(class_path=bioformats.JARS)
        # read in metadata using bioformats and make ararys for t ans z slices
        metadata=v.extract_metadata(filepath,cycle_vm=False,meta_number=meta_number,stage_loop=meta_stage_loop)
        t_slices=np.arange(0,metadata['size_T']-1)
        t_slices=t_slices[::t_sample]
        t_slices_scaled=t_slices*float(metadata['cycle time'])
        z_slices=np.arange(0,metadata['size_Z'])
        mean = np.empty(len(t_slices))
        minimum = np.empty(len(t_slices))
        maximum = np.empty(len(t_slices))
        for i in range(len(t_slices)):
            # show_handler,not_shown,pic_pct=show_controller(i+1,pic_length,*not_shown)
            # read in image and convert to 8 bit
            if len(z_slices)==1:
                image=bioformats.load_image(path=filepath,t=t_slices[i],series=0,rescale=False)
                if len(np.shape(image)>2):
                    image=image[:,:,image_channel]
            else:
                image=max_projection(filepath,t_slices[i],z_slices,image_channel,rescale=False)
            mean[i] = np.mean(image)
            minimum[i] = np.min(image)
            maximum[i] = np.max(image)
    # Otherwise manually read in tifs 
    else: 
        # Boolean process to deterine whether or not tor
        former_path=os.getcwd()
        os.chdir(filepath)
        # Read sort in the tif files of a given pathway
        filenames=sorted(glob.glob('*.tif'))
        # read in the images and get the mean, min, and max then store as dataframe
        mean = np.empty(len(filenames))
        minimum = np.empty(len(filenames))
        maximum = np.empty(len(filenames))
        i=0
        for file in filenames:
            image = io.imread(file)
            mean[i] = np.mean(image)
            minimum[i] = np.min(image)
            maximum[i] = np.max(image)
            i += 1
        os.chdir(former_path)
    # Sometimes it makes sense to specify a index (time point) to subtract background by specifying an image. If no index is specified
    # then I make an autodetect function using a threshold filter
    if zero_index == 'None':
        if threshold_filter:
            zero_index = 0
            for i in range(len(mean)):
                if i == 0:
                    continue
                    if mean[i] > threshold:
                        zero_index = i
                        break
    else:
        zero_index=int(zero_index)
    # cut arrays using zero index value
    mean = mean[zero_index:]
    minimum = minimum[zero_index:]
    maximum = maximum[zero_index:]
    # rescale mean of image if specified
    if rescale!='None':
        mean=mean*rescale
        maximum=maximum*rescale
        minimum=minimum*rescale
    # subtract background. If a specific background in inputted and spec_background = True then it will use the specified value
    if background=='None':
        zero_mean = mean - mean[0]
    else:
        zero_mean = mean - background
    # generate time series data
    if vsi:
            time=t_slices_scaled[zero_index:]
    time = np.linspace(zero_index_time, interval * len(mean), len(mean),endpoint=False)
    # store to dataframe for easy usage
    df = pd.DataFrame({'Time (s)': time, 'Mean': mean, 'Zero Mean': zero_mean, 'Min': minimum, 'Max': maximum})
    # delte rows with saturated values
    df=df[df['Mean']<60000]
    if stats:
        f_metric_options=dict(show_linear=show_linear,interval=interval
                )
        try:
            t_lag_level
        except:
            pass
        else:
            f_metric_options['t_lag_level']=t_lag_level
        F_values = get_F_metrics(df, **f_metric_options)
    
    if stats:
        return df, F_values
    else:
        return df

def max_projection(path,t_slice,z_slices,image_channel):
    count=0
    for z in z_slices:
        if z==np.min(z_slices):
            test_img = bioformats.load_image(path=path, t=t_slice,z=z, series=0)
            shape=np.shape(test_img)
            img_array=np.zeros([shape[0],shape[1],len(z_slices)])
        img_loop=bioformats.load_image(path=path, t=t_slice,z=z, series=0)*65535
        if len(np.shape(img_loop))>2:
            img_loop=img_loop[:,:,image_channel]
        img_array[:,:,count]=img_loop
        count+=1
    max_intensity=np.amax(img_array,axis=2)
    return max_intensity

# controller that tells the return_area function to show a comparison of the thresholded function at set values in the 
# pic thresh list
# controller that tells the return_area function to show a comparison of the thresholded function at set values in the 
# pic thresh list
# controller that tells the return_area function to show a comparison of the thresholded function at set values in the 
# pic thresh list
def show_controller(count_pics,pic_length,*not_shown,pic_thresh=[20,50,90]):
    # cast not shown as list
    not_shown=list(not_shown)
    pic_pct = np.round(count_pics / pic_length * 100, 2)
    if pic_pct >= pic_thresh[0] and not_shown[0] == True:
        show_handler = True
        not_shown[0] = False
    elif pic_pct >= pic_thresh[1] and not_shown[1] == True:
        show_handler = True
        not_shown[1] = False
    elif pic_pct >= pic_thresh[2] and not_shown[2] == True:
        show_handler = True
        not_shown[2] = False
    else:
        show_handler = False
    return show_handler,not_shown,pic_pct

def get_F_metrics(df,show_linear=False,max_t_lag=5*60,t_lag_level=200,interval=1):
    # get min and max surface coverages
    maxSC = np.max(df['Zero Mean'])
    minSC = np.min(df['Zero Mean']) + 0.000001
    # resample series to finer grain for better analysis 
    x,y=interpolate_series(df)
    # get index of 80% val of max
    max_index=get_max_index(x,y)
    # get lag time and update certain indices if need be
    t_lag,t_lag_index,max_index=get_lag_time(x,y,t_lag_level,max_t_lag,max_index,maxSC)
    # Get slope of linear region
    slope,intercept=get_slope(x,y,t_lag_index,max_index)
    if show_linear:
        # Plot raw data 
        plt.figure(1)
        y_linear = x[t_lag_index:max_index] * slope + intercept
        plt.plot(x[t_lag_index:max_index], y_linear, '-r', label='Linear Fit', linewidth=2)
        plt.plot(df['Time (s)'], df['Zero Mean'], 'ob', label='Experimental Data',alpha=0.5)
        plt.axvline(x=t_lag, color='k', linestyle='--', label='T-lag')
        plt.axhline(y=maxSC, color='g', linestyle='--', label='Max')
        plt.xlabel('Time (s)')
        plt.ylabel('Fluorescence Intensity')
        plt.title('Experimental Data and Fitted Kinetic Metrics')
        plt.legend()
        plt.show()
    # Store data using pandas 
    SC_values = pd.DataFrame([{'F T-lag': t_lag, 'F Max': maxSC, 'F Slope': slope}])
    return SC_values

def get_slope(x,y,t_lag_index,max_index,metric='max_slope',min_fit=0.5):
    min_fit_length = int(round(min_fit*(len(x)+1)))
    if max_index>=100:
        iterator=int(round(len(x)/100))
    else:
        iterator=1
    # Now the code chooses the slope on the criteria of either:
    # - getting the largest possible slope of a line that is min_fit*the fitting length
    if metric=='max_slope':
        coefs=max_slope(x,y,t_lag_index,max_index,iterator,min_fit_length)
    # - or by getting the most linear region for the fitting space
    if metric=='best_fit':
        # Search on the interval from t-lag to 2% from the max value
        coefs=chi_squared_min(x,y,t_lag_index,max_index,iterator,min_fit_length)
    slope=coefs[0]
    # no negative slopes!
    if slope<0: slope =0
    intercept=coefs[1]
    return slope,intercept

def max_slope(x,y,t_lag_index,max_index,iterator,min_fit_length):
    slope_best=0
    i_best=t_lag_index
    j_best=max_index
    j_best=max_index
    for i in range(t_lag_index, max_index ,iterator):
        for j in range(i+min_fit_length, max_index,iterator):
            # If the array is zero 
            if not x[i:j].size or not y[i:j].size:
                print ('Array Empty in Loop')
                print (i,j)
            elif not np.absolute(j-i)<min_fit_length:
                coefs_loop = np.polyfit(x[i:j], y[i:j], 1)
                slope_loop=coefs_loop[0]
                if slope_loop>slope_best:
                    i_best = i
                    j_best = j
                    slope_best=slope_loop
    if not x[i_best:j_best].size or not y[i_best:j_best].size:
        print('Array Empty Outside of Loop')
        print(i_best,j_best)
        coefs=[float('nan')]
    else:
          coefs = np.polyfit(x[i_best:j_best], y[i_best:j_best], 1)
    return coefs

def chi_squared_min(x,y,t_lag_index,max_index,iterator,min_fit_length):
    chi_min=1E-3
    for i in range(t_lag_index, max_index ,iterator):
        for j in range(i+min_fit_length, max_index,iterator):
            # If the array is zero 
            if not x[i:j].size or not y[i:j].size:
                print ('Array Empty in Loop')
                print (i,j)
            elif not np.absolute(j-i)<min_fit_length:
                coefs_loop = np.polyfit(x[i:j], y[i:j], 1)
                y_linear = x * coefs_loop[0] + coefs_loop[1]
                chi = 0
                for k in range(i, j):
                    chi += (y_linear[k] - y[k]) ** 2

                if chi < chi_min:
                    i_best = i
                    j_best = j
                    chi_min = chi
                # print 'Chi-min: '+str(chi_min)
                # print 'Chi:'+str(chi)
    if not x[i_best:j_best].size or not y[i_best:j_best].size:
        print('Array Empty Outside of Loop')
        print(i_best,j_best)
        coefs=[float('nan')]
    else:
        coefs = np.polyfit(x[i_best:j_best], y[i_best:j_best], 1)
    return coefs

def get_lag_time(x,y,t_lag_level,max_t_lag,max_index,max_SC):
    for i in range(len(y)):
        # If the surface coverage is greater or equal 5% then the time at that point is stored and the loop ends
        if round(y[i], 2) >= t_lag_level:
            t_lag = x[i]
            t_lag_index = i
            break
    # If the lag time is huge, change the index we pass the slope fitter
    try:
        t_lag_index
    except:
        t_lag_index=0
        t_lag=max_t_lag
    else:
        if t_lag_index>0.6*len(x):
            t_lag_index=0

    return t_lag,t_lag_index,max_index

# Get index for first data point 2% away from the max
def get_max_index(x,y,tolerance=0.2):
    max_index=-1
    max_SC=np.max(y)
    for i in range(len(y)):
        difference = np.absolute((y[i] - max_SC) / max_SC)
        if difference <= tolerance:
            max_index = i
            break
    if max_SC<100:
        max_index=len(y)+1
    return max_index

def interpolate_series(df,interp_interval=1):
    # Get index
    interval=df['Time (s)'][1]-df['Time (s)'][0]
    start=np.min(df['Time (s)'])
    stop=np.max(df['Time (s)'])+interval
    x=np.arange(start,stop,1)
    y=np.interp(x,df['Time (s)'],df['Zero Mean'])
    return x,y
def interpolate_time_series(df,interp_var=['Zero Mean'],x_var='Time (s)',interval=1):
    a_id=df['Assay ID'].unique()
    master_df=pd.DataFrame()
    for i in range(len(a_id)):
        new_indv_assay=pd.DataFrame()
        indv_assay=df[df['Assay ID']==a_id[i]]
        x=indv_assay[x_var].values
        start=np.nanmin(x)
        stop=np.nanmax(x)
        # step=x[1]-x[0]
        # length=int(np.round((stop-start)/interval,0))
        x_new=np.arange(start,stop,interval)
        new_indv_assay[x_var]=x_new
        for j in range(len(interp_var)):
            y=indv_assay[interp_var[j]]
            y_new=np.interp(x_new,x,y)
            new_indv_assay[interp_var[j]]=y_new
        for c in indv_assay.columns:
            if not c in interp_var and c!=x_var:
                info_value=indv_assay[c].iloc[0]
                new_indv_assay[c]=[info_value]*len(new_indv_assay)
        master_df=master_df.append(new_indv_assay,ignore_index=True)
    return master_df
            
            
def pic_sort(filepath,key,channel_key='C000'):
    # Filepath = the path where all the pics are stored
    # Key = list of Strings that you want the subfolders to be called. The order should be in the same order
    # that the files are found in the folder. For example, if I had a folder that showed a brightfield picture first,
    # then FITC, then TRITC, my key list would be ['BF','FITC','TRITC'}
    former_path=os.getcwd()
    os.chdir(filepath)
    
    # Read sort in the tif files of a given pathway
    filenames = sorted(glob.glob('*.tif'))
    # Makes sub directories based on the keys that you put in
    dir=[filepath+'/'+x for x in key]
    do_sort=False
    # Checks to see that there are tif files to be sorted
    for fname in os.listdir(filepath):
        if fname.endswith('.tif'):
            do_sort=True
            break
    # If tiffs are in the folder, then the program will do the sorting
    if do_sort==True:
        # Generates the directories if they haven't been made yet
        for i in range(len(dir)):
            if os.path.isdir(dir[i])==False:
                os.mkdir(dir[i])
        # Goes through the filenames and sorts them by order
        for file in filenames:
            for i in range(len(key)):
                tag=channel_key+str(i+1)
                if tag in file:
                    new_path=dir[i]+'/'+os.path.basename(file)
                    os.rename(file,new_path)
    os.chdir(former_path)

def pic_sort_3d(filepath,key=['BF','TRITC'],t_key='T00'):
    former_path=os.getcwd()
    os.chdir(filepath)
    # determine the number of time steps:
    count=1 # counter for t slices
    # first sort pictures for each channel
    key_names=[]
    for k in key:
        loop_key_name=filepath+'/'+k
        if not os.path.isdir(loop_key_name):
            os.mkdir(loop_key_name)
        key_names.append(loop_key_name)
    pic_sort(filepath,key,channel_key='C00')
#     sort pictures for each time step
    for k in key_names:
        os.chdir(k)
        print(k)
        sub_f_names=sorted(glob.glob('*.tif'))
        # -------- go through folder and determint the number of time steps-------
        find_t=True
        t_string_list=[]
        while(find_t):
            # string that we are looking for in the filename
            t_string=t_key+str(count)
            # check if t_string is in file name: if so add to the counter
            if any(t_string in f for f in sub_f_names):
                t_string_list.append(t_string)
                count+=1
            # otherwise exit the loop
            else:
                find_t=False
        # ----- Now loop through different time steps 
        for t in t_string_list:
            print('\t'+str(t))
            t_path=k+'/'+t
            if not os.path.isdir(t_path):
                os.mkdir(t_path)
            for f in sub_f_names:
                if t in f:
                    newpath=t_path+'/'+os.path.basename(f)
                    os.rename(f,newpath)
    os.chdir(former_path)        