import cv2
from scipy.ndimage.morphology import binary_fill_holes
import numpy as np
import vsi_metadata as v
import bioformats
import javabridge
import pandas as pd
import skimage.filters as f
import skimage 
import skimage.feature as feature
from skimage import exposure
import matplotlib.pyplot as plt



def single_volume_time_series(path,meta_number=None, cycle_vm=True,
                               crop_img=False,print_number=5,image_channel=1,
                               edge=True,stats=False,show_linear=False,
                               t_lag_level=10,zero_index='None',rescale_factor=255,meta_stage_loop=True,
                               t_sample=1,show=False,rebin=True,compute_height=False,canny_range=.99,
                               hist_cutoff=0.01,manual_thresh=False,m_thresh_low=100
                               ):
    if compute_height:
        area_df,height_df=z_and_t_area(path,meta_number=meta_number, cycle_vm=cycle_vm,
                             crop_img=crop_img,print_number=print_number,image_channel=image_channel,
                             edge=edge,rescale_factor=rescale_factor,meta_stage_loop=meta_stage_loop,
                             t_sample=t_sample,show=show,rebin=rebin,compute_height=compute_height,canny_range=canny_range,
                               hist_cutoff=hist_cutoff,manual_thresh=manual_thresh,
                               m_thresh_low=m_thresh_low)
        
    else:
        area_df=z_and_t_area(path,meta_number=meta_number, cycle_vm=cycle_vm,
                             crop_img=crop_img,print_number=print_number,image_channel=image_channel,
                             edge=edge,rescale_factor=rescale_factor,meta_stage_loop=meta_stage_loop,
                             t_sample=t_sample,show=show,rebin=rebin,compute_height=compute_height,canny_range=canny_range,
                               hist_cutoff=hist_cutoff,manual_thresh=manual_thresh,
                              m_thresh_low=m_thresh_low)
    volume_df=thrombus_volume_time_series(area_df)
    if not zero_index == 'None':
        zero_index=int(zero_index)
        volume_df=volume_df.iloc[zero_index:,:]
        volume_df=volume_df.reset_index()
    volume_df['time (s)']=volume_df['time (s)']-np.min(volume_df['time (s)'])
    if stats:
        vol_metrics=get_v_metrics(volume_df,show_linear=show_linear,t_lag_level=t_lag_level)
        if compute_height:
            return volme_df,vol_metrics,height_df
        else:    
            return volume_df, vol_metrics 
    
    else:
        if compute_height:
            return volume_df,height_df
        else:
            return volume_df

# calculates volume versus time for different segments 
def segmented_volume_time_series(area_df, segments=3):
    z_range = area_df['Z loc (um)'].unique()
    # Segment z-stacks into even groups
    z_segments=segment_series(z_range,segments)
    master_df=pd.DataFrame()
    # iterate through these groups
    for i in range(np.shape(z_segments)[0]):
        # isolate the z-stacks for given loop
        z_segments_loop=z_segments[i,:]
        # find parts of area df corresponding to this z-range
        loop_area_df=area_df[(area_df['Z loc (um)']>=z_segments_loop[0])&(area_df['Z loc (um)']<=z_segments_loop[1])]
        # make a string to store the z-range we are looking at 
        z_range_string='{:.2f}-{:.2f}'.format(*z_segments_loop)
        # take volume integral 
        loop_volume_df=thrombus_volume_time_series(loop_area_df)
        # store the z-range and segments 
        loop_volume_df['Z-range (um)']=[z_range_string]*len(loop_volume_df)
        loop_volume_df['Z-segment']=[i+1]*len(loop_volume_df)
        # append all of this to a bigger df and return this bad boy
        master_df=master_df.append(loop_volume_df,ignore_index=True)
    return master_df

# function for making segments out of the z stacks to analyze different regions 
# of the data
def segment_series(series,segments):
    # calculate the space between the different series for a certain number of 
    # segments
    iterator=int(np.rint(len(series)/segments))
    segment_list=[]
    # iterate through the space by stepping though with the iterator
    for i in range(0,len(series)-iterator,iterator):
        # If the data doesn't divide into even segments we will need to handle 
        # this by truncating the last entry
        if i+iterator>len(series) or i>len(series):
            segment_list.append([series[i],series[-1]])
        # Otherwise evenly truncate by segmenting at i+ iterator
        else:
            segment_list.append([series[i],series[i+iterator]])
    segment_list=np.array(segment_list)
    return segment_list

# calculate the thrombus volume as a function of time
def thrombus_volume_time_series(area_df):
    # get the nonrepeating timesteps
    time = area_df['time (s)'].unique()
    time_step_length = (np.max(time) - np.min(time)) / len(time)
    # get nonrepeating z steps
    z_range = area_df['Z loc (um)'].unique()
    # make an interpolated z series 
    z_range_fine = np.linspace(np.min(z_range), np.max(z_range), 50)
    # calculate z_step
    z_step = (np.max(z_range_fine) - np.min(z_range_fine)) / len(z_range_fine)
    # make an empty volume array to be filled
    volume = np.zeros(len(time))
    dV = np.empty(len(time))
    dV[:] = np.nan
    count = 0
    for t in time:
        time_step=area_df[area_df['time (s)']==t]
        # interpolate the z series for a finer integration
        area_fine=interpolate_in_z(time_step,'Thrombus Area (um^2)',50)
        # take the sum of this and multiply by dZ to get the volume integral
        sum_area=np.sum(area_fine)
        volume[count]=sum_area*z_step
        # calculate the volumetric growth rate (dV/dt)
        if count>1:
            dV[count-1]=(volume[count]-volume[count-1])/time_step_length
        count+=1
    # store values in dataframe
    volume_df=pd.DataFrame({'time (s)':time,'Thrombus Volume (um^3)':volume,'dV/dt':dV})
    return volume_df

# function to interpolate the z data more finely to more accurately calculate 
# volume of the images 
def interpolate_in_z(df,variable,steps):
    z_range = df['Z loc (um)']
    z_range_fine = np.linspace(np.min(z_range), np.max(z_range), steps)
    var = df[variable].values
    if len(var)==0:
        var_fine=np.zeros(len(z_range_fine))
    else:
        var_fine = np.interp(z_range_fine, z_range, var)
    return var_fine


# combines in height profile into loops of z and t
def z_and_t_area(path,meta_number=None, cycle_vm=True,crop_img=False,print_number=5,edge=False,binary=False,
                 image_channel=1,rescale_factor=255,meta_stage_loop=True,t_sample=1,show=True,
                 rebin=True,compute_height=False,canny_range=.33,
                               hist_cutoff=0.01,manual_thresh=False,m_thresh_high=1000,m_thresh_low=100):
    
    # start javabridge
    if cycle_vm:
        javabridge.start_vm(class_path=bioformats.JARS)
    # read in metadata using bioformats and make ararys for t ans z slices
    metadata=v.extract_metadata(path,cycle_vm=False,meta_number=meta_number,stage_loop=meta_stage_loop)
    z_slices=np.arange(0,metadata['size_Z'])
    z_slices_scaled=z_slices*float(metadata['relative step width'])
    t_slices=np.arange(0,metadata['size_T']-1)
    t_slices=t_slices[::t_sample]
    t_slices_scaled=t_slices*float(metadata['cycle time'])
    # declare dataframes to store analysis results in
    area_df=pd.DataFrame()
    if compute_height:
        height_df=pd.DataFrame()
    # crop image option to do if specified
    if crop_img:
        # determine cropping bounds
        img=bioformats.load_image(path=path,z=z_slices[3],t=t_slices[-1],series=0)
        img = (img * rescale_factor).astype('uint8')
        bounds = get_bounds(img,scale=0.4)
    print_thresh=np.linspace(0,100,print_number)
    print_count=0
    plot_grid_count=1
    finished_plotting=False
    # iterate through time slices, stop before the last one
    for i in range(len(t_slices)-1):
        percent=t_slices[i]/(len(t_slices)-1)
        if percent>print_thresh[print_count]:
            # print('\tt: {:.2f} S'.format(t_slices_scaled[i]))
            print_count+=1
        # declare surface_area array outside of z-stack loop
        area=np.zeros(len(z_slices))
        for j in range(len(z_slices)):
            # read image
            img=bioformats.load_image(path=path,z=z_slices[j],t=t_slices[i],
                                      rescale=False)
            # if the image consists of multiple channels take just one channel 
            if len(np.shape(img))>2:
                img=img[:,:,image_channel]
            equalize_hist=exposure.equalize_adapthist(img,clip_limit=hist_cutoff,nbins=65535)
            ratio=np.amax(equalize_hist)/65535
            equalize_hist=np.round(equalize_hist/ratio,0).astype(int)
            
            # if rebin is specified rebin image using 2x2 binning
            if rebin:
                equalize_hist=bin2d(img,2)
            
            if crop_img:
                crop=equalize_hist[int(bounds[1]):int(bounds[1] + bounds[3]), int(bounds[0]):int(bounds[0] + bounds[2])]
            else:
                crop=equalize_hist
            # threshold image
            if edge:
                bw=edge_and_fill(crop,canny_range=canny_range)
            else:
                if manual_thresh:
                    bw=thresh_and_fill(img,low=m_thresh_low)
                else:
                    bw=auto_threshold(img,'triangle')
            # if compute height specified and on the firt step of loop make a zeros array
            # the length and width of the image and height of the z_stacks number
            if compute_height:
                if j==0:
                    img_shape=np.shape(img)
                    img_stack=np.zeros([img_shape[0],img_shape[1],len(z_slices)])
                img_stack[:,:,j]=bw
            # subpart of program to show images comparing plots andbin
            if show:
                f1=plt.figure()
                z_count=0
                t_count=0
                count_satisfied,t_count,z_count=show_handler(i,len(t_slices),
                                                             j,len(z_slices),
                                                             t_count,z_count)
                if count_satisfied and not finished_plotting:
                    plt.subplot(6,6,plot_grid_count)
                    plt.imshow(equalize_hist)
                    plt.axis('off')
                    plt.subplot(6,6,plot_grid_count+1)
                    plt.imshow(bw)
                    plt.axis('off')
                    plot_grid_count+=2
                    if plot_grid_count>=36:
                        finished_plotting=True
                f1.show()
            # calculate area of thresholded image
            if rebin:
                scale=metadata['scale']**2
            else:
                scale=metadata['scale']**2
            area[j]=get_area(bw,scale)
            if compute_height:
                average_height,length=calc_height(img_stack,scale,metadata['relative step width'])
                loop_height=pd.DataFrame({
                    'Length (um)':length,
                    'Height (um)':average_height,
                    'time (s)':[t_slices_scaled[i]]*len(length)
                    })
                height_df=height_df.append(loop_height,ignore_index=True)
                
        # store in datafame along with other imporant variables
        loop_area_df=pd.DataFrame({'Thrombus Area (um^2)':area,
                                   'Z loc (um)':z_slices_scaled,
                                   'time (s)':[t_slices_scaled[i]]*len(area)
                                   })
        area_df=area_df.append(loop_area_df,ignore_index=True)
    if cycle_vm:
        javabridge.kill_vm()
    if compute_height:
        return area_df,height_df
    else:
        return area_df

# compute height as a func of length 
def calc_height(img_stack,scale,z_scale):
    dims=np.shape(img_stack)
    height=np.zeros((dims[0],dims[1]))
    avg_height=np.zeros(dims[0])
    for i in range(dims[0]):
        for j in range(dims[1]):
            img_slice=img_stack[i,j,:]
            height[i,j]=np.max(np.where(img_slice))*z_scale
    avg_height=np.mean(height,axis=0)
    total_length=dims[0]*scale
    length_array=np.linspace(0,total_length,dims[0])
    return avg_height,length_array
            
def bin2d(a,K):
    m_bins = a.shape[0]//K
    n_bins = a.shape[1]//K
    new_img=a.reshape(m_bins, K, n_bins, K).sum(3).sum(1)   
    return new_img

def show_handler(i,t_slices,j,z_slices,t_count,z_count,
                 z_regions=[1,50,100],t_points=[25,50,90]):
    count_satisfied=False
    t_percent=i/t_slices*100
    z_percent=j/z_slices*100
    if t_percent>t_points[t_count] and z_percent>z_regions[z_count]:
        count_satisfied=True
        t_count+=1
        z_count+=1
    return count_satisfied, t_count,z_count
    
    

# function to get calculat the area of a thresholded image 
def get_area(bw,scale):
    # get contours from thresholded image
    # bw=bw.astype('uint8')
    # contours, _ = cv2.findContours(bw, cv2.RETR_TREE, cv2.CHAIN_APPROX_NONE)
    # # calculate the surface area of each contour
    # area= [cv2.contourArea(c)*scale**2 for c in contours]
    # take the sum and return
    area_tot=np.sum(bw)*scale
    
    return area_tot

# autothresholding function
# - works by passing in a bw image and the filter_type
# - Filter options are seen below in the if statements
def auto_threshold(image,filter_type):
    if filter_type=='otsu':
        filt= f.threshold_otsu(image) 
        bw= image>filt
        bw=bw.astype(int)
    elif filter_type=='isodata':
        filt= f.threshold_isodata(image) 
        bw= image>filt
        bw=bw.astype(int)
    elif filter_type=='triangle':
        filt= f.threshold_triangle(image) 
        bw= image>filt
        bw=bw.astype(int)
    elif filter_type=='entropy':
        filt= f.threshold_li(image) 
        bw= image>filt
        bw=bw.astype(int)
    elif filter_type=='mixed':
        snr=np.mean(image)/np.std(image)
        if snr>2:
             filt= f.threshold_triangle(image) 
        else:
            filt= f.threshold_otsu(image) 
        bw= image>filt
        bw=bw.astype(int)     
    return bw

# Manual thresholding function
def thresh_and_fill(img,low=125):
    low_pass=img>low
    # high_pass=low_pass<high
    # fill holes
    filled=binary_fill_holes(low_pass)
    # cast as 8 bit agin
    filled = np.array(filled,dtype=np.uint8)
    return filled

def edge_and_fill(img,canny_range=0.99):
    edges=auto_canny(img,sigma=canny_range)
    # Dilate edges to make unconnected lines connected
    # Kernel to dialiate the image later on and fill the holes. Declared outside of the loop to save on processing power
    kernel = np.ones((3, 3), np.uint8)
    dilate = skimage.morphology.dilation(edges, selem=kernel)
    # Fill them holes
    fill = binary_fill_holes(dilate).astype(int)
    # Lines of code I found online to fill in stuff around the edges better
    labels = skimage.morphology.label(fill)
    labelCount = np.bincount(labels.ravel())
    background = np.argmax(labelCount)
    fill[labels != background] = 1
    return fill
    
def auto_canny(image, sigma=0.33,std_eval=False,correct_fact=1):
    float_img=np.float32(image)
    blurred = cv2.GaussianBlur(float_img, (3, 3),0)
    
    #compute the median of the single channel pixel intensities
    v = np.mean(blurred)
    std=np.std(blurred)
    mean_calc=v
    std_calc=std
    snr=mean_calc/std_calc
    # if snr>3:
    #     correct_fact=1
	# apply automatic Canny edge detection using the computed median
    if std_eval:
        lower = int((v-2*std)*correct_fact)
        upper = int((v+2*std)*correct_fact*correct_fact)
    else:
        lower = int(max(0, (1.0 - sigma) * v))
        upper = int(min(65535, (1.0 + sigma) * v))
    edged =feature.canny(blurred, low_threshold=lower, high_threshold=upper)
    return edged  


# Function that uses crop function to get bounds of the image
def get_bounds(img,scale=0.4):
    # resize image so it fits on screen when cropping
    dims = list(np.shape(img))
    dims = [float(d) for d in dims]
    adj_dims = [int(np.rint(scale * d)) for d in dims]
    adj_dims = tuple(adj_dims)
    img_scale = cv2.resize(img, dsize=adj_dims, interpolation=cv2.INTER_CUBIC)
    # crop image and return the bounds
    bounds = crop_img(img_scale, get_bounds=True)
    bounds=list(bounds) # convert to list for operations
    # rescale bounds based on rescale factor 
    for i in range(len(bounds)):
        bounds [i] = int(np.rint(bounds[i]/scale))
    return bounds
 
# Function for cropping image 
def crop_img(img,get_bounds=False,scale=1):
    # Select ROI
    r = cv2.selectROI(img)
    cv2.destroyAllWindows()
    if get_bounds==True:
        return r
    else:
        # Crop image
        im_crop = img[int(r[1]):int(r[1] + r[3]), int(r[0]):int(r[0] + r[2])]
        return im_crop
    
def get_v_metrics(df,show_linear=False,max_t_lag=5*60,t_lag_level=10):
    # get min and max surface coverages
    maxSC = np.max(df['Thrombus Volume (um^3)'])
    minSC = np.min(df['Thrombus Volume (um^3)']) + 0.000001
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
        f2=plt.figure(2)
        y_linear = x[t_lag_index:max_index] * slope + intercept
        plt.plot(x[t_lag_index:max_index], y_linear, '-r', label='Linear Fit', linewidth=2)
        plt.plot(df['time (s)'], df['Thrombus Volume (um^3)'], 'ob', label='Experimental Data',alpha=0.5)
        plt.axvline(x=t_lag, color='k', linestyle='--', label='T-lag')
        plt.axhline(y=maxSC, color='g', linestyle='--', label='Max')
        plt.xlabel('time (s)')
        plt.ylabel('Thrombus Volume (um^3)')
        plt.title('Experimental Data and Fitted Kinetic Metrics')
        plt.legend()
        plt.show()
    # Store data using pandas 
    SC_values = pd.DataFrame([{'V T-lag': t_lag, 'V Max': maxSC, 'V Slope': slope}])
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
    try:
        coefs[1]
    except:
        intercept=0
    else:
        intercept=coefs[1]
    return slope,intercept

def max_slope(x,y,t_lag_index,max_index,iterator,min_fit_length=25):
    x_resample=np.arange(t_lag_index,max_index)
    y_resample=np.interp(x_resample,x,y)
    grad=np.gradient(y_resample,x_resample)
    slope=np.max(grad)
    intercept=y[t_lag_index]
    coefs=np.array([slope,intercept])
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

# Get index for first data point 2% away from the max
def get_lag_time(x,y,t_lag_level,max_t_lag,max_index,maxSC):
    y=y-y[0]
    t_lag_index=None
    for i in range(len(y)):
        # If the surface coverage is greater or equal 5% then the time at that point is stored and the loop ends
        if round(y[i], 2) >= t_lag_level:
            t_lag = x[i]
            t_lag_index = i
            break
    # if a lag time wasn't picked up auto assing the maximum
    if t_lag_index == None:
        t_lag=max_t_lag
        t_lag_index=0
    # If the lag time is huge, change the index we pass the slope fitter
    if t_lag_index>0.6*len(x):
        t_lag_index=0

    return t_lag,t_lag_index,max_index

# Get index for first data point 20% away from the max
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
    interval=df['time (s)'][1]-df['time (s)'][0]
    start=np.min(df['time (s)'])
    stop=np.max(df['time (s)'])+interval
    x=np.arange(start,stop,1)
    y=np.interp(x,df['time (s)'],df['Thrombus Volume (um^3)'])
    return x,y
