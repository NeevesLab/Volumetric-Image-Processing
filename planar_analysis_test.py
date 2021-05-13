import planar_analysis as p
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
path='F:/2021-2-13 Het Patterning\Homo Patterning Assay 2 - Homo-1kpa_01.vsi'
path=path.replace('\\','/')
df=p.yz_timepoint(path)

timepoints=df['time (s)'].unique()
for t in timepoints:
    plt.figure()
    sub_df=df[df['time (s)']==t]
    sns.lineplot(x='y-z length (um)',y='height (um)',hue='y-z selection',data=sub_df)
    # plt.title('time = ',t)
    