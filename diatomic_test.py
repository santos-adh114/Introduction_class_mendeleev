import numpy as np
import pandas as pd

def string_number_separator(string):
    import re
    temp = re.compile("([a-zA-Z]+)([0-9]+)")
    res = temp.match(string).groups()
    return res

def list_of_atoms_and_numbers(string):
    atoms=[]
    atoms_number=[]
    str1=''
    while len(str1)<len(string):
        string1=str(string)[len(str1):]
        x,y=string_number_separator(string1)
        str1=str1+str(x)+str(y)
        atoms.append(str(x))
        atoms_number.append(int(y))
    return atoms, atoms_number

df=pd.read_csv('pbe_krr_lap.csv')
def return_dataframe(df):
    problematic_list=[]
    for entries in df.compound.unique().tolist():
        a,b=list_of_atoms_and_numbers(entries)
        if len(a)>=3:
            problematic_list.append(entries)
        else:
            pass
    df_new=df[df['compound'].isin(problematic_list)]
    return df_new
#df_new.to_csv('ternary_and_binary_summary.csv')
#df_features=pd.read_csv('pbe_features_ef_ed.csv')
#df_ter_quater=return_dataframe(df_features)
#df_ter_quater=df_ter_quater[['compound','pbe_hf','expt_hf','pbe_hf_err']]
#print(df_ter_quater.head())

def error_analysis(df):
    error_oxide=[]
    error_non_oxide=[]
    for entries in df.compound.unique().tolist():
        a,b=list_of_atoms_and_numbers(entries)
        value=1000*df[df.compound==entries].pbe_hf_err.tolist()[0]
        if 'O' in a:
            error_oxide.append(value)
        else:
            error_non_oxide.append(value)
    oxides=np.array(error_oxide);non_oxides=np.array(error_non_oxide)
    return np.sort(oxides),np.sort(non_oxides)
#oxides,non_oxides=error_analysis(df_ter_quater)
#print(len(a),len(b),a,b)
def statistics(lst):
    lst=np.sort(lst)
    max_err=lst.max()
    min_err=lst.min()
    diff=abs(max_err-min_err)
    abs_avg=abs(lst).mean()
    avg=lst.mean()
    med=np.percentile(lst,50)
    std=lst.std()
    q1=np.percentile(lst,25)
    q3=np.percentile(lst,75)
    return len(lst),abs_avg,avg,std,min_err,q1,med,q3,max_err
#print(oxides.sort())
#print(statistics(oxides))
#print(statistics(non_oxides))

#import plotly.express as px

def return_updated_dataframe(df,xc):
    oxide=[]
    class_type=[]
    for entries in df.compound.unique().tolist():
        a,b=list_of_atoms_and_numbers(entries)
        if len(a)==3:
            class_type.append('Ternary')
        elif len(a)==4:
            class_type.append('Quaternary')
        else:
            class_type.append('Binary')
        if ('O' in a) or ('F' in a) or ('N' in a) or ('Cl' in a) or ('H' in a):
            oxide.append('Diatomic')
        else:
            oxide.append('No Diatomic')
    df_new=df
    df_new[[str(xc)+'_hf_err',str(xc)+'_hf','expt_hf']]=1000*df_new[[str(xc)+'_hf_err',str(xc)+'_hf','expt_hf']]
    df_new['diatomic']=oxide
    df_new['class']=class_type
    return df_new

#df_new.to_csv('ternary_and_binary_summary.csv')
def violin_plot(df_features,xc,x,y,hue,plotname): # x='class',y='pbe_hf_err',hue='oxide'
    import seaborn as sns
    import matplotlib.pyplot as plt
    fontsize=12
    df_features[str(xc)+'_hf_err']=abs(df_features[str(xc)+'_hf_err'])
    df_get=return_updated_dataframe(df_features,xc)
    my_col = {"No Diatomic": "mediumslateblue", "Diatomic": "red"}
    sns.violinplot(data=df_get, x=x, y=y, hue=hue,palette=my_col)
    #sns.set_color_codes("bright")
    ylabel='|$\Delta$H$_{f,'+str(xc).upper()+' }$ error| (meV/at)'
    #ylabel='Ionicity'
    xlabel='Class'
    plt.ylabel(ylabel, fontsize=fontsize)
    plt.xlabel(xlabel, fontsize=fontsize)
    plt.legend(fontsize=fontsize-2,frameon=False)
    axis_width=1.5
    plt.subplots_adjust(left=0.2,bottom=0.2)
    plt.tick_params('both', length = 6, width = axis_width, which = 'major',right=True,top=True)
    plt.yticks(fontsize=fontsize)
    plt.xticks(fontsize=fontsize)
    plt.savefig(str(xc)+'_'+plotname+'.pdf',dpi=80) # plotname='binary_vs_ternary_violin'
    plt.show()
    plt.close()

xc='pbe'
df_features=pd.read_csv(str(xc)+'_features_ef_ed.csv')
#violin_plot(df_features,xc,'class','ionicity','oxide','ionicity_binary_vs_ternary_violin')
#violin_plot(df_features,xc,'class',str(xc)+'_hf_err','oxide','hf_err_binary_vs_ternary_violin')
violin_plot(df_features,xc,'class',str(xc)+'_hf_err','diatomic','abs_hf_err_binary_vs_ternary_violin')
'''
df_updated=return_updated_dataframe(df_features,xc)
final_list=[]
for entries in ['Binary','Ternary','Quaternary']:
    for types in ['Oxygen','No Oxygen']:
        #print(entries,types)
        list_interest=np.array(df_updated[(df_updated['class']==entries) & (df_updated['oxide']==types)][str(xc)+'_hf_err'].tolist())
        #print(list_interest)
        #list_interest=np.array(df_updated[(df_updated['class']==entries) & (df_updated['oxide']==types)]['ionicity'].tolist())
        count,abs_avg,avg,std,min_err,q1,med,q3,max_err=statistics(list_interest)
        #print('%d\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f\t %6.2f'%(count,abs_avg,avg,std,min_err,q1,med,q3,max_err))
        final_list.append([entries,types,count,abs_avg,avg,std,min_err,q1,med,q3,max_err])
#df_tot=pd.DataFrame(final_list)
#df_tot.to_csv(str(xc)+'summary.csv')
#print()
'''
pass
