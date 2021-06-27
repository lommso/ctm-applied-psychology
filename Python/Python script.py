#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import seaborn as sns
get_ipython().run_line_magic('matplotlib', 'inline')
from scipy import stats


# ### 1. Extract relevant columns

# List of manifest and latent variables:

# In[2]:


manifests = {'pa18i1':  {'max': 5,  'min': 1, 'inverted': True,  'desc': 'Partner finds it all right if I stand up for my own interests'}, 
             'pa18i2':  {'max': 5,  'min': 1, 'inverted': False, 'desc': 'Sometimes I am afraid that partner would rather spend time with others'}, 
             'pa18i4':  {'max': 5,  'min': 1, 'inverted': False, 'desc': 'Partner clings to me so much that I feel like I am suffocating'}, 
             'pa18i6':  {'max': 5,  'min': 1, 'inverted': True,  'desc': 'I can settle my personal matters by myself without conflicts'}, 
             'pa18i7':  {'max': 5,  'min': 1, 'inverted': False, 'desc': 'I have the feeling that I like partner more than he/she likes me'}, 
             'pa18i10': {'max': 5,  'min': 1, 'inverted': False, 'desc': 'Sometimes not sure if partner enjoys being with me as much as I'}, 
             'pa18i11': {'max': 5,  'min': 1, 'inverted': True,  'desc': 'I can usually do what I want'}, 
             'pa18i12': {'max': 5,  'min': 1, 'inverted': False, 'desc': 'Afraid partner will think I am silly/stupid if I make a mistake'}, 
             'pa18i14': {'max': 5,  'min': 1, 'inverted': False, 'desc': 'Partner clings to me so tightly that I cannot do what I want '}, 
             'pa18i15': {'max': 5,  'min': 1, 'inverted': False, 'desc': 'When I disappoint/annoy partner, I am afraid he/she will not like me'}, 
             'pa18i16': {'max': 5,  'min': 1, 'inverted': True,  'desc': 'I can follow own interests without partner getting upset'}, 
             'per1i2':  {'max': 5,  'min': 1, 'inverted': False, 'desc': 'Sometimes I believe that I am worthless'}, 
             'per1i7':  {'max': 5,  'min': 1, 'inverted': True,  'desc': 'I like myself just the way I am'}, 
             'per1i13': {'max': 5,  'min': 1, 'inverted': True,  'desc': 'All in all, I am pleased with myself'},                      
             
             'per1i6':  {'max': 5,  'min': 1, 'inverted': False, 'desc': 'I feel lonely'}, 
             'sat6':    {'max': 10, 'min': 0, 'inverted': False, 'desc': 'General satisfaction with life'}}

latents   = {'attAvd':  {'manifests': ['pa18i4', 'pa18i14', 'pa18i1', 'pa18i6', 'pa18i11', 'pa18i16'], 'desc': 'Attachment Anxiety'},
             'attAnx':  {'manifests': ['pa18i7', 'pa18i10', 'pa18i2', 'pa18i12', 'pa18i15', 'per1i2', 'per1i7', 'per1i13'], 'desc': 'Attachment Avoidance'}}


# Mapping, which variables were observed at the different measurement occasions:

# In[3]:


aval_cols = {1:  ['sat6', 'per1i6'] + latents['attAvd']['manifests'] + latents['attAnx']['manifests'],
             2:  ['sat6'],
             3:  ['sat6']           + latents['attAvd']['manifests'] + latents['attAnx']['manifests'],
             4:  ['sat6', 'per1i6'],
             5:  ['sat6', 'per1i6'] + latents['attAvd']['manifests'] + latents['attAnx']['manifests'],
             6:  ['sat6'],
             7:  ['sat6', 'per1i6'] + latents['attAvd']['manifests'] + latents['attAnx']['manifests'],
             8:  ['sat6', 'per1i6'],
             9:  ['sat6', 'per1i6'] + latents['attAvd']['manifests'] + latents['attAnx']['manifests'],
             10: ['sat6', 'per1i6'],
             11: ['sat6', 'per1i6'] + latents['attAvd']['manifests'] + latents['attAnx']['manifests']}


# Load the data for the relevant columns:

# In[4]:


data_anchor = {}
for wave in [1,2,3,4,5,6,7,8,9,10,11]:
    filepath = '../../data/pairfam_v11/Data/Stata/anchor'+str(wave)+'.dta'
    columns = ['sample', 'wave', 'id', 'pid', 'sex_gen', 'original_doby'] + aval_cols[wave]
    data_anchor[wave] = pd.read_stata(filepath, columns=columns, convert_categoricals=False)
    
data_partner = {}
for wave in [1,2,3,4,5,6,7,8,9,10,11]:
    filepath = '../../data/pairfam_v11/Data/Stata/partner'+str(wave)+'.dta'
    columns = ['sample', 'wave', 'id', 'pid', 'psex', 'pdoby'] + ['p'+col for col in aval_cols[wave]]
    data_partner[wave] = pd.read_stata(filepath, columns=columns, convert_categoricals=False)


# ### 2. Merge data from different waves

# In[5]:


data_anchor  = pd.concat(data_anchor.values(), ignore_index=True)
data_partner = pd.concat(data_partner.values(), ignore_index=True)


# ### 3. Keep only the main pairfam sample

# In[6]:


data_anchor  = data_anchor[data_anchor['sample']==1]
data_anchor.drop('sample', axis=1, inplace=True)
data_partner = data_partner[data_partner['sample']==1]
data_partner.drop('sample', axis=1, inplace=True)


# In[7]:


data_partner[(data_partner['pper1i6']>0) & (data_partner['wave']==11)]


# In[8]:


data_anchor[data_anchor.id==734104000]


# Analyze sample size:

# In[9]:


summary = pd.DataFrame(columns=['w1', 'w2', 'w3', 'w4', 'w5', 'w6', 'w7', 'w8', 'w9', 'w10', 'w11'])

for wave in range(1,12):
    N_all = sum(data_anchor.wave==wave)
    N_inRelationship = sum((data_anchor.wave==wave) & (data_anchor.pid > 0))
    N_partnerData = sum(data_anchor[data_anchor.wave==wave].pid.isin(data_partner[data_partner.wave==wave].pid))
    summary['w'+str(wave)] = [N_all, N_inRelationship, N_partnerData]

summary.index = ['Total number of subjects', 'among those: in a relationship', 'among those: partners participated']
summary


# ### 4. Exclude singles

# In[10]:


data_anchor.dropna(subset=['pid'], inplace=True)


# ### 5. Mask missing values

# In[11]:


data_anchor[data_anchor<0] = np.nan
data_partner[data_partner<0] = np.nan


# Analyze descriptives for anchors:

# In[12]:


desc_anchor = pd.DataFrame()

for col in manifests:
    N = data_anchor[[col, 'wave']].groupby('wave').count().T.add_prefix('N_')
    m = np.round(data_anchor[col].mean(),3)
    sd = np.round(data_anchor[col].std(),3)
    desc_anchor = desc_anchor.append(pd.DataFrame({'desc': manifests[col]['desc'], 
                                                   'min': manifests[col]['min'], 'max': manifests[col]['max'], 
                                                   'mean': m, 'sd': sd}, index=[col]).join(N))   
desc_anchor


# Analyze descriptives for partners:

# In[13]:


desc_partner = pd.DataFrame()

for manifest in manifests:
    col='p'+manifest
    N = data_partner[[col, 'wave']].groupby('wave').count().T.add_prefix('w')
    m = np.round(data_partner[col].mean(),3)
    sd = np.round(data_partner[col].std(),3)
    desc_partner = desc_partner.append(pd.DataFrame({'desc': manifests[manifest]['desc'], 
                                                     'min': manifests[manifest]['min'], 'max': manifests[manifest]['max'], 
                                                     'mean': m, 'sd': sd}, index=[col]).join(N))   
    
desc_partner


# ### 6. Invert variable scales

# In[14]:


for manifest in manifests: 
    if(manifests[manifest]['inverted']): 
        data_anchor[manifest] = manifests[manifest]['max'] - data_anchor[manifest] + manifests[manifest]['min']
        data_partner['p'+manifest] = manifests[manifest]['max'] - data_partner['p'+manifest] + manifests[manifest]['min']


# Compute Cronbach alpha's for the different scales:

# In[15]:


scales = {'Engulfment anxiety     ': ['pa18i4', 'pa18i14'],
          'Autonomy               ': ['pa18i1', 'pa18i6', 'pa18i11', 'pa18i16'],
          'Ambivalence            ': ['pa18i7', 'pa18i10'],
          'Fear of love withdrawal': ['pa18i2', 'pa18i12', 'pa18i15'],
          'Self-esteem            ': ['per1i2', 'per1i7', 'per1i13'],
          
          'Attachment Anxiety     ': ['pa18i4', 'pa18i14', 'pa18i1', 'pa18i6', 'pa18i11', 'pa18i16'],
          'Attachment Avoidance   ': ['pa18i7', 'pa18i10', 'pa18i2', 'pa18i12', 'pa18i15', 'per1i2', 'per1i7', 'per1i13']}


# In[16]:


for scale in scales:
    N = len(scales[scale])
    mean_r = np.mean((data_anchor[scales[scale]].corr().sum()-1)/(N-1))
    cronbach_alpha = (N * mean_r) / (1 + (N - 1) * mean_r)
    print(scale + ':\t' + str(round(cronbach_alpha,3)))


# ### 7. Increase sample size:

# Check t-statistic:

# In[27]:


summary = pd.DataFrame(columns = ['ID', 'Description', 't-value', 'p-value', 'mean (anchor)', 'mean (partner)', 'sd (anchor)', 'sd (partner)'])

for manifest in manifests:
    val_a = data_anchor[manifest]
    val_p = data_partner['p'+manifest]
    t     = np.round(stats.ttest_ind(val_a, val_p, nan_policy='omit')[0],7)
    p     = np.round(stats.ttest_ind(val_a, val_p, nan_policy='omit')[1],7)
    m_a   = np.round(val_a.mean(),3)
    m_p   = np.round(val_p.mean(),3)
    sd_a  = np.round(val_a.std(),3)
    sd_p  = np.round(val_p.std(),3)
    summary.loc[len(summary)] = [manifest, manifests[manifest]['desc'], t, p, m_a, m_p, sd_a, sd_p]

pd.set_option('display.max_colwidth', None)
summary


# In[21]:


manifests[manifest]


# Make labels in both datasets match each other:

# In[47]:


data_partner = data_partner.rename(columns={'pid': 'pid', 'id': 'ppid'}) # Rename IDs so that next row works for them too
data_partner.columns = data_partner.columns.str.replace('^p', '') # Remove prefix 'p'
data_anchor = data_anchor.rename(columns={'sex_gen': 'sex', 'original_doby': 'doby'}) # Rename column labels so they match


# Standardize variables:

# In[48]:


data_anchor[list(manifests.keys())] = (data_anchor[manifests]-data_anchor[manifests].mean())/data_anchor[manifests].std()
data_partner[list(manifests.keys())] = (data_partner[manifests]-data_partner[manifests].mean())/data_partner[manifests].std()


# Append the two dataframes:

# In[49]:


data_anchor['subj_type'] = 'anchor'
data_partner['subj_type'] = 'partner'
data = pd.concat([data_anchor, data_partner], ignore_index=True)


# ### 8. Calculate latent variables:

# In[50]:


for lat in latents:
    data[lat] = data[latents[lat]['manifests']].mean(axis=1, skipna=False)
    #data.drop(latents[lat]['manifests'], axis=1, inplace=True)


# ### 9. Calculate age column

# In[51]:


data[data.id==907000].T


# In[52]:


data['age'] = 2009-data['doby']+data['wave']
data.drop('doby', axis=1, inplace=True)


# ### 10. Combine anchor and partner data:

# In[53]:


data_join = data.copy()
for col in data_join.columns.difference(['wave','id','pid']):
    data_join.rename(columns={col:'p'+col}, inplace=True)
    
data = data.merge(data_join, left_on=['wave','id'], right_on=['wave','pid'], how='left', suffixes=['', '_p']) 
data.drop(['id_p', 'pid_p'], axis=1, inplace=True)


# ### Export data

# In[55]:


data.to_csv('../../data/samples/data8.csv', index=False)


# In[ ]:




