
import pandas as pd
import numpy as np

df = pd.read_csv("DNA_step_library_frame.csv")

origin0 = np.array(map(float,df.ix[0]['origin'].split(',')))
M_axis0 = np.array(map(float,df.ix[0]['M_axis'].split(','))).reshape(3,3)

origin1 = np.array(map(float,df.ix[1]['origin'].split(',')))
M_axis1 = np.array(map(float,df.ix[1]['M_axis'].split(','))).reshape(3,3)

print origin0
print M_axis0
print origin1
print M_axis1

distance_score = np.sqrt(np.sum((origin0-origin1)**2)) + 2*np.sqrt(np.sum((M_axis0-M_axis1)**2))
print distance_score



