import json
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import seaborn as sns
import random

infile = "TSP/uber_data/uber.csv"
outfile_filtered = "TSP/uber_data/uber_filtered.csv"
outfile_condensed = "TSP/uber_condensed.csv"
trips_raw = pd.read_csv(infile)

#Time information is presented in two columns (key and pickup_datetime).
#data in the first column (Unnamed:0) will not be helpful to our analytics.
#removing the column 'Unnamed:0' and the column named 'key'.
del trips_raw['Unnamed: 0']
del trips_raw['key']
del trips_raw['passenger_count']
del trips_raw['pickup_datetime']


lat1 = trips_raw['pickup_latitude']
lon1 = trips_raw['pickup_longitude']
lat2 = trips_raw['dropoff_latitude']
lon2 = trips_raw['dropoff_longitude']
trips_raw['distance'] = np.sqrt((lat1 - lat2)**2 + (lon1 - lon2)**2)

trips_raw = trips_raw.round({'pickup_longitude':1, 'pickup_latitude': 1, 
                            'dropoff_latitude':1,'dropoff_longitude':1})

#remove outliers in data
fare_low = trips_raw['fare_amount'].quantile(0.05)
fare_hi = trips_raw['fare_amount'].quantile(0.95)
dist_low = trips_raw['distance'].quantile(0.05)
dist_hi = trips_raw['distance'].quantile(0.95)

trip_outliers = trips_raw[(trips_raw["distance"] < dist_hi) & 
    (trips_raw["distance"] > dist_low) & (trips_raw['fare_amount'] < fare_hi)
    & (trips_raw['fare_amount'] > fare_low)]

trips_filtered = trip_outliers
trips_filtered = trips_filtered.groupby(['pickup_longitude','pickup_latitude','dropoff_latitude','dropoff_longitude']).mean()
trips_filtered['edge'] = range(1, len(trips_filtered) + 1)
trips_condensed = trips_filtered[['edge','distance']]
trips_condensed = trips_condensed.reset_index()



#convert columns to lists 
pickup_longitude = trips_condensed['pickup_longitude'].tolist()
pickup_latitude = trips_condensed['pickup_latitude'].tolist()
dropoff_longitude = trips_condensed['dropoff_longitude'].tolist()
dropoff_latitude = trips_condensed['dropoff_latitude'].tolist()

node_dict = {}



#give unique ids to every node, count should be the total amount of unique nodes at the end of the loops
count = 0 
for i in range(len(pickup_longitude)):
    index = node_dict.get((pickup_longitude[i],pickup_latitude[i]),-1)
    if (index == -1):
        node_dict[(pickup_longitude[i],pickup_latitude[i])] = count
        count+=1

for i in range(len(dropoff_longitude)):
    index = node_dict.get((dropoff_longitude[i],dropoff_latitude[i]),-1)
    if (index == -1):
        node_dict[(dropoff_longitude[i],dropoff_latitude[i])] = count
        count+=1



#apply dictionary get function for every row to apply node enumeration
trips_condensed['dropoff_node'] = trips_condensed.apply(lambda row: node_dict.get((row['dropoff_longitude'],row['dropoff_latitude']),-1),axis = 1)
trips_condensed['pickup_node'] = trips_condensed.apply(lambda row: node_dict.get((row['pickup_longitude'],row['pickup_latitude']),-1),axis = 1)

nodes = 64
for drop_node in range(nodes):
    for pickup_node in range(nodes):
        if not (((trips_condensed['dropoff_node']==drop_node) & (trips_condensed['pickup_node']==pickup_node)).any()):
            print(f'not included dropnode: {drop_node} pickupnode: {pickup_node}')
            filler_distance = random.uniform(0.006,0.09)
            new_row = {'dropoff_node': drop_node, 'pickup_node':pickup_node,'distance':filler_distance}

            print(trips_condensed.columns)

            df_new_row = pd.DataFrame({ 'pickup_longitude': [0], 'pickup_latitude': [0],'dropoff_latitude':[0],'dropoff_longitude':[0],'edge':[0],'distance':[filler_distance],'dropoff_node':[drop_node],'pickup_node':[pickup_node]})
            trips_condensed = pd.concat([trips_condensed, df_new_row])               
           








print(node_dict)

trips_condensed = trips_condensed[['dropoff_node','pickup_node','distance','dropoff_longitude','dropoff_latitude','pickup_longitude','pickup_latitude']]
print(trips_condensed)
trips_filtered.to_csv(outfile_filtered)
trips_condensed.to_csv(outfile_condensed)

