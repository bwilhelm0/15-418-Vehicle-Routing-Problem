import json
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
#import matplotlib.pyplot as plt
import seaborn as sns

infile = "uber_data/uber.csv"
outfile_filtered = "uber_data/uber_filtered.csv"
outfile_condensed = "uber_data/uber_condensed.csv"
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


# trips_condensed['pickup_id'] = trips_condensed.groupby(['pickup_longitude','pickup_latitude']).ngroup()
# trips_condensed['dropoff_id'] = trips_condensed.groupby(['dropoff_longitude','dropoff_latitude']).ngroup()

pickup_longitude = trips_condensed['pickup_longitude'].tolist()
pickup_latitude = trips_condensed['pickup_latitude'].tolist()

dropoff_longitude = trips_condensed['dropoff_longitude'].tolist()
dropoff_latitude = trips_condensed['dropoff_latitude'].tolist()

node_dict = {}

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

# for row in trips_condensed:
#     row['pickup_node'] = node_dict[(row['pickup_longitude'],row['pickup_latitude'])] 
#     row['dropoff_node'] = node_dict[(row['dropoff_longitude'],row['dropoff_latitude'])]

trips_condensed['dropoff_node'] = trips_condensed.apply(lambda row: node_dict.get((row['dropoff_longitude'],row['dropoff_latitude']),-1),axis = 1)
trips_condensed['pickup_node'] = trips_condensed.apply(lambda row: node_dict.get((row['pickup_longitude'],row['pickup_latitude']),-1),axis = 1)


print(node_dict)
# for row in trips_condensed:
#     print(int(row['pickup_latitude']))

trips_condensed = trips_condensed[['dropoff_node','pickup_node','distance','dropoff_longitude','dropoff_latitude','pickup_longitude','pickup_latitude']]
print(trips_condensed)
trips_filtered.to_csv(outfile_filtered)
trips_condensed.to_csv(outfile_condensed)

