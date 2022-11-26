import json
import numpy as np # linear algebra
import pandas as pd # data processing, CSV file I/O (e.g. pd.read_csv)
import matplotlib.pyplot as plt
import seaborn as sns

infile = "uber_data/uber.csv"
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



trips_filtered = trips_filtered.groupby(['pickup_longitude','pickup_latitude','dropoff_latitude','dropoff_longitude'])


#['distance'].apply(', '.join).reset_index()
print(trips_filtered.size())
print(trips_filtered.head(10))
