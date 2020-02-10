PayPal Assignment
Submitted by Nomi Hadar
September 2019
#################################################

Things I would improve with more time:

1. Handling the missing values.
In this dataset, 3 out of 11 features had missing values.
There are several approaches to handle missing values.
The simple one is just to remove rows with missing data.
However, in this case we would lose ~50% of the data.
Instead, in 2 of the 3 features I filled the missing values with the mean value.
This can be improved by filling it with the mean plus the standard deviation.
(Median can be better if there are outliers, but from quick glance data doesn't contain outliers).
However, in the feature of "age", it doesn't make sense to fill it with mean,
so I filled the missing values with zero, and add a bool feature "is_null",
indicating if there is a missing value. 
The disadventage of this is that we add a feature. 
To sum up, it's important to know the disadventages of each approach,
and try to "play" with it more.

2. Features
Several features were pre-prccessed, such that they will be useful. 
For example - since I think the most important component in timestamp is the hour, I used hour only.
And, because time is cyclic, and 23 should be in an 1 distance from 00, 
a sinus cosinus manipulation was needed.
However, if timestamp is not from a single time zone, then this feature is probably useless.
Also, I splited the category into main category (like clothing )and subcategory. 
Not sure in this split, should be examined. 


3. Model
I used the Gradient Boosting model.
If I would have more time, I would think more thoroughly 
which model is most suitable for data.
Also, I would "play" with the GB model parameters more.

4. Feature importance
I would do fetautres importance to know which features contribute the most,
ang get rid of those which not at all.

5. Visualization
I would do more data visualization, to get to know data better.



Thanks!
Nomi