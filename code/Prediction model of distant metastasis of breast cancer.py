import numpy as np
import pandas as pd
import sklearn
import matplotlib as mlp
import matplotlib.pyplot as plt
import seaborn as sns
import time
import re, pip, conda
from sklearn.metrics import *
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.model_selection import cross_validate, KFold, GridSearchCV,train_test_split

data = pd.read_csv(r"D:\model\data.csv",index_col=0)

X = data.iloc[:,:-1]
y = data.iloc[:,-1]
X_train,X_test,y_train,y_test=train_test_split(X, y, test_size=0.3,random_state=15)

param_grid_simple = {
                     "criterion": ["gini","entropy"]
                     ,'n_estimators': [*range(60,100,5)]
                     , 'max_depth': [*range(2,8,1)]
                    ,"max_features": ["log2",16,21,"auto"]
                    }

reg = RFC(random_state=1412,verbose=True,n_jobs=8)
cv = KFold(n_splits=5,shuffle=True,random_state=1412)
search = GridSearchCV(estimator=reg
                     ,param_grid=param_grid_simple
                     ,scoring = "roc_auc"
                     ,verbose = True
                     ,cv = cv
                     ,n_jobs=8)
search.fit(X_train,y_train)

search.best_estimator_