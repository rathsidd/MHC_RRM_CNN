import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import utils

convolved_mhc = utils.load_data('./CNN_with_EIIP.csv', arr_cols=['conv_length'], names=['convolved'])

from lazypredict.Supervised import LazyRegressor
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(pd.DataFrame(convolved_mhc['convolved'].to_list()), convolved_mhc['pIC50'], test_size=0.33, random_state=42)

reg = LazyRegressor(verbose=0,ignore_warnings=False, custom_metric=None )
models,predictions = reg.fit(X_train, X_test, y_train, y_test)