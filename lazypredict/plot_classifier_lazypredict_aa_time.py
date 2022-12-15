import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

models =df = pd.DataFrame({
    'Model':['NuSVC','SVC','XGBClassifier', 'LGBMClassifier', 'AdaBoostClassifier', 'LinearDiscriminantAnalysis',
             'CalibratedClassifierCV', 'BaggingClassifier', 'KNeighborsClassifier', 'QuadraticDiscriminantAnalysis', 'LinearSVC', 'DecisionTreeClassifier',
             'RidgeClassifierCV', 'PassiveAggressiveClassifier', 'Perceptron', 'RidgeClassifier', 'SGDClassifier', 'LabelSpreading',
             'LabelPropagation', 'ExtraTreesClassifier', 'RandomForestClassifier', 'GaussianNB', 'BernoulliNB', 'NearestCentroid', 'ExtraTreeClassifier', 'DummyClassifier'
             ],
    'Time Taken':[9977.59,7198.56,7130.53,3828.11,1290.92,1262.56,1232.80,1058.16, 641.09,600.81,345.87,246.68,181.45,161.53,140.64,139.79,138.73,137.71,137.37,
                  100.00,92.28,76.63,71.14,60.63,58.52,54.95]
})


print(models)
#sns.set(font_scale=0.5)
plt.figure(figsize=(5, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y='Model', x='Time Taken', data=models)

plt.show()

plt.figure(figsize=(10, 5))
sns.set_theme(style="whitegrid")
ax = sns.barplot(x='Model', y='Time Taken', data=models)
plt.xticks(rotation=90)

plt.show()

