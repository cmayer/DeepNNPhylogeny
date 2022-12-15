import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

models =df = pd.DataFrame({
    'Model':['NuSVC','SVC','CalibratedClassifierCV','BaggingClassifier','LinearSVC','RandomForestClassifier',
             'AdaBoostClassifier','XGBClassifier','LinearDiscriminantAnalysis','KNeighborsClassifier','ExtraTreesClassifier','PassiveAggressiveClassifier',
            'RidgeClassifierCV','DecisionTreeClassifier', 'LGBMClassifier', 'SGDClassifier', 'Perceptron', 'QuadraticDiscriminantAnalysis',
             'LogisticRegression', 'NearestCentroid', 'GaussianNB', 'RidgeClassifier', 'ExtraTreeClassifier', 'BernoulliNB', 'DummyClassifier'
             ],
    'Time Taken':[158537.34,24110.06,4252.85,2349.07, 1413.49,1029.35,1009.68,955.87,724.17,470.67,465.56,385.54,351.46,345.67,176.52,175.87,
                  126.85,87.17,50.28,39.26,37.73,19.01,14.25,12.93,9.26]
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

