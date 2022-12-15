import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

models =df = pd.DataFrame({
    'Model':['SVC','XGBClassifier', 'LGBMClassifier','LogisticRegression','LinearSVC','CalibratedClassifierCV',
             'SGDClassifier','ExtraTreesClassifier',
             'RandomForestClassifier','PassiveAggressiveClassifier'
        ,'Perceptron','BaggingClassifier','NuSVC', 'AdaBoostClassifier',
             'QuadraticDiscriminantAnalysis', 'DecisionTreeClassifier', 'RidgeClassifier', 'RidgeClassifierCV',
             'LinearDiscriminantAnalysis', 'NearestCentroid','ExtraTreeClassifier','BernoulliNB','KNeighborsClassifier',
             'GaussianNB', 'DummyClassifier'],
    'Accuracy':[0.96,0.93,0.91,0.91,0.91,0.91,0.90,0.90,0.90,0.89,0.87,0.86,0.83,0.81,0.80,0.80,0.79,0.79,0.79,0.77,0.76,0.75,0.74,0.70,0.33]
})
#KNeighborsClassifier               0.74               0.74    None      0.74      470.67
#GaussianNB                         0.70               0.70    None      0.70       37.73
#DummyClassifier                    0.33               0.33    None      0.16        9.26

print(models)
#sns.set(font_scale=0.5)
plt.figure(figsize=(5, 10))
sns.set_theme(style="whitegrid")
ax = sns.barplot(y='Model', x='Accuracy', data=models)

plt.show()

plt.figure(figsize=(10, 5))
sns.set_theme(style="whitegrid")
ax = sns.barplot(x='Model', y='Accuracy', data=models)
plt.xticks(rotation=90)

plt.show()

