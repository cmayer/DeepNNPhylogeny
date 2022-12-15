import matplotlib.pyplot as plt
import seaborn as sns

import pandas as pd

models =df = pd.DataFrame({
    'Model':['LGBMClassifier','XGBClassifier','AdaBoostClassifier','BaggingClassifier','RandomForestClassifier','ExtraTreesClassifier','PassiveAggressiveClassifier',
             'LinearSVC','CalibratedClassifierCV','SGDClassifier','Perceptron','DecisionTreeClassifier', 'BernoulliNB', 'SVC', 'NearestCentroid',
             'NuSVC','GaussianNB', 'RidgeClassifier','RidgeClassifierCV', 'LinearDiscriminantAnalysis', 'ExtraTreeClassifier',
             'KNeighborsClassifier', 'LabelPropagation', 'LabelSpreading', 'DummyClassifier', 'QuadraticDiscriminantAnalysis'],
    'Accuracy':[0.98,0.97,0.95,0.95,0.95,0.94,0.93,0.93,0.93,0.92,0.90,0.90,0.84,0.83,0.82,0.81,0.78,0.78,0.78,0.78,0.75,0.66,0.41,0.41,0.35,0.35]
})



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

