# Number of succ. JTT
JTT_b1 = 944 + 923 + 931 
JTT_b1 # 2798
JTT_b2 = 937 + 950 + 948 
JTT_b2 # 2835 
JTT_b3 =  968 + 941 + 956 
JTT_b3  # 2865 
JTT_b10 = 870 + 885 + 829 
JTT_b10 # 2584 
JTT_u = 931 + 934 + 944 
JTT_u  # 2809 
JTT_cu = 920 + 933 + 936 
JTT_cu  # 2789
JTT_cud = 952 + 925 + 936 
JTT_cud   # 2813 
JTT_NJ = 966 + 907 + 904 
JTT_NJ # 2777
JTT_ML = 1000 + 999 + 999
JTT_ML # 2998 

JTT_b1/30 # 93.27%
JTT_b2/30 # 94.5%
JTT_b3/30 # 95.5%
JTT_b10/30 # 86.13%
JTT_u/30 # 93.63%
JTT_cu/30 # 92.97%
JTT_cud/30 # 93.77%
JTT_NJ/30 # 92.57%
JTT_ML/30 # 99.93%

# Binomal test JTT  

prop.test(c(JTT_b1,JTT_b2),c(3000,3000)) # almost diff. p-value = 0.05245
prop.test(c(JTT_b1,JTT_b3),c(3000,3000)) # sign diff p-value = 0.000215
prop.test(c(JTT_b1,JTT_b10),c(3000,3000)) # sign diff  p-value < 2.2e-16
prop.test(c(JTT_b1,JTT_u),c(3000,3000)) # no diff p-value = 0.6018
prop.test(c(JTT_b1,JTT_cu),c(3000,3000)) # no diff p-value = 0.6833
prop.test(c(JTT_b1,JTT_cud),c(3000,3000)) # no diff p-value = 0.4629
prop.test(c(JTT_b1,JTT_NJ),c(3000,3000)) # no diff p-value = 0.3142
prop.test(c(JTT_b1,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_b2,JTT_b3),c(3000,3000)) # almost diff p-value = 0.08583
prop.test(c(JTT_b2,JTT_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_b2,JTT_u),c(3000,3000)) # no diff p-value = 0.1719
prop.test(c(JTT_b2,JTT_cu),c(3000,3000)) # sign diff p-value = 0.01653
prop.test(c(JTT_b2,JTT_cud),c(3000,3000)) # no diff p-value = 0.2486
prop.test(c(JTT_b2,JTT_NJ),c(3000,3000)) # sign diff p-value = 0.002771
prop.test(c(JTT_b2,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_b3,JTT_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_b3,JTT_u),c(3000,3000)) # sign diff p-value = 0.001734
prop.test(c(JTT_b3,JTT_cu),c(3000,3000)) # sign diff p-value = 3.273e-05
prop.test(c(JTT_b3,JTT_cud),c(3000,3000)) # sign diff p-value = 0.003482
prop.test(c(JTT_b3,JTT_NJ),c(3000,3000)) # sign diff p-value = 2.119e-06
prop.test(c(JTT_b3,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_b10,JTT_u),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_b10,JTT_cu),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_b10,JTT_cud),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_b10,JTT_NJ),c(3000,3000)) # sign diff p-value = 9.333e-16
prop.test(c(JTT_b10,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_u,JTT_cu),c(3000,3000)) # no diff p-value = 0.3266
prop.test(c(JTT_u,JTT_cud),c(3000,3000)) # no diff p-value = 0.8733
prop.test(c(JTT_u,JTT_NJ),c(3000,3000)) # no diff  p-value = 0.1143
prop.test(c(JTT_u,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_cu,JTT_cud),c(3000,3000)) # no diff p-value = 0.2328
prop.test(c(JTT_cu,JTT_NJ),c(3000,3000)) # no diff p-value = 0.5835
prop.test(c(JTT_cu,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_cud,JTT_NJ),c(3000,3000)) # almost diff p-value = 0.07333
prop.test(c(JTT_cud,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_NJ,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# Number of suc. LG 
LG_b1 = 923 + 944 + 943 
LG_b1 # 2810
LG_b2 = 921 + 947 + 939 
LG_b2 # 2807
LG_b3 = 939 + 959 + 947 
LG_b3 # 2845
LG_b10 = 881 + 862 + 881 
LG_b10 # 2626
LG_u = 916 + 933 + 956 
LG_u # 2805
LG_cu = 931 + 930 + 952 
LG_cu # 2813
LG_cud = 934 + 927 + 937 
LG_cud #   2798
LG_NJ = 970 + 910 + 910
LG_NJ # 2790
LG_ML = 1000 + 999 + 998 
LG_ML # 2997

LG_b1/30 # 93.67%
LG_b2/30 # 93.57%
LG_b3/30 # 94.83%
LG_b10/30 # 87.47%
LG_u/30 # 93.5%
LG_cu/30 # 93.77%
LG_cud/30 # 93.27%
LG_NJ/30 # 93.0%
LG_ML/30 # 99.9%

# Binomal test LG 

prop.test(c(LG_b1,LG_b2),c(3000,3000)) # no diff  p-value = 0.9159
prop.test(c(LG_b1,LG_b3),c(3000,3000)) # almost diff p-value = 0.05936
prop.test(c(LG_b1,LG_b10),c(3000,3000)) # sign diff p-value = 3.056e-16
prop.test(c(LG_b1,LG_u),c(3000,3000)) # no diff p-value = 0.8331
prop.test(c(LG_b1,LG_cu),c(3000,3000)) # no diff p-value = 0.9153
prop.test(c(LG_b1,LG_cud),c(3000,3000)) # no diff p-value = 0.5655
prop.test(c(LG_b1,LG_NJ),c(3000,3000)) # no diff p-value = 0.3254
prop.test(c(LG_b1,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_b2,LG_b3),c(3000,3000)) # sign diff p-value = 0.041
prop.test(c(LG_b2,LG_b10),c(3000,3000)) # sign diff p-value = 1.061e-15
prop.test(c(LG_b2,LG_u),c(3000,3000)) # no diff p-value = 0.9581
prop.test(c(LG_b2,LG_cu),c(3000,3000)) # no diff p-value = 0.791
prop.test(c(LG_b2,LG_cud),c(3000,3000)) # no diff p-value = 0.6771
prop.test(c(LG_b2,LG_NJ),c(3000,3000)) # no diff p-value = 0.4093
prop.test(c(LG_b2,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_b3,LG_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_b3,LG_u),c(3000,3000)) # sign diff p-value = 0.0317
prop.test(c(LG_b3,LG_cu),c(3000,3000)) # almost diff p-value = 0.08431
prop.test(c(LG_b3,LG_cud),c(3000,3000)) # sign diff p-value = 0.01206
prop.test(c(LG_b3,LG_NJ),c(3000,3000)) # sign diff p-value = 0.003539
prop.test(c(LG_b3,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_b10,LG_u),c(3000,3000)) # sign diff p-value = 2.394e-15
prop.test(c(LG_b10,LG_cu),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_b10,LG_cud),c(3000,3000)) # sign diff p-value = 3.741e-14
prop.test(c(LG_b10,LG_NJ),c(3000,3000)) # sign diff p-value = 7.205e-13
prop.test(c(LG_b10,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_u,LG_cu),c(3000,3000)) # no diff p-value = 0.7113
prop.test(c(LG_u,LG_cud),c(3000,3000)) # no diff p-value = 0.7553
prop.test(c(LG_u,LG_NJ),c(3000,3000)) # no diff p-value = 0.4713
prop.test(c(LG_u,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_cu,LG_cud),c(3000,3000)) # no diff p-value = 0.4629
prop.test(c(LG_cu,LG_NJ),c(3000,3000)) # no diff p-value = 0.2532
prop.test(c(LG_cu,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_cud,LG_NJ),c(3000,3000)) # no diff p-value = 0.7208
prop.test(c(LG_cud,LG_ML),c(3000,3000)) # sign diff  p-value < 2.2e-16

prop.test(c(LG_NJ,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16


# Number of suc. WAG 
WAG_b1 = 925 + 933 + 936 
WAG_b1 # 2794
WAG_b2 = 940 + 926 + 953 
WAG_b2 # 2819
WAG_b3 = 955 + 918 + 952 
WAG_b3 # 2825
WAG_b10 = 856 + 850 + 863 
WAG_b10 # 2569
WAG_u = 954 + 888 + 930 
WAG_u # 2772
WAG_cu = 943 + 908 + 930 
WAG_cu # 2781
WAG_cud = 938 + 904 + 937 
WAG_cud #   2779
WAG_NJ = 964 + 892 + 891 
WAG_NJ # 2747
WAG_ML = 1000 + 999 + 1000
WAG_ML # 2999

WAG_b1/30 # 93.13%
WAG_b2/30 # 93.97%
WAG_b3/30 # 94.17%
WAG_b10/30 # 85.63%
WAG_u/30 # 92.4%
WAG_cu/30 # 92.7%
WAG_cud/30 # 92.63%
WAG_NJ/30 # 91.57%
WAG_ML/30 # 99.97%

# Binomial test WAG

prop.test(c(WAG_b1,WAG_b2),c(3000,3000)) # no diff p-value = 0.2072
prop.test(c(WAG_b1,WAG_b3),c(3000,3000)) # no diff p-value = 0.1122
prop.test(c(WAG_b1,WAG_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_b1,WAG_u),c(3000,3000)) # no diff p-value = 0.2953
prop.test(c(WAG_b1,WAG_cu),c(3000,3000)) # no diff p-value = 0.5459
prop.test(c(WAG_b1,WAG_cud),c(3000,3000)) # no diff p-value = 0.4821
prop.test(c(WAG_b1,WAG_NJ),c(3000,3000)) # sign diff p-value = 0.02547
prop.test(c(WAG_b1,WAG_ML),c(3000,3000)) #sign diff p-value < 2.2e-16

prop.test(c(WAG_b2,WAG_b3),c(3000,3000)) # no diff p-value = 0.7847
prop.test(c(WAG_b2,WAG_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_b2,WAG_u),c(3000,3000)) # sign diff p-value = 0.01846
prop.test(c(WAG_b2,WAG_cu),c(3000,3000)) # almost sign diff p-value = 0.0555
prop.test(c(WAG_b2,WAG_cud),c(3000,3000)) # sign diff p-value = 0.04403
prop.test(c(WAG_b2,WAG_NJ),c(3000,3000)) # sign diff p-value = 0.0004024
prop.test(c(WAG_b2,WAG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_b3,WAG_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_b3,WAG_u),c(3000,3000)) # sign diff p-value = 0.00732
prop.test(c(WAG_b3,WAG_cu),c(3000,3000)) # sign diff p-value = 0.02502
prop.test(c(WAG_b3,WAG_cud),c(3000,3000)) # sign diff  p-value = 0.01929
prop.test(c(WAG_b3,WAG_NJ),c(3000,3000)) # sign diff p-value = 0.0001124
prop.test(c(WAG_b3,WAG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_b10,WAG_u),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_b10,WAG_cu),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_b10,WAG_cud),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_b10,WAG_NJ),c(3000,3000)) # sign diff p-value = 6.48e-13
prop.test(c(WAG_b10,WAG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_u,WAG_cu),c(3000,3000)) # no diff p-value = 0.6941
prop.test(c(WAG_u,WAG_cud),c(3000,3000)) # no diff p-value = 0.7685
prop.test(c(WAG_u,WAG_NJ),c(3000,3000)) #no diff p-value = 0.2539
prop.test(c(WAG_u,WAG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_cud,WAG_cud),c(3000,3000)) # no diff p-value = 0.9605
prop.test(c(WAG_cu,WAG_NJ),c(3000,3000)) # no diff p-value = 0.1135
prop.test(c(WAG_cu,WAG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_cud,WAG_NJ),c(3000,3000)) # no diff p-value = 0.1379
prop.test(c(WAG_cud,WAG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_NJ,WAG_ML),c(3000,3000)) # sign diff  p-value < 2.2e-16

# Number of suc. DAY

DAY_b1 = 932 + 933 + 943 
DAY_b1 # 2808
DAY_b2 = 932 + 936 + 952 
DAY_b2 # 2820
DAY_b3 = 936 + 946 + 939 
DAY_b3 # 2821
DAY_b10 = 849 + 824 + 902
DAY_b10 # 2575
DAY_u = 938 + 924 + 934 
DAY_u # 2796
DAY_cu = 923 + 923 + 930 
DAY_cu # 2776
DAY_cud = 938 + 911 + 929 
DAY_cud # 2778
DAY_NJ = 961 + 900 + 914
DAY_NJ # 2775
DAY_ML = 999 + 999 +999
DAY_ML #   2997

DAY_b1/30 # 93.6%
DAY_b2/30 # 94.0%
DAY_b3/30 # 94.03%
DAY_b10/30 # 85.83%
DAY_u/30 # 93.2%
DAY_cu/30 # 92.53%
DAY_cud/30 # 92.6%
DAY_NJ/30 # 92.5%
DAY_ML/30 # 99.9%

# Binomial test DAY

prop.test(c(DAY_b1,DAY_b2),c(3000,3000)) # no diff p-value = 0.5559
prop.test(c(DAY_b1,DAY_b3),c(3000,3000)) # no diff p-value = 0.5201
prop.test(c(DAY_b1,DAY_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_b1,DAY_u),c(3000,3000)) # no diff p-value = 0.5673
prop.test(c(DAY_b1,DAY_cu),c(3000,3000)) # no diff p-value = 0.1151
prop.test(c(DAY_b1,DAY_cud),c(3000,3000)) # no diff p-value = 0.1396
prop.test(c(DAY_b1,DAY_NJ),c(3000,3000)) # no diff p-value = 0.1043
prop.test(c(DAY_b1,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_b2,DAY_b3),c(3000,3000)) # no diff p-value = 1
prop.test(c(DAY_b2,DAY_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_b2,DAY_u),c(3000,3000)) # no diff p-value = 0.2251
prop.test(c(DAY_b2,DAY_cu),c(3000,3000)) # sign diff p-value = 0.02675
prop.test(c(DAY_b2,DAY_cud),c(3000,3000)) # sign diff p-value = 0.03426
prop.test(c(DAY_b2,DAY_NJ),c(3000,3000)) # sign diff p-value = 0.02357
prop.test(c(DAY_b2,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_b3,DAY_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_b3,DAY_u),c(3000,3000)) # no diff p-value = 0.205
prop.test(c(DAY_b3,DAY_cu),c(3000,3000)) # sign diff p-value = 0.02325
prop.test(c(DAY_b3,DAY_cud),c(3000,3000)) # sign diff p-value = 0.02992
prop.test(c(DAY_b3,DAY_NJ),c(3000,3000)) # sign diff p-value = 0.02044
prop.test(c(DAY_b3,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_b10,DAY_u),c(3000,3000)) # sign diff  p-value < 2.2e-16
prop.test(c(DAY_b10,DAY_cu),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_b10,DAY_cud),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_b10,DAY_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_b10,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_u,DAY_cu),c(3000,3000)) # no diff p-value = 0.3406
prop.test(c(DAY_u,DAY_cud),c(3000,3000)) # no diff p-value = 0.3928
prop.test(c(DAY_u,DAY_NJ),c(3000,3000)) # no diff p-value = 0.3163
prop.test(c(DAY_u,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_cu,DAY_cud),c(3000,3000)) # no diff p-value = 0.9607
prop.test(c(DAY_cu,DAY_NJ),c(3000,3000)) # no diff p-value = 1
prop.test(c(DAY_cu,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_cud,DAY_NJ),c(3000,3000)) # no diff p-value = 0.9217
prop.test(c(DAY_cud,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_NJ,DAY_ML),c(3000,3000)) # sign diff  p-value < 2.2e-16
