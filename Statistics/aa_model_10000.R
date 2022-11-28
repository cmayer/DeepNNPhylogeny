# model 10000 JTT
JTT_mod_10000_b1 = 954 + 973 + 970
JTT_mod_10000_b1 # 2897
JTT_mod_10000_b2 = 966 + 983 + 976 
JTT_mod_10000_b2 # 2925
JTT_mod_10000_b3 = 983 + 983 + 981 
JTT_mod_10000_b3 # 2947
JTT_mod_10000_b10 = 894 + 940 + 926 
JTT_mod_10000_b10 # 2760
JTT_mod_10000_u = 974 + 975 + 977
JTT_mod_10000_u # 2926
JTT_mod_10000_cu = 976 + 966 + 975 
JTT_mod_10000_cu # 2917
JTT_mod_10000_cud = 942 + 980 + 964
JTT_mod_10000_cud # 2886
JTT_NJ = 966 + 907 + 904 
JTT_NJ # 2777
JTT_ML = 1000 + 999 + 999
JTT_ML # 2998 

JTT_mod_10000_b1/30 # 96.57%
JTT_mod_10000_b2/30 # 97.5%
JTT_mod_10000_b3/30 # 98.23%
JTT_mod_10000_b10/30 # 92%
JTT_mod_10000_u/30 # 97.53%
JTT_mod_10000_cu/30 # 97.23%
JTT_mod_10000_cud/30 # 96.2%
JTT_NJ/30 # 92.57%
JTT_ML/30 # 99.93%

# JTT model 10000 binomial test 

prop.test(c(JTT_mod_10000_b1 ,JTT_mod_10000_b2),c(3000,3000)) # sign diff p-value = 0.03993
prop.test(c(JTT_mod_10000_b1 ,JTT_mod_10000_b3),c(3000,3000)) # sign diff p-value = 7.033e-05
prop.test(c(JTT_mod_10000_b1 ,JTT_mod_10000_b10),c(3000,3000)) # sign diff p-value = 3.949e-14
prop.test(c(JTT_mod_10000_b1 ,JTT_mod_10000_u),c(3000,3000)) # sign diff p-value = 0.03265
prop.test(c(JTT_mod_10000_b1 ,JTT_mod_10000_cu),c(3000,3000)) # no diff p-value = 0.157
prop.test(c(JTT_mod_10000_b1 ,JTT_mod_10000_cud),c(3000,3000)) # no diff p-value = 0.4893
prop.test(c(JTT_mod_10000_b1 ,JTT_NJ),c(3000,3000)) # sign diff p-value = 1.223e-11
prop.test(c(JTT_mod_10000_b1 ,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_mod_10000_b2 ,JTT_mod_10000_b3),c(3000,3000)) # almost diff p-value = 0.06062
prop.test(c(JTT_mod_10000_b2 ,JTT_mod_10000_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_mod_10000_b2 ,JTT_mod_10000_u),c(3000,3000)) # no diff  p-value = 1
prop.test(c(JTT_mod_10000_b2 ,JTT_mod_10000_cu),c(3000,3000)) # no diff p-value = 0.5725
prop.test(c(JTT_mod_10000_b2 ,JTT_mod_10000_cud),c(3000,3000)) # sign diff p-value = 0.004974
prop.test(c(JTT_mod_10000_b2 ,JTT_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_mod_10000_b2 ,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_mod_10000_b3 ,JTT_mod_10000_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_mod_10000_b3 ,JTT_mod_10000_u),c(3000,3000)) # almost diff p-value = 0.07285
prop.test(c(JTT_mod_10000_b3 ,JTT_mod_10000_cu),c(3000,3000)) # sign diff p-value = 0.01189
prop.test(c(JTT_mod_10000_b3 ,JTT_mod_10000_cud),c(3000,3000)) # sign diff p-value = 2.49e-06
prop.test(c(JTT_mod_10000_b3 ,JTT_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_mod_10000_b3 ,JTT_ML),c(3000,3000)) # sign diff p-value = 1.26e-11

prop.test(c(JTT_mod_10000_b10 ,JTT_mod_10000_u),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_mod_10000_b10 ,JTT_mod_10000_cu),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_mod_10000_b10 ,JTT_mod_10000_cud),c(3000,3000)) # sign diff p-value = 7.448e-12
prop.test(c(JTT_mod_10000_b10 ,JTT_NJ),c(3000,3000)) # no diff p-value = 0.4389
prop.test(c(JTT_mod_10000_b10 ,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_mod_10000_u ,JTT_mod_10000_cu),c(3000,3000)) # no diff p-value = 0.5176
prop.test(c(JTT_mod_10000_u ,JTT_mod_10000_cud),c(3000,3000)) # sign diff p-value = 0.003852
prop.test(c(JTT_mod_10000_u ,JTT_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_mod_10000_u ,JTT_ML),c(3000,3000)) # sign diff p-value = 2.478e-16

prop.test(c(JTT_mod_10000_cu ,JTT_mod_10000_cud),c(3000,3000)) # sign diff p-value = 0.02975
prop.test(c(JTT_mod_10000_cu ,JTT_NJ),c(3000,3000)) # sign diff p-value = 3.44e-16
prop.test(c(JTT_mod_10000_cu ,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_mod_10000_cud ,JTT_NJ),c(3000,3000)) # sign diff p-value = 1.398e-09
prop.test(c(JTT_mod_10000_cud ,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_NJ ,JTT_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# model 10000 WAG
WAG_mod_10000_b1 = 970 + 964 + 958
WAG_mod_10000_b1 # 2892 
WAG_mod_10000_b2 = 988 + 980 + 966
WAG_mod_10000_b2 # 2934
WAG_mod_10000_b3 = 988 + 974 + 985
WAG_mod_10000_b3 # 2947
WAG_mod_10000_b10 = 902 + 923 + 912 
WAG_mod_10000_b10 # 2737
WAG_mod_10000_u = 960 + 982 + 971 
WAG_mod_10000_u # 2913
WAG_mod_10000_cu = 964 + 960 + 965
WAG_mod_10000_cu # 2889
WAG_mod_10000_cud = 937 + 969 + 968
WAG_mod_10000_cud # 2874
WAG_NJ = 964 + 892 + 891 
WAG_NJ # 2747
WAG_ML = 1000 + 999 + 1000
WAG_ML # 2999

WAG_mod_10000_b1/30 # 96.4%
WAG_mod_10000_b2/30 # 97.8%
WAG_mod_10000_b3/30 # 98.23%
WAG_mod_10000_b10/30 # 91.23%
WAG_mod_10000_u/30 # 97.1%
WAG_mod_10000_cu/30 # 96.3%
WAG_mod_10000_cud/30 # 95.8%
WAG_NJ/30 # 91.57%
WAG_ML/30 # 99.97%

# JTT model 10000 binomial test 

prop.test(c(WAG_mod_10000_b1 ,WAG_mod_10000_b2),c(3000,3000)) # sign diff p-value = 0.001609
prop.test(c(WAG_mod_10000_b1 ,WAG_mod_10000_b3),c(3000,3000)) # sign diff p-value = 0.001609
prop.test(c(WAG_mod_10000_b1 ,WAG_mod_10000_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_mod_10000_b1 ,WAG_mod_10000_u),c(3000,3000)) # no diff p-value = 0.1454
prop.test(c(WAG_mod_10000_b1 ,WAG_mod_10000_cu),c(3000,3000)) # no diff  p-value = 0.8905
prop.test(c(WAG_mod_10000_b1 ,WAG_mod_10000_cud),c(3000,3000)) # no diff p-value = 0.2569
prop.test(c(WAG_mod_10000_b1 ,WAG_NJ),c(3000,3000)) # sign diff p-value = 5.376e-15
prop.test(c(WAG_mod_10000_b1 ,WAG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_mod_10000_b2 ,WAG_mod_10000_b3),c(3000,3000)) # no diff p-value = 0.2665
prop.test(c(WAG_mod_10000_b2 ,WAG_mod_10000_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_mod_10000_b2 ,WAG_mod_10000_u),c(3000,3000)) # no diff p-value = 0.1014
prop.test(c(WAG_mod_10000_b2 ,WAG_mod_10000_cu),c(3000,3000)) # sign diff p-value = 0.0007876
prop.test(c(WAG_mod_10000_b2 ,WAG_mod_10000_cud),c(3000,3000)) # sign diff p-value = 1.506e-05
prop.test(c(WAG_mod_10000_b2 ,WAG_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_mod_10000_b2 ,WAG_ML),c(3000,3000)) # sign diff p-value = 3.754e-15

prop.test(c(WAG_mod_10000_b3 ,WAG_mod_10000_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_mod_10000_b3 ,WAG_mod_10000_u),c(3000,3000)) # sign diff p-value = 0.004771
prop.test(c(WAG_mod_10000_b3 ,WAG_mod_10000_cu),c(3000,3000)) # sign diff p-value = 6.39e-06
prop.test(c(WAG_mod_10000_b3 ,WAG_mod_10000_cud),c(3000,3000)) # sign diff p-value = 4.664e-08
prop.test(c(WAG_mod_10000_b3 ,WAG_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_mod_10000_b3 ,WAG_ML),c(3000,3000)) # sign diff p-value = 3.132e-12

prop.test(c(WAG_mod_10000_b10 ,WAG_mod_10000_u),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_mod_10000_b10 ,WAG_mod_10000_cu),c(3000,3000)) # sign diff p-value = 7.422e-16
prop.test(c(WAG_mod_10000_b10 ,WAG_mod_10000_cud),c(3000,3000)) # sign diff p-value = 1e-12
prop.test(c(WAG_mod_10000_b10 ,WAG_NJ),c(3000,3000)) # no diff p-value = 0.6786
prop.test(c(WAG_mod_10000_b10 ,WAG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_mod_10000_u ,WAG_mod_10000_cu),c(3000,3000)) # almost diff p-value = 0.09647
prop.test(c(WAG_mod_10000_u ,WAG_mod_10000_cud),c(3000,3000)) # sign diff p-value = 0.008021
prop.test(c(WAG_mod_10000_u ,WAG_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_mod_10000_u ,WAG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_mod_10000_cu ,WAG_mod_10000_cud),c(3000,3000)) # no diff p-value = 0.3535
prop.test(c(WAG_mod_10000_cu ,WAG_NJ),c(3000,3000)) # sign diff p-value = 2.434e-14
prop.test(c(WAG_mod_10000_cu ,WAG_ML),c(3000,3000)) # sign diff  p-value < 2.2e-16

prop.test(c(WAG_mod_10000_cud ,WAG_NJ),c(3000,3000)) # sign diff p-value = 2.281e-11
prop.test(c(WAG_mod_10000_cud ,WAG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_NJ ,WAG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# model 10000 LG
LG_mod_10000_b1 = 953 + 964 + 974 
LG_mod_10000_b1 # 2891
LG_mod_10000_b2 = 978 + 968 + 972 
LG_mod_10000_b2 # 2918
LG_mod_10000_b3 = 988 + 975 + 987 
LG_mod_10000_b3 # 2950
LG_mod_10000_b10 = 907 + 933 + 916
LG_mod_10000_b10 # 2756
LG_mod_10000_u = 986 + 973 + 970 
LG_mod_10000_u # 2929
LG_mod_10000_cu = 961 + 971 + 968
LG_mod_10000_cu # 2900
LG_mod_10000_cud = 968 + 961 + 976 
LG_mod_10000_cud # 2905
LG_NJ = 970 + 910 + 910
LG_NJ # 2790
LG_ML = 1000 + 999 + 998 
LG_ML # 2997

LG_mod_10000_b1/30 # 96.37%
LG_mod_10000_b2/30 # 97.27%
LG_mod_10000_b3/30 # 98.33%
LG_mod_10000_b10/30 # 91.87%
LG_mod_10000_u/30 # 97.63%
LG_mod_10000_cu/30 # 96.67%
LG_mod_10000_cud/30 # 96.83%
LG_NJ/30 # 93%
LG_ML/30 # 99.9%

# LG model 10000 binomial test 

prop.test(c(LG_mod_10000_b1 ,LG_mod_10000_b2),c(3000,3000)) # almost sign p-value = 0.05588
prop.test(c(LG_mod_10000_b1 ,LG_mod_10000_b3),c(3000,3000)) # sign diff p-value = 3.133e-06
prop.test(c(LG_mod_10000_b1 ,LG_mod_10000_b10),c(3000,3000)) # sign diff p-value = 1.958e-13
prop.test(c(LG_mod_10000_b1 ,LG_mod_10000_u),c(3000,3000)) # sign diff p-value = 0.005108
prop.test(c(LG_mod_10000_b1 ,LG_mod_10000_cu),c(3000,3000)) # no diff p-value = 0.5733
prop.test(c(LG_mod_10000_b1 ,LG_mod_10000_cud),c(3000,3000)) # no diff p-value = 0.3544
prop.test(c(LG_mod_10000_b1 ,LG_NJ),c(3000,3000)) # sign diff p-value = 8.717e-09
prop.test(c(LG_mod_10000_b1 ,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_mod_10000_b2 ,LG_mod_10000_b3),c(3000,3000)) # sign diff p-value = 0.006365
prop.test(c(LG_mod_10000_b2 ,LG_mod_10000_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_mod_10000_b2 ,LG_mod_10000_u),c(3000,3000)) # no diff p-value = 0.4128
prop.test(c(LG_mod_10000_b2 ,LG_mod_10000_cu),c(3000,3000)) # no diff p-value = 0.2007
prop.test(c(LG_mod_10000_b2 ,LG_mod_10000_cud),c(3000,3000)) # no diff p-value = 0.3599
prop.test(c(LG_mod_10000_b2 ,LG_NJ),c(3000,3000)) # sign diff p-value = 2.54e-14
prop.test(c(LG_mod_10000_b2 ,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_mod_10000_b3 ,LG_mod_10000_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_mod_10000_b3 ,LG_mod_10000_u),c(3000,3000)) # almost diff p-value = 0.06624
prop.test(c(LG_mod_10000_b3 ,LG_mod_10000_cu),c(3000,3000)) # sign diff p-value = 5.082e-05
prop.test(c(LG_mod_10000_b3 ,LG_mod_10000_cud),c(3000,3000)) # sign diff p-value = 0.0002165
prop.test(c(LG_mod_10000_b3 ,LG_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_mod_10000_b3 ,LG_ML),c(3000,3000)) # sign diff p-value = 2.2e-10

prop.test(c(LG_mod_10000_b10 ,LG_mod_10000_u),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_mod_10000_b10 ,LG_mod_10000_cu),c(3000,3000)) # sign diff p-value = 2.005e-15
prop.test(c(LG_mod_10000_b10 ,LG_mod_10000_cud),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_mod_10000_b10 ,LG_NJ),c(3000,3000)) # no diff  p-value = 0.1072
prop.test(c(LG_mod_10000_b10 ,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_mod_10000_u ,LG_mod_10000_cu),c(3000,3000)) # sign diff p-value = 0.02983
prop.test(c(LG_mod_10000_u ,LG_mod_10000_cud),c(3000,3000)) # almost diff p-value = 0.07024
prop.test(c(LG_mod_10000_u ,LG_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_mod_10000_u ,LG_ML),c(3000,3000)) # sign diff p-value = 4.612e-15

prop.test(c(LG_mod_10000_cu ,LG_mod_10000_cud),c(3000,3000)) # no diff p-value = 0.7709
prop.test(c(LG_mod_10000_cu ,LG_NJ),c(3000,3000)) # sign diff p-value = 2.055e-10
prop.test(c(LG_mod_10000_cu ,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_mod_10000_cud ,LG_NJ),c(3000,3000)) # sign diff p-value = 2.082e-11
prop.test(c(LG_mod_10000_cud ,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_NJ ,LG_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# model 10000 DAY
DAY_mod_10000_b1 =  960 + 959 + 959 
DAY_mod_10000_b1 # 2878
DAY_mod_10000_b2 =  968 + 977 + 961 
DAY_mod_10000_b2 # 2906
DAY_mod_10000_b3 =  977 + 984 + 965 
DAY_mod_10000_b3 # 2926
DAY_mod_10000_b10 =  902 + 908 + 918 
DAY_mod_10000_b10 # 2728
DAY_mod_10000_u =  963 + 980 + 965 
DAY_mod_10000_u # 2908
DAY_mod_10000_cu =  943 + 982 + 951 
DAY_mod_10000_cu # 2876
DAY_mod_10000_cud =  955 + 976 + 947 
DAY_mod_10000_cud # 2878
DAY_NJ = 961 + 900 + 914
DAY_NJ # 2775
DAY_ML = 999 + 999 +999
DAY_ML #   2997

DAY_mod_10000_b1/30 # 95.93%
DAY_mod_10000_b2/30 # 96.87%
DAY_mod_10000_b3/30 # 97.53%
DAY_mod_10000_b10/30 #  90.93%
DAY_mod_10000_u/30 # 96.93%
DAY_mod_10000_cu/30 # 95.87%
DAY_mod_10000_cud/30 # 95.93%
DAY_NJ/30 # 92.5%
DAY_ML/30 # 99.9%

# DAY model 10000 binomial test 

prop.test(c(DAY_mod_10000_b1 ,DAY_mod_10000_b2),c(3000,3000)) # almost diff p-value = 0.06133
prop.test(c(DAY_mod_10000_b1 ,DAY_mod_10000_b3),c(3000,3000)) # sign diff p-value = 0.0006417
prop.test(c(DAY_mod_10000_b1 ,DAY_mod_10000_b10),c(3000,3000)) # sign diff p-value = 8.112e-15
prop.test(c(DAY_mod_10000_b1 ,DAY_mod_10000_u),c(3000,3000)) # sign diff p-value = 0.04352
prop.test(c(DAY_mod_10000_b1 ,DAY_mod_10000_cu),c(3000,3000)) # no diff p-value = 0.9481
prop.test(c(DAY_mod_10000_b1 ,DAY_mod_10000_cud),c(3000,3000)) # no diff p-value = 1
prop.test(c(DAY_mod_10000_b1 ,DAY_NJ),c(3000,3000)) # sign diff p-value = 1.689e-08
prop.test(c(DAY_mod_10000_b1 ,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_mod_10000_b2 ,DAY_mod_10000_b3),c(3000,3000)) # no diff p-value = 0.1371
prop.test(c(DAY_mod_10000_b2 ,DAY_mod_10000_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_mod_10000_b2 ,DAY_mod_10000_u),c(3000,3000)) # no diff p-value = 0.9406
prop.test(c(DAY_mod_10000_b2 ,DAY_mod_10000_cu),c(3000,3000)) # sign diff p-value = 0.04541
prop.test(c(DAY_mod_10000_b2 ,DAY_mod_10000_cud),c(3000,3000)) # almost diff p-value = 0.06133
prop.test(c(DAY_mod_10000_b2 ,DAY_NJ),c(3000,3000)) # sign diff p-value = 7.423e-14
prop.test(c(DAY_mod_10000_b2 ,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_mod_10000_b3 ,DAY_mod_10000_b10),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_mod_10000_b3 ,DAY_mod_10000_u),c(3000,3000)) # no diff p-value = 0.1809
prop.test(c(DAY_mod_10000_b3 ,DAY_mod_10000_cu),c(3000,3000)) # sign diff p-value = 0.0003983
prop.test(c(DAY_mod_10000_b3 ,DAY_mod_10000_cud),c(3000,3000)) # sign diff p-value = 0.0006417
prop.test(c(DAY_mod_10000_b3 ,DAY_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_mod_10000_b3 ,DAY_ML),c(3000,3000)) # sign diff p-value = 9.833e-16

prop.test(c(DAY_mod_10000_b10 ,DAY_mod_10000_u),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_mod_10000_b10 ,DAY_mod_10000_cu),c(3000,3000)) # sign diff p-value = 2.113e-14
prop.test(c(DAY_mod_10000_b10 ,DAY_mod_10000_cud),c(3000,3000)) # sign diff p-value = 8.112e-15
prop.test(c(DAY_mod_10000_b10 ,DAY_NJ),c(3000,3000)) # sign diff p-value = 0.0312
prop.test(c(DAY_mod_10000_b10 ,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_mod_10000_u ,DAY_mod_10000_cu),c(3000,3000)) # sign diff p-value = 0.03169
prop.test(c(DAY_mod_10000_u ,DAY_mod_10000_cud),c(3000,3000)) # sign diff p-value = 0.04352
prop.test(c(DAY_mod_10000_u ,DAY_NJ),c(3000,3000)) # sign diff p-value = 2.58e-14
prop.test(c(DAY_mod_10000_u ,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_mod_10000_cu ,DAY_mod_10000_cud),c(3000,3000)) # no diff p-value = 0.9481
prop.test(c(DAY_mod_10000_cu ,DAY_NJ),c(3000,3000)) # sign diff p-value = 3.474e-08
prop.test(c(DAY_mod_10000_cu ,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_mod_10000_cud ,DAY_NJ),c(3000,3000)) # sign diff p-value = 1.689e-08
prop.test(c(DAY_mod_10000_cud ,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_NJ ,DAY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16
