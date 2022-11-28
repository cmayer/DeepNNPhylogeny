# JTT 10000 succ. rate
JTT_10000_b1 = 998 + 994 + 994 
JTT_10000_b1 # 2986
JTT_10000_b2 = 999 + 994 + 996
JTT_10000_b2 # 2989
JTT_10000_b3 = 1000 + 988 + 997 
JTT_10000_b3 # 2985 
JTT_10000_b10 = 993 + 989 + 980
JTT_10000_b10 # 2962
JTT_10000_u = 996 + 990 + 994 
JTT_10000_u # 2980
JTT_10000_cu = 994 + 997 + 987 
JTT_10000_cu # 2978
JTT_10000_cud = 996 + 994 + 987 
JTT_10000_cud # 2977
JTT_10000_NJ = 971 + 915 + 922 
JTT_10000_NJ #  2808
JTT_10000_ML = 1000 + 1000 + 1000
JTT_10000_ML #  3000

JTT_10000_b1/30 # 99.53%
JTT_10000_b2/30 # 99.63%
JTT_10000_b3/30 # 99.5%
JTT_10000_b10/30 # 98.73%
JTT_10000_u/30 # 99.33%
JTT_10000_cu/30 # 99.27%
JTT_10000_cud/30 # 99.23%
JTT_10000_NJ/30 # 93.6%
JTT_10000_ML/30 # 100%

# JTT 10000 binomial test 

prop.test(c(JTT_10000_b1 ,JTT_10000_b2),c(3000,3000)) # no diff p-value = 0.6885
prop.test(c(JTT_10000_b1 ,JTT_10000_b3),c(3000,3000)) # no diff p-value = 1
prop.test(c(JTT_10000_b1 ,JTT_10000_b10),c(3000,3000)) # sign diff p-value = 0.001358
prop.test(c(JTT_10000_b1 ,JTT_10000_u),c(3000,3000)) # no diff p-value = 0.3898
prop.test(c(JTT_10000_b1 ,JTT_10000_cu),c(3000,3000)) # no diff p-value = 0.2419
prop.test(c(JTT_10000_b1 ,JTT_10000_cud),c(3000,3000)) # no diff p-value = 0.1871
prop.test(c(JTT_10000_b1 ,JTT_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_10000_b1 ,JTT_10000_ML),c(3000,3000)) # sign diff p-value = 0.0005043

prop.test(c(JTT_10000_b2 ,JTT_10000_b3),c(3000,3000)) # no diff p-value = 0.5554
prop.test(c(JTT_10000_b2 ,JTT_10000_b10),c(3000,3000)) # sign diff p-value = 0.0001918
prop.test(c(JTT_10000_b2 ,JTT_10000_u),c(3000,3000)) # no diff p-value = 0.1497
prop.test(c(JTT_10000_b2 ,JTT_10000_cu),c(3000,3000)) # almost diff p-value = 0.08088
prop.test(c(JTT_10000_b2 ,JTT_10000_cud),c(3000,3000)) # almost diff p-value = 0.05851
prop.test(c(JTT_10000_b2 ,JTT_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_10000_b2 ,JTT_10000_ML),c(3000,3000)) # sign diff p-value = 0.002545

prop.test(c(JTT_10000_b3 ,JTT_10000_b10),c(3000,3000)) # sign diff p-value = 0.002402
prop.test(c(JTT_10000_b3 ,JTT_10000_u),c(3000,3000)) # no diff p-value = 0.4977
prop.test(c(JTT_10000_b3 ,JTT_10000_cu),c(3000,3000)) # no diff p-value = 0.3224
prop.test(c(JTT_10000_b3 ,JTT_10000_cud),c(3000,3000)) # no diff p-value = 0.2546
prop.test(c(JTT_10000_b3 ,JTT_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_10000_b3 ,JTT_10000_ML),c(3000,3000)) # sign diff p-value = 0.0002954

prop.test(c(JTT_10000_b10 ,JTT_10000_u),c(3000,3000)) # sign diff p-value = 0.02489
prop.test(c(JTT_10000_b10 ,JTT_10000_cu),c(3000,3000)) # almost diff  p-value = 0.05163
prop.test(c(JTT_10000_b10 ,JTT_10000_cud),c(3000,3000)) # almost diff p-value = 0.07159
prop.test(c(JTT_10000_b10 ,JTT_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_10000_b10 ,JTT_10000_ML),c(3000,3000)) # sign diff p-value = 1.73e-09

prop.test(c(JTT_10000_u ,JTT_10000_cu),c(3000,3000)) # no diff p-value = 0.8769
prop.test(c(JTT_10000_u ,JTT_10000_cud),c(3000,3000)) # no diff p-value = 0.7595
prop.test(c(JTT_10000_u ,JTT_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_10000_u ,JTT_10000_ML),c(3000,3000)) # sign diff p-value = 2.085e-05

prop.test(c(JTT_10000_cu ,JTT_10000_cud),c(3000,3000)) # no diff  p-value = 1
prop.test(c(JTT_10000_cu ,JTT_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_10000_cu ,JTT_10000_ML),c(3000,3000)) # sign diff p-value = 7.276e-06

prop.test(c(JTT_10000_cud ,JTT_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_10000_cud ,JTT_10000_ML),c(3000,3000)) # sign diff p-value = 4.304e-06

prop.test(c(JTT_10000_NJ ,JTT_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# LG 10000 succ. rate
LG_10000_b1 = 995 + 993 + 994 
LG_10000_b1 #
LG_10000_b2 = 999 + 992 + 994 
LG_10000_b2 #
LG_10000_b3 = 999 + 992 + 993 
LG_10000_b3 # 2984
LG_10000_b10 = 996 + 991 + 990
LG_10000_b10 # 2977
LG_10000_u = 998 + 992 + 996 
LG_10000_u # 2986
LG_10000_cu = 991 + 993 + 995 
LG_10000_cu # 2979
LG_10000_cud = 998 + 990 + 989
LG_10000_cud # 2977
LG_10000_NJ = 978 + 935 + 929 
LG_10000_NJ # 2842 
LG_10000_ML = 1000 + 1000 + 1000
LG_10000_ML #  3000

LG_10000_b1/30 # 99.4%
LG_10000_b2/30 # 99.5%
LG_10000_b3/30 # 99.47%
LG_10000_b10/30 # 99.23%
LG_10000_u/30 # 99.53%
LG_10000_cu/30 # 99.3%
LG_10000_cud/30 # 99.23%
LG_10000_NJ/30 # 94.73%
LG_10000_ML/30 # 100%

models <- c(LG_10000_b3,
      LG_10000_b10,
      LG_10000_u,
      LG_10000_cu,
      LG_10000_cud)
# LG 10000 binomial test 

prop.test(c(LG_10000_b1 ,LG_10000_b2),c(3000,3000)) # no diff  p-value = 0.727
prop.test(c(LG_10000_b1 ,LG_10000_b3),c(3000,3000)) # no diff p-value = 0.8634
prop.test(c(LG_10000_b1 ,LG_10000_b10),c(3000,3000)) # no diff p-value = 0.5308
prop.test(c(LG_10000_b1 ,LG_10000_u),c(3000,3000)) # no diff p-value = 0.5949
prop.test(c(LG_10000_b1 ,LG_10000_cu),c(3000,3000)) # no diff p-value = 0.748
prop.test(c(LG_10000_b1 ,LG_10000_cud),c(3000,3000)) # no diff p-value = 0.5308
prop.test(c(LG_10000_b1 ,LG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_10000_b1 ,LG_10000_ML),c(3000,3000)) # sign diff p-value = 5.996e-05

prop.test(c(LG_10000_b2 ,LG_10000_b3),c(3000,3000)) # no diff p-value = 1
prop.test(c(LG_10000_b2 ,LG_10000_b10),c(3000,3000)) # no diff p-value = 0.2546
prop.test(c(LG_10000_b2 ,LG_10000_u),c(3000,3000)) # no diff p-value = 1
prop.test(c(LG_10000_b2 ,LG_10000_cu),c(3000,3000)) # no diff p-value = 0.4032
prop.test(c(LG_10000_b2 ,LG_10000_cud),c(3000,3000)) # no diff p-value = 0.2546
prop.test(c(LG_10000_b2 ,LG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_10000_b2 ,LG_10000_ML),c(3000,3000)) # sign diff p-value = 0.0002954

prop.test(c(LG_10000_b3 ,LG_10000_b10),c(3000,3000)) # no diff p-value = 0.3351
prop.test(c(LG_10000_b3 ,LG_10000_u),c(3000,3000)) # no diff p-value = 0.8548
prop.test(c(LG_10000_b3 ,LG_10000_cu),c(3000,3000)) # no diff p-value = 0.5095
prop.test(c(LG_10000_b3 ,LG_10000_cud),c(3000,3000)) # no diff p-value = 0.3351
prop.test(c(LG_10000_b3 ,LG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_10000_b3 ,LG_10000_ML),c(3000,3000)) # sign diff p-value = 0.0001733

prop.test(c(LG_10000_b10 ,LG_10000_u),c(3000,3000)) # no diff p-value = 0.1871
prop.test(c(LG_10000_b10 ,LG_10000_cu),c(3000,3000)) # no diff p-value = 0.8797
prop.test(c(LG_10000_b10 ,LG_10000_cud),c(3000,3000)) # no diff p-value = 1
prop.test(c(LG_10000_b10 ,LG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_10000_b10 ,LG_10000_ML),c(3000,3000)) # sign diff p-value = 4.304e-06

prop.test(c(LG_10000_u ,LG_10000_cu),c(3000,3000)) # no diff p-value = 0.3091
prop.test(c(LG_10000_u ,LG_10000_cud),c(3000,3000)) # no diff  p-value = 0.1871
prop.test(c(LG_10000_u ,LG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_10000_u ,LG_10000_ML),c(3000,3000)) # sign diff p-value = 0.0005043

prop.test(c(LG_10000_cu ,LG_10000_cud),c(3000,3000)) # no diff p-value = 0.8797
prop.test(c(LG_10000_cu ,LG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_10000_cu ,LG_10000_ML),c(3000,3000)) # sign diff p-value = 1.231e-05

prop.test(c(LG_10000_cud ,LG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_10000_cud ,LG_10000_ML),c(3000,3000)) # sign diff p-value = 4.304e-06

prop.test(c(LG_10000_NJ ,LG_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# WAG 10000 succ. rate
WAG_10000_b1 = 988 + 990 + 996
WAG_10000_b1 # 2974
WAG_10000_b2 = 992 + 993 + 997
WAG_10000_b2 # 2982
WAG_10000_b3 = 997 + 990 + 993 
WAG_10000_b3 # 2980
WAG_10000_b10 = 979 + 988 + 985
WAG_10000_b10 # 2952
WAG_10000_u = 993 + 990 + 990 
WAG_10000_u # 2973
WAG_10000_cu = 991 + 989 + 987
WAG_10000_cu # 2967
WAG_10000_cud = 989 + 983 + 985
WAG_10000_cud # 2957
WAG_10000_NJ = 965 + 913 + 914 
WAG_10000_NJ # 2792 
WAG_10000_ML = 1000 + 1000 + 1000
WAG_10000_ML #  3000

WAG_10000_b1/30 # 99.13%
WAG_10000_b2/30 # 99.4%
WAG_10000_b3/30 # 99.33%
WAG_10000_b10/30 # 98.4%
WAG_10000_u/30 # 99.1%
WAG_10000_cu/30 # 98.9%
WAG_10000_cud/30 # 98.57%
WAG_10000_NJ/30 # 93.07%
WAG_10000_ML/30 # 100% 

# WAG 10000 binomial test 

prop.test(c(WAG_10000_b1 ,WAG_10000_b2),c(3000,3000)) # no diff p-value = 0.2895
prop.test(c(WAG_10000_b1 ,WAG_10000_b3),c(3000,3000)) # no diff p-value = 0.4593
prop.test(c(WAG_10000_b1 ,WAG_10000_b10),c(3000,3000)) # sign diff p-value = 0.01403
prop.test(c(WAG_10000_b1 ,WAG_10000_u),c(3000,3000)) # no diff p-value = 1
prop.test(c(WAG_10000_b1 ,WAG_10000_cu),c(3000,3000)) # no diff p-value = 0.4325
prop.test(c(WAG_10000_b1 ,WAG_10000_cud),c(3000,3000)) # almost diff p-value = 0.0527
prop.test(c(WAG_10000_b1 ,WAG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_10000_b1 ,WAG_10000_ML),c(3000,3000)) # sign diff p-value = 8.944e-07

prop.test(c(WAG_10000_b2 ,WAG_10000_b3),c(3000,3000)) # no diff p-value = 0.8707
prop.test(c(WAG_10000_b2 ,WAG_10000_b10),c(3000,3000)) # sign diff p-value = 0.0003314
prop.test(c(WAG_10000_b2 ,WAG_10000_u),c(3000,3000)) # no diff p-value = 0.2313
prop.test(c(WAG_10000_b2 ,WAG_10000_cu),c(3000,3000)) # sign diff p-value = 0.04898
prop.test(c(WAG_10000_b2 ,WAG_10000_cud),c(3000,3000)) # sign diff p-value = 0.002011
prop.test(c(WAG_10000_b2 ,WAG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_10000_b2 ,WAG_10000_ML),c(3000,3000)) # sign diff p-value = 5.996e-05

prop.test(c(WAG_10000_b3 ,WAG_10000_b10),c(3000,3000)) # sign diff p-value = 0.0009914
prop.test(c(WAG_10000_b3 ,WAG_10000_u),c(3000,3000)) # no diff p-value = 0.3796
prop.test(c(WAG_10000_b3 ,WAG_10000_cu),c(3000,3000)) # almost diff p-value = 0.09779
prop.test(c(WAG_10000_b3 ,WAG_10000_cud),c(3000,3000)) # sign diff p-value = 0.00533
prop.test(c(WAG_10000_b3 ,WAG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_10000_b3 ,WAG_10000_ML),c(3000,3000)) # sign diff p-value = 2.085e-05

prop.test(c(WAG_10000_b10 ,WAG_10000_u),c(3000,3000)) # sign diff p-value = 0.02013
prop.test(c(WAG_10000_b10 ,WAG_10000_cu),c(3000,3000)) # no diff p-value = 0.1173
prop.test(c(WAG_10000_b10 ,WAG_10000_cud),c(3000,3000)) # no diff p-value = 0.6726
prop.test(c(WAG_10000_b10 ,WAG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_10000_b10 ,WAG_10000_ML),c(3000,3000)) # sign diff p-value = 9.681e-12

prop.test(c(WAG_10000_u ,WAG_10000_cu),c(3000,3000)) # no diff p-value = 0.5165
prop.test(c(WAG_10000_u ,WAG_10000_cud),c(3000,3000)) # almost diff p-value = 0.07133
prop.test(c(WAG_10000_u ,WAG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_10000_u ,WAG_10000_ML),c(3000,3000)) # sign diff p-value = 5.303e-07

prop.test(c(WAG_10000_cu ,WAG_10000_cud),c(3000,3000)) # no diff p-value = 0.2988
prop.test(c(WAG_10000_cu ,WAG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_10000_cu ,WAG_10000_ML),c(3000,3000)) # sign diff p-value = 2.325e-08

prop.test(c(WAG_10000_cud ,WAG_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_10000_cud ,WAG_10000_ML),c(3000,3000)) # sign diff p-value = 1.293e-10

prop.test(c(WAG_10000_NJ ,WAG_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# DAY 10000 succ. rate
DAY_10000_b1 = 994 + 997 + 996 
DAY_10000_b1 # 2987
DAY_10000_b2 = 997 + 991 + 996
DAY_10000_b2 # 2984 
DAY_10000_b3 = 997 + 995 + 998 
DAY_10000_b3 # 2990
DAY_10000_b10 = 986 + 986 + 993 
DAY_10000_b10 # 2965
DAY_10000_u = 994 + 990 + 994 
DAY_10000_u # 2978
DAY_10000_cu = 999 + 993 + 989
DAY_10000_cu # 2981
DAY_10000_cud = 995 + 990 + 997 
DAY_10000_cud #  2982
DAY_10000_NJ = 973 + 921 + 920 
DAY_10000_NJ # 2814 
DAY_10000_ML = 1000 + 1000 + 1000
DAY_10000_ML #  3000

DAY_10000_b1/30 # 99.57%
DAY_10000_b2/30 # 99.47%
DAY_10000_b3/30 # 99.67%
DAY_10000_b10/30 # 98.83%
DAY_10000_u/30 # 99.27%
DAY_10000_cu/30 # 99.37%
DAY_10000_cud/30 # 99.4%
DAY_10000_NJ/30 # 93.8%
DAY_10000_ML/30 # 100%

# DAY 10000 binomial test 
prop.test(c(DAY_10000_b1 ,DAY_10000_b2),c(3000,3000)) # no diff p-value = 0.7097
prop.test(c(DAY_10000_b1 ,DAY_10000_b3),c(3000,3000)) # no diff p-value = 0.6761
prop.test(c(DAY_10000_b1 ,DAY_10000_b10),c(3000,3000)) # sign diff p-value = 0.00234
prop.test(c(DAY_10000_b1 ,DAY_10000_u),c(3000,3000)) # no diff p-value = 0.175
prop.test(c(DAY_10000_b1 ,DAY_10000_cu),c(3000,3000)) # no diff p-value = 0.3755
prop.test(c(DAY_10000_b1 ,DAY_10000_cud),c(3000,3000)) # no diff p-value = 0.4713
prop.test(c(DAY_10000_b1 ,DAY_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_10000_b1 ,DAY_10000_ML),c(3000,3000)) # sign diff p-value = 0.0008628

prop.test(c(DAY_10000_b2 ,DAY_10000_b3),c(3000,3000)) # no diff p-value = 0.3257
prop.test(c(DAY_10000_b2 ,DAY_10000_b10),c(3000,3000)) # sign diff p-value = 0.01136
prop.test(c(DAY_10000_b2 ,DAY_10000_u),c(3000,3000)) # no diff p-value = 0.4158
prop.test(c(DAY_10000_b2 ,DAY_10000_cu),c(3000,3000)) # no diff p-value = 0.7346
prop.test(c(DAY_10000_b2 ,DAY_10000_cud),c(3000,3000)) # no diff p-value = 0.8634
prop.test(c(DAY_10000_b2 ,DAY_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_10000_b2 ,DAY_10000_ML),c(3000,3000)) # sign diff p-value = 0.0001733

prop.test(c(DAY_10000_b3 ,DAY_10000_b10),c(3000,3000)) # sign diff p-value = 0.0003292
prop.test(c(DAY_10000_b3 ,DAY_10000_u),c(3000,3000)) # almost diff p-value = 0.05121
prop.test(c(DAY_10000_b3 ,DAY_10000_cu),c(3000,3000)) # no diff p-value = 0.1364
prop.test(c(DAY_10000_b3 ,DAY_10000_cud),c(3000,3000)) # no diff p-value = 0.1848
prop.test(c(DAY_10000_b3 ,DAY_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_10000_b3 ,DAY_10000_ML),c(3000,3000)) # sign diff p-value = 0.004394

prop.test(c(DAY_10000_b10 ,DAY_10000_u),c(3000,3000)) # no diff p-value = 0.1103
prop.test(c(DAY_10000_b10 ,DAY_10000_cu),c(3000,3000)) # sign diff p-value = 0.04032
prop.test(c(DAY_10000_b10 ,DAY_10000_cud),c(3000,3000)) # sign diff p-value = 0.02728
prop.test(c(DAY_10000_b10 ,DAY_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_10000_b10 ,DAY_10000_ML),c(3000,3000)) # sign diff p-value = 8.22e-09

prop.test(c(DAY_10000_u ,DAY_10000_cu),c(3000,3000)) # no diff p-value = 0.754
prop.test(c(DAY_10000_u ,DAY_10000_cud),c(3000,3000)) # no diff p-value = 0.6341
prop.test(c(DAY_10000_u ,DAY_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_10000_u ,DAY_10000_ML),c(3000,3000)) # sign diff p-value = 7.276e-06

prop.test(c(DAY_10000_cu ,DAY_10000_cud),c(3000,3000)) # no diff  p-value = 1
prop.test(c(DAY_10000_cu ,DAY_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_10000_cu ,DAY_10000_ML),c(3000,3000)) # sign diff p-value = 3.534e-05

prop.test(c(DAY_10000_cud ,DAY_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_10000_cud ,DAY_10000_ML),c(3000,3000)) # sign diff p-value = 5.996e-05

prop.test(c(DAY_10000_NJ ,DAY_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

