# Number of succ. JC
JC_b1 = 993 + 993 + 995 
JC_b1 # 2981 99.37%
JC_b2 = 991 + 991 + 996 
JC_b2 # 2978 99.27%
JC_b3 = 988 + 995 + 997 
JC_b3 # 2980 99.33%
JC_b10 = 991 + 998 + 995 
JC_b10 # 2984 99.5%
JC_u = 992 + 998 + 990
JC_u # 2980 99.33%
JC_cu = 989 + 994 + 992 
JC_cu # 2975 99.17%
JC_cud = 994 + 992 + 994 
JC_cud # 2980 99.33%
JC_NJ = 852 + 649 + 673 
JC_NJ # 2174 72.47%
JC_ML = 993 + 995 + 998 
JC_ML # 2986 99.5%

# Binomial test JC 

prop.test(c(JC_b1,JC_b2),c(3000,3000)) # no diff p-value = 0.754
prop.test(c(JC_b1,JC_b3),c(3000,3000)) # no diff  p-value = 1
prop.test(c(JC_b1,JC_b10),c(3000,3000)) # no diff p-value = 0.7346
prop.test(c(JC_b1,JC_u),c(3000,3000)) # no diff  p-value = 1
prop.test(c(JC_b1,JC_cu),c(3000,3000)) # no diff p-value = 0.4493
prop.test(c(JC_b1,JC_cud),c(3000,3000)) # no diff  p-value = 1
prop.test(c(JC_b1,JC_NJ),c(3000,3000)) # sign diff  p-value < 2.2e-16
prop.test(c(JC_b1,JC_ML),c(3000,3000)) # no diff p-value = 0.485

prop.test(c(JC_b2,JC_b3),c(3000,3000)) # no diff p-value = 0.8769
prop.test(c(JC_b2,JC_b10),c(3000,3000)) # no diff  p-value = 0.4158
prop.test(c(JC_b2,JC_u),c(3000,3000)) # no diff p-value = 0.8769
prop.test(c(JC_b2,JC_cu),c(3000,3000)) # no diff p-value = 0.7696
prop.test(c(JC_b2,JC_cud),c(3000,3000)) # no diff  p-value = 0.8769
prop.test(c(JC_b2,JC_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JC_b2,JC_ML),c(3000,3000)) # no diff p-value = 0.2419

prop.test(c(JC_b3,JC_b10),c(3000,3000)) # no diff p-value = 0.616
prop.test(c(JC_b3,JC_u),c(3000,3000)) # no diff p-value = 1
prop.test(c(JC_b3,JC_cu),c(3000,3000)) # no diff p-value = 0.5495
prop.test(c(JC_b3,JC_cud),c(3000,3000)) # no diff p-value = 1
prop.test(c(JC_b3,JC_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JC_b3,JC_ML),c(3000,3000)) # no diff p-value = 0.3898

prop.test(c(JC_b10,JC_u),c(3000,3000)) # no diff p-value = 0.616
prop.test(c(JC_b10,JC_cu),c(3000,3000)) # no diff p-value = 0.21
prop.test(c(JC_b10,JC_cud),c(3000,3000)) # no diff p-value = 0.616
prop.test(c(JC_b10,JC_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JC_b10,JC_ML),c(3000,3000)) # no diff p-value = 0.8548

prop.test(c(JC_u,JC_cu),c(3000,3000)) # no diff p-value = 0.5495
prop.test(c(JC_u,JC_cud),c(3000,3000)) # no diff p-value = 1
prop.test(c(JC_u,JC_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JC_u,JC_ML),c(3000,3000)) # no diff p-value = 0.3898

prop.test(c(JC_cu,JC_cud),c(3000,3000)) # no diff p-value = 0.5495
prop.test(c(JC_cu,JC_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JC_cu,JC_ML),c(3000,3000)) # no diff  p-value = 0.1082

prop.test(c(JC_cud,JC_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JC_cud,JC_ML),c(3000,3000)) # no diff p-value = 0.3898

prop.test(c(JC_NJ,JC_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# Number of succ. K2P 
K2P_b1 = 991 + 996 + 991 
K2P_b1 # 2978
K2P_b1/30 # 99.27%
K2P_b2 = 991 + 996 + 990 
K2P_b2 # 2977
K2P_b2/30 # 99.23%
K2P_b3 = 996 + 993 + 989 
K2P_b3 # 2978
K2P_b3/30 # 99.27%
K2P_b10 = 994 + 993 + 990
K2P_b10 # 2977
K2P_b10/30 # 99.23%
K2P_u = 989 + 986 + 989 
K2P_u # 2964 
K2P_u/30 # 98.8 %
K2P_cu = 995 + 990 + 989 
K2P_cu # 2974 
K2P_cu/30 # 99.13%
K2P_cud = 992 + 989 + 992 
K2P_cud # 2973
K2P_cud/30 #
K2P_NJ = 924 + 862 + 849 
K2P_NJ #   2635
K2P_NJ/30 # 87.83%
K2P_ML = 991 + 996 + 995 
K2P_ML # 2982 
K2P_ML/30 # 99.4% 

# Binomial test K2P 

prop.test(c(K2P_b1,K2P_b2),c(3000,3000)) # no diff p-value = 1
prop.test(c(K2P_b1,K2P_b3),c(3000,3000)) # no diff,  p-value = 1
prop.test(c(K2P_b1,K2P_b10),c(3000,3000)) # no diff, p-value = 1
prop.test(c(K2P_b1,K2P_u),c(3000,3000)) # no diff,   p-value = 0.08629 (close)
prop.test(c(K2P_b1,K2P_cu),c(3000,3000)) # no diff, p-value = 0.6637
prop.test(c(K2P_b1,K2P_cud),c(3000,3000)) # no diff, p-value = 0.5661
prop.test(c(K2P_b1,K2P_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(K2P_b1,K2P_ML),c(3000,3000)) # no diff, p-value = 0.6341

prop.test(c(K2P_b2,K2P_b3),c(3000,3000)) # no diff p-value = 1
prop.test(c(K2P_b2,K2P_b10),c(3000,3000)) # no diff p-value = 1
prop.test(c(K2P_b2,K2P_u),c(3000,3000)) # no diff p-value = 0.1164
prop.test(c(K2P_b2,K2P_cu),c(3000,3000)) # no diff p-value = 0.7742
prop.test(c(K2P_b2,K2P_cud),c(3000,3000)) # no diff p-value = 0.6701
prop.test(c(K2P_b2,K2P_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(K2P_b2,K2P_ML),c(3000,3000)) # no diff p-value = 0.5308

prop.test(c(K2P_b3,K2P_b10),c(3000,3000)) # no diff, p-value = 1
prop.test(c(K2P_b3,K2P_u),c(3000,3000)) # no diff, p-value = 0.08629 (close)
prop.test(c(K2P_b3,K2P_cu),c(3000,3000)) # no diff, p-value = 0.6637
prop.test(c(K2P_b3,K2P_cud),c(3000,3000)) # no diff, p-value = 0.5661
prop.test(c(K2P_b2,K2P_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(K2P_b3,K2P_ML),c(3000,3000)) # no diff, p-value = 0.6341

prop.test(c(K2P_b10,K2P_u),c(3000,3000)) # no diff, p-value = 0.1164
prop.test(c(K2P_b10,K2P_cu),c(3000,3000)) # no diff, p-value = 0.7742
prop.test(c(K2P_b10,K2P_cud),c(3000,3000)) # no diff, p-value = 0.6701
prop.test(c(K2P_b1,K2P_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(K2P_b10,K2P_ML),c(3000,3000)) # no diff,  p-value = 0.5308

prop.test(c(K2P_u,K2P_cu),c(3000,3000)) # no diff, p-value = 0.2506
prop.test(c(K2P_u,K2P_cud),c(3000,3000)) # no diff, p-value = 0.3109
prop.test(c(K2P_b1,K2P_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(K2P_u,K2P_ML),c(3000,3000)) # sign diff, p-value = 0.02013

prop.test(c(K2P_cu,K2P_cud),c(3000,3000)) # no diff, p-value = 1
prop.test(c(K2P_b1,K2P_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(K2P_cu,K2P_ML),c(3000,3000)) # no diff, 0.2895

prop.test(c(K2P_b1,K2P_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(K2P_cud,K2P_ML),c(3000,3000)) # no diff, p-value = 0.2313

prop.test(c(K2P_b1,K2P_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16

# Number of succ. F81. EPOCH = 1000
F81_b1 = 993 + 998 + 997 
F81_b1 # 2988
F81_b2 = 996 + 993 + 994 
F81_b2 # 2983    
F81_b3 = 995 + 997 + 995 
F81_b3 # 2987 
F81_b10 = 995 + 988 + 988 
F81_b10 # 2971 
F81_u = 994 + 997 + 991 
F81_u #  2982
F81_cu = 987 + 987 + 995 
F81_cu #  2969 
F81_cud = 987 + 993 + 999
F81_cud # 2979 
F81_NJ = 866 + 676 + 667
F81_NJ #   2209 
F81_ML = 996 + 994 + 998
F81_ML #  2988
F81_b1/30 # 99.6%
F81_b2/30 # 99.43%
F81_b3/30 # 99.57%
F81_b10/30 # 99.03%
F81_u/30 # 99.4%
F81_cu/30 # 98.97%
F81_cud/30 # 99.3%
F81_NJ/30 # 73.63%
F81_ML/30 # 99.6% 
# Binomial test F81 

prop.test(c(F81_b1,F81_b2),c(3000,3000)) # no diff, p-value = 0.4565
prop.test(c(F81_b1,F81_b3),c(3000,3000)) # no diff, p-value = 1
prop.test(c(F81_b1,F81_b10),c(3000,3000)) # sign diff, p-value = 0.01216
prop.test(c(F81_b1,F81_u),c(3000,3000)) # no diff, p-value = 0.3601
prop.test(c(F81_b1,F81_cu),c(3000,3000)) # sign diff, p-value = 0.005872
prop.test(c(F81_b1,F81_cud),c(3000,3000)) # no diff, p-value = 0.1626
prop.test(c(F81_b1,F81_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F81_b1,F81_ML),c(3000,3000)) # no diff, p-value = 1

prop.test(c(F81_b2,F81_b3),c(3000,3000)) # no diff,  p-value = 0.5829
prop.test(c(F81_b2,F81_b10),c(3000,3000)) # no diff, p-value = 0.1035
prop.test(c(F81_b2,F81_u),c(3000,3000)) # no diff,  p-value = 1
prop.test(c(F81_b2,F81_cu),c(3000,3000)) # almost sign diff, p-value = p-value = 0.05957
prop.test(c(F81_b2,F81_cud),c(3000,3000)) # no diff,  p-value = 0.6254
prop.test(c(F81_b2,F81_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F81_b2,F81_ML),c(3000,3000)) # no diff, p-value = 0.4565

prop.test(c(F81_b3,F81_b10),c(3000,3000)) # sign diff, p-value = 0.0202
prop.test(c(F81_b3,F81_u),c(3000,3000)) # no diff, p-value = 0.4713
prop.test(c(F81_b3,F81_cu),c(3000,3000)) # sign diff, p-value = 0.0101
prop.test(c(F81_b3,F81_cud),c(3000,3000)) # no diff, p-value = 0.2286
prop.test(c(F81_b3,F81_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F81_b3,F81_ML),c(3000,3000)) # no diff, p-value = 1

prop.test(c(F81_b10,F81_u),c(3000,3000)) # no diff, p-value = 0.1431
prop.test(c(F81_b10,F81_cu),c(3000,3000)) # no diff, p-value = 0.8968
prop.test(c(F81_b10,F81_cud),c(3000,3000)) # no diff, p-value = 0.3202
prop.test(c(F81_b10,F81_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F81_b10,F81_ML),c(3000,3000)) # sign diff,  p-value = 0.01216

prop.test(c(F81_u,F81_cu),c(3000,3000)) # no diff, p-value = 0.08519
prop.test(c(F81_u,F81_cud),c(3000,3000)) # no diff, p-value = 0.748
prop.test(c(F81_u,F81_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F81_u,F81_ML),c(3000,3000)) # no diff, p-value = 0.3601

prop.test(c(F81_cu,F81_cud),c(3000,3000)) # no diff, p-value = 0.21
prop.test(c(F81_cu,F81_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F81_cu,F81_ML),c(3000,3000)) # sign diff,  p-value = 0.005872

prop.test(c(F81_cud,F81_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F81_cud,F81_ML),c(3000,3000)) # no diff,  p-value = 0.1626

prop.test(c(F81_NJ,F81_ML),c(3000,3000)) # sign diff p-value < 2.2e-16


# Number of succ. F84 
F84_b1 =  996 + 994 + 990 
F84_b1 # 2980
F84_b2 =  997 + 994 + 991 
F84_b2 #  2982      
F84_b3 =  995 + 996 + 991 
F84_b3 #  2982 
F84_b10 =  998 + 994 + 991 
F84_b10 #  2983 
F84_u = 995 + 992 + 990 
F84_u #  2977
F84_cu = 998 + 992 + 987 
F84_cu #   2977
F84_cud = 997 + 998 + 992 
F84_cud #  2987
F84_NJ = 918 + 855 + 846 
F84_NJ #   2619 
F84_ML = 993 + 997 + 994 
F84_ML #  2984 
F84_b1/30 # 99.33%
F84_b2/30 # 99.4%
F84_b3/30 # 99.4%
F84_b10/30 # 99.43%
F84_u/30 # 99.23%
F84_cu/30 # 99.23%
F84_cud/30 # 99.57%
F84_NJ/30 # 87.3%
F84_ML/30 # 99.47%
# Binomial test F84 

prop.test(c(F84_b1,F84_b2),c(3000,3000)) # no diff, p-value = 0.8707
prop.test(c(F84_b1,F84_b3),c(3000,3000)) # no diff, p-value = 0.8707
prop.test(c(F84_b1,F84_b10),c(3000,3000)) # no diff, p-value = 0.7415
prop.test(c(F84_b1,F84_u),c(3000,3000)) # no diff, p-value = 0.7595
prop.test(c(F84_b1,F84_cu),c(3000,3000)) # no diff, p-value = 0.7595
prop.test(c(F84_b1,F84_cud),c(3000,3000)) # no diff, p-value = 0.2949
prop.test(c(F84_b1,F84_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F84_b1,F84_ML),c(3000,3000)) # no diff, p-value = 0.616

prop.test(c(F84_b2,F84_b3),c(3000,3000)) # no diff,  p-value = 1
prop.test(c(F84_b2,F84_b10),c(3000,3000)) # no diff, p-value = 1
prop.test(c(F84_b2,F84_u),c(3000,3000)) # no diff,  p-value = 0.5308
prop.test(c(F84_b2,F84_cu),c(3000,3000)) # no diff, p-value = 0.5308
prop.test(c(F84_b2,F84_cud),c(3000,3000)) # no diff,  p-value = 0.4713
prop.test(c(F84_b2,F84_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F84_b2,F84_ML),c(3000,3000)) # no diff, p-value = 0.8634

prop.test(c(F84_b3,F84_b10),c(3000,3000)) # no diff, p-value = 1 
prop.test(c(F84_b3,F84_u),c(3000,3000)) # no diff, p-value = 0.5308
prop.test(c(F84_b3,F84_cu),c(3000,3000)) # no diff, p-value = 0.5308
prop.test(c(F84_b3,F84_cud),c(3000,3000)) # no diff,  p-value = 0.4713
prop.test(c(F84_b3,F84_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F84_b3,F84_ML),c(3000,3000)) # no diff, p-value = 0.8634

prop.test(c(F84_b10,F84_u),c(3000,3000)) # no diff, p-value = 0.4277
prop.test(c(F84_b10,F84_cu),c(3000,3000)) # no diff, p-value = 0.4277
prop.test(c(F84_b10,F84_cud),c(3000,3000)) # no diff, p-value = 0.5829
prop.test(c(F84_b10,F84_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F84_b10,F84_ML),c(3000,3000)) # no diff,  p-value = 1 

prop.test(c(F84_u,F84_cu),c(3000,3000)) # no diff, p-value = 1
prop.test(c(F84_u,F84_cud),c(3000,3000)) # no diff, p-value = 0.1324
prop.test(c(F84_u,F84_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F84_u,F84_ML),c(3000,3000)) # no diff, p-value = 0.3351

prop.test(c(F84_cu,F84_cud),c(3000,3000)) # no diff, p-value = 0.1324
prop.test(c(F84_cu,F84_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F84_cu,F84_ML),c(3000,3000)) # no diff,  p-value = 0.3351

prop.test(c(F84_cud,F84_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(F84_cud,F84_ML),c(3000,3000)) # no diff,  p-value = 0.7097

prop.test(c(F84_NJ,F84_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# Number of suc. GTR 

GTR_b1 = 991 + 995 + 992 
GTR_b1 # 2978
GTR_b2 = 985 + 995 + 997 
GTR_b2 # 2977
GTR_b3 = 992 + 997 + 990
GTR_b3 # 2979
GTR_b10 = 992 + 994 + 996 
GTR_b10 # 2982
GTR_u = 990 + 997 + 992 
GTR_u # 2979
GTR_cu = 993 + 994 + 996 
GTR_cu # 2983
GTR_cud = 988 + 992 + 998
GTR_cud # 2978
GTR_NJ = 887 + 754 + 755 
GTR_NJ # 2396
GTR_ML =991 + 997 + 997 
GTR_ML #   2985 

GTR_b1/30 # 99.27%
GTR_b2/30 # 99.23%
GTR_b3/30 # 99.3%
GTR_b10/30 # 99.4%
GTR_u/30 # 99.3%
GTR_cu/30 # 99.43%
GTR_cud/30 # 99.27%
GTR_NJ/30 # 79.87%
GTR_ML/30 # 99.5%

# Binomial test GTR 

prop.test(c(GTR_b1,GTR_b2),c(3000,3000)) # no diff  p-value = 1
prop.test(c(GTR_b1,GTR_b3),c(3000,3000)) # no diff p-value = 1
prop.test(c(GTR_b1,GTR_b10),c(3000,3000)) # no diff p-value = 0.6341
prop.test(c(GTR_b1,GTR_u),c(3000,3000)) # no diff  p-value = 1
prop.test(c(GTR_b1,GTR_cu),c(3000,3000)) # no diff p-value = 0.5205
prop.test(c(GTR_b1,GTR_cud),c(3000,3000)) # no diff p-value = 1
prop.test(c(GTR_b1,GTR_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(GTR_b1,GTR_ML),c(3000,3000)) # no diff p-value = 0.3224

prop.test(c(GTR_b2,GTR_b3),c(3000,3000)) # no diff p-value = 0.8797
prop.test(c(GTR_b2,GTR_b10),c(3000,3000)) # no diff p-value = 0.5308
prop.test(c(GTR_b2,GTR_u),c(3000,3000)) # no diff p-value = 0.8797
prop.test(c(GTR_b2,GTR_cu),c(3000,3000)) # no diff p-value = 0.4277
prop.test(c(GTR_b2,GTR_cud),c(3000,3000)) # no diff p-value = 1
prop.test(c(GTR_b2,GTR_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(GTR_b2,GTR_ML),c(3000,3000)) # no diff p-value = 0.2546

prop.test(c(GTR_b3,GTR_b10),c(3000,3000)) # no diff p-value = 0.748
prop.test(c(GTR_b3,GTR_u),c(3000,3000)) # no diff p-value = 1
prop.test(c(GTR_b3,GTR_cu),c(3000,3000)) # no diff p-value = 0.6254
prop.test(c(GTR_b3,GTR_cud),c(3000,3000)) # no diff p-value = 1
prop.test(c(GTR_b3,GTR_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(GTR_b3,GTR_ML),c(3000,3000)) # no diff  p-value = 0.4032

prop.test(c(GTR_b10,GTR_u),c(3000,3000)) # no diff p-value = 0.748
prop.test(c(GTR_b10,GTR_cu),c(3000,3000)) # no diff p-value = 1
prop.test(c(GTR_b10,GTR_cud),c(3000,3000)) # no diff p-value = 0.6341
prop.test(c(GTR_b10,GTR_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(GTR_b10,GTR_ML),c(3000,3000)) # no diff p-value = 0.727

prop.test(c(GTR_u,GTR_cu),c(3000,3000)) # no diff p-value = 0.6254
prop.test(c(GTR_u,GTR_cud),c(3000,3000)) # no diff p-value = 1
prop.test(c(GTR_u,GTR_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(GTR_u,GTR_ML),c(3000,3000)) # no diff p-value = 0.4032

prop.test(c(GTR_cu,GTR_cud),c(3000,3000)) # no diff p-value = 0.5205
prop.test(c(GTR_cu,GTR_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(GTR_cu,GTR_ML),c(3000,3000)) # no diff p-value = 0.8593

prop.test(c(GTR_cud,GTR_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(GTR_cud,GTR_ML),c(3000,3000)) # no diff p-value = 0.3224

prop.test(c(GTR_NJ,GTR_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# Number of suc. HKY
HKY_b1 = 991 + 993 + 995 
HKY_b1 # 2979
HKY_b2 = 994 + 996 + 992 
HKY_b2 # 2982
HKY_b3 = 995 + 995 + 993 
HKY_b3 # 2983
HKY_b10 = 991 + 997 + 991 
HKY_b10 # 2979
HKY_u = 989 + 990 + 994 
HKY_u # 2973
HKY_cu = 995 + 998 + 992 
HKY_cu # 2985
HKY_cud = 994 + 992 + 994 
HKY_cud # 2980
HKY_NJ = 924 + 864 + 860 
HKY_NJ # 2648
HKY_ML = 994 + 998 + 993
HKY_ML # 2985

HKY_b1/30 # 99.3% 
HKY_b2/30 # 99.4%
HKY_b3/30 # 99.43%
HKY_b10/30 # 99.3%
HKY_u/30 # 99.1%
HKY_cu/30 # 99.5%
HKY_cud/30 # 99.33%
HKY_NJ/30 # 88.27%
HKY_ML/30 # 99.5%

# Binomial test HKY

prop.test(c(HKY_b1,HKY_b2),c(3000,3000)) # no diff p-value = 0.748
prop.test(c(HKY_b1,HKY_b3),c(3000,3000)) #no diff p-value = 0.6254
prop.test(c(HKY_b1,HKY_b10),c(3000,3000)) # no diff  p-value = 1
prop.test(c(HKY_b1,HKY_u),c(3000,3000)) # no diff p-value = 0.4687
prop.test(c(HKY_b1,HKY_cu),c(3000,3000)) # no diff p-value = 0.4032
prop.test(c(HKY_b1,HKY_cud),c(3000,3000)) # no diff p-value = 1
prop.test(c(HKY_b1,HKY_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(HKY_b1,HKY_ML),c(3000,3000)) # no diff p-value = 0.4032

prop.test(c(HKY_b2,HKY_b3),c(3000,3000)) # no diff p-value = 1
prop.test(c(HKY_b2,HKY_b10),c(3000,3000)) # no diff p-value = 0.748
prop.test(c(HKY_b2,HKY_u),c(3000,3000)) # no diff p-value = 0.2313
prop.test(c(HKY_b2,HKY_cu),c(3000,3000)) # no diff p-value = 0.727
prop.test(c(HKY_b2,HKY_cud),c(3000,3000)) # no diff p-value = 0.8707
prop.test(c(HKY_b2,HKY_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(HKY_b2,HKY_ML),c(3000,3000)) # no diff  p-value = 0.727

prop.test(c(HKY_b3,HKY_b10),c(3000,3000)) # no diff p-value = 0.6254
prop.test(c(HKY_b3,HKY_u),c(3000,3000)) # no diff p-value = 0.1733
prop.test(c(HKY_b3,HKY_cu),c(3000,3000)) # no diff p-value = 0.8593
prop.test(c(HKY_b3,HKY_cud),c(3000,3000)) # no diff p-value = 0.7415
prop.test(c(HKY_b3,HKY_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(HKY_b3,HKY_ML),c(3000,3000)) # no diff p-value = 0.8593

prop.test(c(HKY_b10,HKY_u),c(3000,3000)) # no diff p-value = 0.4687
prop.test(c(HKY_b10,HKY_cu),c(3000,3000)) # no diff p-value = 0.4032
prop.test(c(HKY_b10,HKY_cud),c(3000,3000)) # no diff  p-value = 1
prop.test(c(HKY_b10,HKY_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(HKY_b10,HKY_ML),c(3000,3000)) # no diff p-value = 0.4032
 
prop.test(c(HKY_u,HKY_cu),c(3000,3000)) # almost diff p-value = 0.08851
prop.test(c(HKY_u,HKY_cud),c(3000,3000)) # no diff p-value = 0.3796
prop.test(c(HKY_u,HKY_NJ),c(3000,3000)) # sign diff  p-value < 2.2e-16
prop.test(c(HKY_u,HKY_ML),c(3000,3000)) # no diff p-value = 0.08851

prop.test(c(HKY_cu,HKY_cud),c(3000,3000)) # no diff p-value = 0.4977
prop.test(c(HKY_cu,HKY_NJ),c(3000,3000)) # sign diff  p-value < 2.2e-16
prop.test(c(HKY_cu,HKY_ML),c(3000,3000)) # no diff  p-value = 1

prop.test(c(HKY_cud,HKY_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(HKY_cud,HKY_ML),c(3000,3000)) # no diff p-value = 0.4977

prop.test(c(HKY_NJ,HKY_ML),c(3000,3000)) # sign diff p-value < 2.2e-16
