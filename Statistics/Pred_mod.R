# Number of suc DNA pred mod 

NN_b1 = 988 + 980 + 828 + 834 + 968 
NN_b1 # 4598
NN_b2 = 984 + 979 + 833 + 851 + 979 
NN_b2 # 4626
NN_b3 = 981 + 973 + 829 + 838 + 958 
NN_b3 # 4579
NN_b10 = 982 + 980 + 845 + 849 + 975 
NN_b10 # 4631
NN_u = 992 + 975 + 755 + 752 + 921 
NN_u # 4395
NN_cu = 980 + 974 + 849 + 859 + 968
NN_cu # 4630
NN_cud = 974 + 980 + 844 + 841 + 971
NN_cud # 4610
ML_mod = 999 + 1000 + 689 + 681 + 985 
ML_mod #   4354

NN_b1/50 # 91.96%
NN_b2/50 # 92.52%
NN_b3/50 # 91.58% 
NN_b10/50 # 92.62%
NN_u/50 # 87.9%
NN_cu/50 # 92.6%
NN_cud/50 # 92.2%
ML_mod/50 # 87.08%

# Binomial test DNA pred mod 

prop.test(c(NN_b1,NN_b2),c(5000,5000)) # no diff p-value = 0.3129
prop.test(c(NN_b1,NN_b3),c(5000,5000)) # no diff p-value = 0.5125
prop.test(c(NN_b1,NN_b10),c(5000,5000)) # no diff p-value = 0.2303
prop.test(c(NN_b1,NN_u),c(5000,5000)) # sign diff p-value = 1.913e-11
prop.test(c(NN_b1,NN_cu),c(5000,5000)) # no diff p-value = 0.2455
prop.test(c(NN_b1,NN_cud),c(5000,5000)) # no diff p-value = 0.6838
prop.test(c(NN_b1,ML_mod),c(5000,5000)) # sign diff p-value = 2.13e-15

prop.test(c(NN_b2,NN_b3),c(5000,5000)) # almost diff p-value = 0.08905
prop.test(c(NN_b2,NN_b10),c(5000,5000)) # no diff p-value = 0.8788
prop.test(c(NN_b2,NN_u),c(5000,5000)) # sign diff p-value = 9.986e-15
prop.test(c(NN_b2,NN_cu),c(5000,5000)) # no diff p-value = 0.909
prop.test(c(NN_b2,NN_cud),c(5000,5000)) # no diff p-value = 0.5723
prop.test(c(NN_b2,ML_mod),c(5000,5000)) # sign diff p-value < 2.2e-16

prop.test(c(NN_b3,NN_b10),c(5000,5000)) # almost diff p-value = 0.05866
prop.test(c(NN_b3,NN_u),c(5000,5000)) # sign diff p-value = 1.63e-09
prop.test(c(NN_b3,NN_cu),c(5000,5000)) # almost diff p-value = 0.06394
prop.test(c(NN_b3,NN_cud),c(5000,5000)) # no diff p-value = 0.2718
prop.test(c(NN_b3,ML_mod),c(5000,5000)) # sign diff p-value = 4.002e-13

prop.test(c(NN_b10,NN_u),c(5000,5000)) # sign diff p-value = 2.268e-15
prop.test(c(NN_b10,NN_cu),c(5000,5000)) # no diff p-value = 1
prop.test(c(NN_b10,NN_cud),c(5000,5000)) # no diff p-value = 0.4501
prop.test(c(NN_b10,ML_mod),c(5000,5000)) # sign diff p-value < 2.2e-16

prop.test(c(NN_u,NN_cu),c(5000,5000)) # sign diff p-value = 3.06e-15
prop.test(c(NN_u,NN_cud),c(5000,5000)) # sign diff p-value = 8.725e-13
prop.test(c(NN_u,ML_mod),c(5000,5000)) # no diff  p-value = 0.2266

prop.test(c(NN_cu,NN_cud),c(5000,5000)) # no diff p-value = 0.4734
prop.test(c(NN_cu,ML_mod),c(5000,5000)) # sign diff p-value < 2.2e-16

prop.test(c(NN_cud,ML_mod),c(5000,5000)) # sign diff p-value < 2.2e-16


# Number of suc Amino acids pred mod 

NN_b1_aa = 988 + 993 + 989 + 997
NN_b1_aa #   3967
NN_b2_aa = 1000 + 993 + 992 + 995
NN_b2_aa # 3980
NN_b3_aa = 984 + 989 + 992 + 994 
NN_b3_aa # 3959
NN_b10_aa = 986 + 992 + 989 + 1000
NN_b10_aa # 3967
NN_u_aa = 995 + 992 + 990 + 996 
NN_u_aa # 3973
NN_cu_aa = 971 + 985 + 976 + 993
NN_cu_aa # 3925
NN_cud_aa = 959 + 997 + 973 + 979 
NN_cud_aa # 3908
ML_aa = 999 + 999 + 999 + 1000
ML_aa # 3997

NN_b1_aa/40 # 99.18%
NN_b2_aa/40 # 99.5%
NN_b3_aa/40 # 98.98%
NN_b10_aa/40 # 99.18%
NN_u_aa/40 # 99.33%
NN_cu_aa/40 # 98.16%
NN_cud_aa/40 # 97.7%
ML_aa/40 # 99.93%

# Binomial test AA pred mod 

prop.test(c(NN_b1_aa,NN_b2_aa),c(4000,4000)) # no diff p-value = 0.09817
prop.test(c(NN_b1_aa,NN_b3_aa),c(4000,4000)) # no diff p-value = 0.4136
prop.test(c(NN_b1_aa,NN_b10_aa),c(4000,4000)) # no diff  p-value = 1
prop.test(c(NN_b1_aa,NN_u_aa),c(4000,4000)) # no diff p-value = 0.517
prop.test(c(NN_b1_aa,NN_cu_aa),c(4000,4000)) # sign diff p-value = 7.123e-05
prop.test(c(NN_b1_aa,NN_cud_aa),c(4000,4000)) # sign diff p-value = 1.707e-07
prop.test(c(NN_b1_aa,ML_aa),c(4000,4000)) # sign diff p-value = 1.271e-06

prop.test(c(NN_b2_aa,NN_b3_aa),c(4000,4000)) # sign diff p-value = 0.01015
prop.test(c(NN_b2_aa,NN_b10_aa),c(4000,4000)) # no diff p-value = 0.09817
prop.test(c(NN_b2_aa,NN_u_aa),c(4000,4000)) # no diff p-value = 0.3801
prop.test(c(NN_b2_aa,NN_cu_aa),c(4000,4000)) # sign diff p-value = 2.497e-08
prop.test(c(NN_b2_aa,NN_cud_aa),c(4000,4000)) # sign diff p-value = 1.415e-11
prop.test(c(NN_b2_aa,ML_aa),c(4000,4000)) # sign diff p-value = 0.0008347

prop.test(c(NN_b3_aa,NN_b10_aa),c(4000,4000)) # no diff p-value = 0.1134
prop.test(c(NN_b3_aa,NN_u_aa),c(4000,4000)) # no diff  p-value = 0.1134
prop.test(c(NN_b3_aa,NN_cu_aa),c(4000,4000)) # sign diff p-value = 0.002026
prop.test(c(NN_b3_aa,NN_cud_aa),c(4000,4000)) # sign diff p-value = 1.231e-05
prop.test(c(NN_b3_aa,ML_aa),c(4000,4000)) # sign diff  p-value = 2.227e-08

prop.test(c(NN_b10_aa,NN_u_aa),c(4000,4000)) # no diff p-value = 0.517
prop.test(c(NN_b10_aa,NN_cu_aa),c(4000,4000)) # sign diff  p-value = 7.123e-05
prop.test(c(NN_b10_aa,NN_cud_aa),c(4000,4000)) # sign diff p-value = 1.707e-07
prop.test(c(NN_b10_aa,ML_aa),c(4000,4000)) # sign diff p-value = 1.271e-06

prop.test(c(NN_u_aa,NN_cu_aa),c(4000,4000)) # sign diff p-value = 2.818e-06
prop.test(c(NN_u_aa,NN_cud_aa),c(4000,4000)) # sign diff p-value = 3.4e-09
prop.test(c(NN_u_aa,ML_aa),c(4000,4000)) # sign diff p-value = 2.587e-05

prop.test(c(NN_cu_aa,NN_cud_aa),c(4000,4000)) # no diff p-value = 0.2108
prop.test(c(NN_cu_aa,ML_aa),c(4000,4000)) # sign diff p-value = 6.549e-16

prop.test(c(NN_cud_aa,ML_aa),c(4000,4000)) # sign diff p-value < 2.2e-16

