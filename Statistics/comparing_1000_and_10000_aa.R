# comparing 1000 and 10000 trained aa NN 
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

# JTT model 1000 & 10000 binomial test 

prop.test(c(JTT_b1 ,JTT_mod_10000_b1),c(3000,3000)) # sign diff p-value = 8.423e-09
prop.test(c(JTT_b2 ,JTT_mod_10000_b2),c(3000,3000)) # sign diff p-value = 4.535e-09
prop.test(c(JTT_b3 ,JTT_mod_10000_b3),c(3000,3000)) # sign diff p-value = 1.945e-09
prop.test(c(JTT_b10 ,JTT_mod_10000_b10),c(3000,3000)) # sign diff p-value = 4.492e-13
prop.test(c(JTT_u ,JTT_mod_10000_u),c(3000,3000)) # sign diff p-value = 3.132e-13
prop.test(c(JTT_cu ,JTT_mod_10000_cu),c(3000,3000)) # sign diff p-value = 3.073e-14
prop.test(c(JTT_cud ,JTT_mod_10000_cud),c(3000,3000)) # sign diff p-value = 2.061e-05

# Diff in accuracy % 
JTT_b1_diff = (JTT_mod_10000_b1/3000 - JTT_b1/3000)*100  
JTT_b1_diff #   3.3%
JTT_b2_diff = (JTT_mod_10000_b2/3000 - JTT_b2/3000)*100  
JTT_b2_diff # 3%
JTT_b3_diff = (JTT_mod_10000_b3/3000 - JTT_b3/3000)*100  
JTT_b3_diff # 2.733%
JTT_b10_diff = (JTT_mod_10000_b10/3000 - JTT_b10/3000)*100 
JTT_b10_diff # 5.866%
JTT_u_diff = (JTT_mod_10000_u/3000 - JTT_u/3000)*100 
JTT_u_diff # 3.9%
JTT_cu_diff = (JTT_mod_10000_cu/3000 - JTT_cu/3000)*100 
JTT_cu_diff # 4.266%
JTT_cud_diff = (JTT_mod_10000_cud/3000 - JTT_cud/3000)*100  
JTT_cud_diff # 2.4%

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

# LG model 1000 & 10000 binomial test 

prop.test(c(LG_b1 ,LG_mod_10000_b1),c(3000,3000)) # sign diff p-value = 2.072e-06
prop.test(c(LG_b2 ,LG_mod_10000_b2),c(3000,3000)) # sign diff p-value = 1.116e-11
prop.test(c(LG_b3 ,LG_mod_10000_b3),c(3000,3000)) # sign diff p-value = 1.457e-13
prop.test(c(LG_b10 ,LG_mod_10000_b10),c(3000,3000)) # sign diff p-value = 2.761e-08
prop.test(c(LG_u ,LG_mod_10000_u),c(3000,3000)) # sign diff p-value = 1.214e-14
prop.test(c(LG_cu ,LG_mod_10000_cu),c(3000,3000)) # sign diff p-value = 1.968e-07
prop.test(c(LG_cud ,LG_mod_10000_cud),c(3000,3000)) # sign diff  p-value = 2.811e-10

# Diff in accuracy % 
LG_b1_diff = (LG_mod_10000_b1/3000 - LG_b1/3000)*100  
LG_b1_diff #   2.7%
LG_b2_diff = (LG_mod_10000_b2/3000 - LG_b2/3000)*100  
LG_b2_diff # 3.7%
LG_b3_diff = (LG_mod_10000_b3/3000 - LG_b3/3000)*100  
LG_b3_diff # 3.5%
LG_b10_diff = (LG_mod_10000_b10/3000 - LG_b10/3000)*100 
LG_b10_diff # 4.4%
LG_u_diff = (LG_mod_10000_u/3000 - LG_u/3000)*100 
LG_u_diff # 4.13%
LG_cu_diff = (LG_mod_10000_cu/3000 - LG_cu/3000)*100 
LG_cu_diff # 2.9%
LG_cud_diff = (LG_mod_10000_cud/3000 - LG_cud/3000)*100  
LG_cud_diff # 3.56%

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

# WAG model 1000 & 10000 binomial test 

prop.test(c(WAG_b1 ,WAG_mod_10000_b1),c(3000,3000)) # sign diff p-value = 1.875e-08
prop.test(c(WAG_b2 ,WAG_mod_10000_b2),c(3000,3000)) # sign diff p-value = 1.285e-13
prop.test(c(WAG_b3 ,WAG_mod_10000_b3),c(3000,3000)) # sign diff p-value = 3.08e-16
prop.test(c(WAG_b10 ,WAG_mod_10000_b10),c(3000,3000)) # sign diff p-value = 1.572e-11
prop.test(c(WAG_u ,WAG_mod_10000_u),c(3000,3000)) # sign diff p-value = 5.331e-16
prop.test(c(WAG_cu ,WAG_mod_10000_cu),c(3000,3000)) # sign diff p-value = 1.369e-09
prop.test(c(WAG_cud ,WAG_mod_10000_cud),c(3000,3000)) # sign diff p-value = 2.006e-07

# Diff in accuracy % 
WAG_b1_diff = (WAG_mod_10000_b1/3000 - WAG_b1/3000)*100  
WAG_b1_diff #   3.26%
WAG_b2_diff = (WAG_mod_10000_b2/3000 - WAG_b2/3000)*100  
WAG_b2_diff # 3.833%
WAG_b3_diff = (WAG_mod_10000_b3/3000 - WAG_b3/3000)*100  
WAG_b3_diff # 4.06%
WAG_b10_diff = (WAG_mod_10000_b10/3000 - WAG_b10/3000)*100 
WAG_b10_diff # 5.6%
WAG_u_diff = (WAG_mod_10000_u/3000 - WAG_u/3000)*100 
WAG_u_diff # 4.7%
WAG_cu_diff = (WAG_mod_10000_cu/3000 - WAG_cu/3000)*100 
WAG_cu_diff # 3.6%
WAG_cud_diff = (WAG_mod_10000_cud/3000 - WAG_cud/3000)*100  
WAG_cud_diff # 3.16%


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

# DAY model 1000 & 10000 binomial test 

prop.test(c(DAY_b1 ,DAY_mod_10000_b1),c(3000,3000)) # sign diff p-value = 6.335e-05
prop.test(c(DAY_b2 ,DAY_mod_10000_b2),c(3000,3000)) # sign diff p-value = 1.469e-07
prop.test(c(DAY_b3 ,DAY_mod_10000_b3),c(3000,3000)) # sign diff p-value = 2.377e-11
prop.test(c(DAY_b10 ,DAY_mod_10000_b10),c(3000,3000)) # sign diff p-value = 9.12e-10
prop.test(c(DAY_u ,DAY_mod_10000_u),c(3000,3000)) # sign diff p-value = 3.665e-11
prop.test(c(DAY_cu ,DAY_mod_10000_cu),c(3000,3000)) # sign diff p-value = 4.554e-08
prop.test(c(DAY_cud ,DAY_mod_10000_cud),c(3000,3000)) # sign diff p-value = 3.849e-08

# Diff in accuracy % 
DAY_b1_diff = (DAY_mod_10000_b1/3000 - DAY_b1/3000)*100  
DAY_b1_diff # 2.33%
DAY_b2_diff = (DAY_mod_10000_b2/3000 - DAY_b2/3000)*100  
DAY_b2_diff # 
DAY_b3_diff = (DAY_mod_10000_b3/3000 - DAY_b3/3000)*100  
DAY_b3_diff # 2.86%
DAY_b10_diff = (DAY_mod_10000_b10/3000 - DAY_b10/3000)*100  
DAY_b10_diff # 5.1%
DAY_u_diff = (DAY_mod_10000_u/3000 - DAY_u/3000)*100  
DAY_u_diff # 3.73%
DAY_cu_diff = (DAY_mod_10000_cu/3000 - DAY_cud/3000)*100  
DAY_cu_diff # 3.26%
DAY_cud_diff = (DAY_mod_10000_cud/3000 - DAY_cud/3000)*100  
DAY_cud_diff # 3.33%
