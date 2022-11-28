# Number of suc. JTT_sb 
JTT_b1_sb = 438 + 338 + 351 
JTT_b1_sb  # 1127
JTT_b2_sb = 382 + 401 + 384 
JTT_b2_sb # 1167
JTT_b3_sb = 449 + 371 + 350 
JTT_b3_sb # 1170
JTT_b10_sb = 383 + 400 + 294 
JTT_b10_sb # 1077
JTT_u_sb = 400 + 347 + 395 
JTT_u_sb # 1142
JTT_cu_sb = 366 + 373 + 359 
JTT_cu_sb # 1098
JTT_cud_sb = 427 + 372 + 338 
JTT_cud_sb # 1137
JTT_NJ_sb = 495 + 464 + 428 
JTT_NJ_sb # 1387
JTT_ML_sb = 605 + 596 + 605 
JTT_ML_sb # 1806  

JTT_b1_sb/30 # 37.57%
JTT_b2_sb/30 # 38.9%
JTT_b3_sb/30 # 39.0%
JTT_b10_sb/30 # 35.9%
JTT_u_sb/30 # 38.07%
JTT_cu_sb/30 # 36.6%
JTT_cud_sb/30 # 37.9%
JTT_NJ_sb/30 # 37.9%
JTT_ML_sb/30 # 46.23%

# JTT_sb binomial test 

prop.test(c(JTT_b1_sb,JTT_b2_sb),c(3000,3000)) # no diff -value = 0.3002
prop.test(c(JTT_b1_sb,JTT_b3_sb),c(3000,3000)) # no diff  p-value = 0.2646
prop.test(c(JTT_b1_sb,JTT_b10_sb),c(3000,3000)) # no diff p-value = 0.1894
prop.test(c(JTT_b1_sb,JTT_u_sb),c(3000,3000)) # no diff p-value = 0.7094
prop.test(c(JTT_b1_sb,JTT_cu_sb),c(3000,3000)) # no diff p-value = 0.4542
prop.test(c(JTT_b1_sb,JTT_cud_sb),c(3000,3000)) # no diff p-value = 0.8106
prop.test(c(JTT_b1_sb,JTT_NJ_sb),c(3000,3000)) # sign diff p-value = 1.228e-11
prop.test(c(JTT_b1_sb,JTT_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_b2_sb,JTT_b3_sb),c(3000,3000)) # no diff p-value = 0.9578
prop.test(c(JTT_b2_sb,JTT_b10_sb),c(3000,3000)) # sign diff p-value = 0.01757
prop.test(c(JTT_b2_sb,JTT_u_sb),c(3000,3000)) # no diff p-value = 0.5243
prop.test(c(JTT_b2_sb,JTT_cu_sb),c(3000,3000)) # almost diff p-value = 0.07015
prop.test(c(JTT_b2_sb,JTT_cud_sb),c(3000,3000)) # no diff p-value = 0.4414
prop.test(c(JTT_b2_sb,JTT_NJ_sb),c(3000,3000)) # sign diff  p-value = 1.077e-08
prop.test(c(JTT_b2_sb,JTT_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_b3_sb,JTT_b10_sb),c(3000,3000)) # sign diff p-value = 0.01413
prop.test(c(JTT_b3_sb,JTT_u_sb),c(3000,3000)) # no diff p-value = 0.4739
prop.test(c(JTT_b3_sb,JTT_cu_sb),c(3000,3000)) # almost diff  p-value = 0.05871
prop.test(c(JTT_b3_sb,JTT_cud_sb),c(3000,3000)) # no diff p-value = 0.3958
prop.test(c(JTT_b3_sb,JTT_NJ_sb),c(3000,3000)) # sign diff p-value = 1.711e-08
prop.test(c(JTT_b3_sb,JTT_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_b10_sb,JTT_u_sb),c(3000,3000)) # almost diff p-value = 0.08699
prop.test(c(JTT_b10_sb,JTT_cu_sb),c(3000,3000)) # no diff p-value = 0.5912 
prop.test(c(JTT_b10_sb,JTT_cud_sb),c(3000,3000)) # no diff p-value = 0.1144
prop.test(c(JTT_b10_sb,JTT_NJ_sb),c(3000,3000)) # sign diff p-value = 5.111e-16
prop.test(c(JTT_b10_sb,JTT_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_u_sb,JTT_cu_sb),c(3000,3000)) # no diff p-value = 0.2511
prop.test(c(JTT_u_sb,JTT_cud_sb),c(3000,3000)) # no diff p-value = 0.9153
prop.test(c(JTT_u_sb,JTT_NJ_sb),c(3000,3000)) # sign diff p-value = 1.781e-10
prop.test(c(JTT_u_sb,JTT_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_cu_sb,JTT_cud_sb),c(3000,3000)) # no diff p-value = 0.3102
prop.test(c(JTT_cu_sb,JTT_NJ_sb),c(3000,3000)) # sign diff p-value = 4.414e-14
prop.test(c(JTT_cu_sb,JTT_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_cud_sb,JTT_NJ_sb),c(3000,3000)) # sign diff p-value = 7.433e-11
prop.test(c(JTT_cud_sb,JTT_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_NJ_sb,JTT_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

# Number of suc. LG_sb

LG_b1_sb = 380 + 397 + 362 
LG_b1_sb # 1139
LG_b2_sb = 368 + 392 + 388 
LG_b2_sb # 1148
LG_b3_sb = 394 + 423 + 345 
LG_b3_sb # 1162
LG_b10_sb = 399 + 336 + 356 
LG_b10_sb # 1091
LG_u_sb = 330 + 398 + 417 
LG_u_sb # 1145
LG_cu_sb = 379 + 390 + 376 
LG_cu_sb # 1145
LG_cud_sb = 383 + 380 + 402 
LG_cud_sb # 1165
LG_NJ_sb = 499 + 465 + 459 
LG_NJ_sb # 1423
LG_ML_sb = 599 + 597 + 609 
LG_ML_sb #  1805

LG_b1_sb/30 # 37.97%
LG_b2_sb/30 # 38.27%
LG_b3_sb/30 # 38.73%
LG_b10_sb/30 # 36.37%
LG_u_sb/30 # 38.17%
LG_cu_sb/30 # 38.17%
LG_cud_sb/30 # 38.83%
LG_NJ_sb/30 # 47.43%
LG_ML_sb/30 # 60.17%

# LG_sb binomial test 

prop.test(c(LG_b1_sb,LG_b2_sb),c(3000,3000)) # no diff  p-value = 0.8316
prop.test(c(LG_b1_sb,LG_b3_sb),c(3000,3000)) # no diff p-value = 0.5591
prop.test(c(LG_b1_sb,LG_b10_sb),c(3000,3000)) # no diff p-value = 0.2093
prop.test(c(LG_b1_sb,LG_u_sb),c(3000,3000)) # no diff p-value = 0.8942
prop.test(c(LG_b1_sb,LG_cu_sb),c(3000,3000)) # no diff p-value = 0.8942
prop.test(c(LG_b1_sb,LG_cud_sb),c(3000,3000)) # no diff  p-value = 0.5069
prop.test(c(LG_b1_sb,LG_NJ_sb),c(3000,3000)) # sign diff p-value = 1.511e-13
prop.test(c(LG_b1_sb,LG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_b2_sb,LG_b3_sb),c(3000,3000)) # no diff p-value = 0.7302
prop.test(c(LG_b2_sb,LG_b10_sb),c(3000,3000)) # no diff p-value = 0.135
prop.test(c(LG_b2_sb,LG_u_sb),c(3000,3000)) # no diff p-value = 0.9576
prop.test(c(LG_b2_sb,LG_cu_sb),c(3000,3000)) # no diff p-value = 0.9576
prop.test(c(LG_b2_sb,LG_cud_sb),c(3000,3000)) # no diff p-value = 0.6713
prop.test(c(LG_b2_sb,LG_NJ_sb),c(3000,3000)) # sign diff p-value = 8.798e-13
prop.test(c(LG_b2_sb,LG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_b3_sb,LG_b10_sb),c(3000,3000)) # almost diff p-value = 0.06202
prop.test(c(LG_b3_sb,LG_u_sb),c(3000,3000)) # no diff p-value = 0.6711
prop.test(c(LG_b3_sb,LG_cu_sb),c(3000,3000)) # no diff p-value = 0.6711
prop.test(c(LG_b3_sb,LG_cud_sb),c(3000,3000)) # no diff p-value = 0.9577 
prop.test(c(LG_b3_sb,LG_NJ_sb),c(3000,3000)) # sign diff p-value = 1.216e-11
prop.test(c(LG_b3_sb,LG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_b10_sb,LG_u_sb),c(3000,3000)) # no diff p-value = 0.157
prop.test(c(LG_b10_sb,LG_cu_sb),c(3000,3000)) #no diff p-value = 0.157
prop.test(c(LG_b10_sb,LG_cud_sb),c(3000,3000)) # almost diff p-value = 0.0517
prop.test(c(LG_b10_sb,LG_NJ_sb),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_b10_sb,LG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_u_sb,LG_cu_sb),c(3000,3000)) # no diff  p-value = 1
prop.test(c(LG_u_sb,LG_cud_sb),c(3000,3000)) # no diff p-value = 0.6142
prop.test(c(LG_u_sb,LG_NJ_sb),c(3000,3000)) # sign diff p-value = 4.922e-13
prop.test(c(LG_u_sb,LG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_cu_sb,LG_cud_sb),c(3000,3000)) # no diff p-value = 0.6142
prop.test(c(LG_cu_sb,LG_NJ_sb),c(3000,3000)) # sign diff p-value = 4.922e-13
prop.test(c(LG_cu_sb,LG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_cud_sb,LG_NJ_sb),c(3000,3000)) # sign diff  p-value = 2.096e-11
prop.test(c(LG_cud_sb,LG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_NJ_sb,LG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

# Number of suc. WAG_sb
WAG_b1_sb = 350 + 430 + 368
WAG_b1_sb # 1148 
WAG_b2_sb = 388 + 401 + 398 
WAG_b2_sb # 1187
WAG_b3_sb = 416 + 372 + 413 
WAG_b3_sb # 1201
WAG_b10_sb = 380 + 388 + 356 
WAG_b10_sb # 1124
WAG_u_sb = 430 + 342 + 394 
WAG_u_sb # 1166
WAG_cu_sb = 417 + 406 + 356 
WAG_cu_sb # 1179
WAG_cud_sb = 402 + 353 + 417 
WAG_cud_sb # 1172
WAG_NJ_sb = 511 + 459 + 432
WAG_NJ_sb # 1402
WAG_ML_sb = 593 + 611 + 598 
WAG_ML_sb #   1802

WAG_b1_sb/30 # 38.27%
WAG_b2_sb/30 # 39.57%
WAG_b3_sb/30 # 40.03%
WAG_b10_sb/30 # 37.47%
WAG_u_sb/30 # 38.87%
WAG_cu_sb/30 # 39.3%
WAG_cud_sb/30 # 39.07%
WAG_NJ_sb/30 # 46.73%
WAG_ML_sb/30 # 60.07%

# WAG_sb binomial test 

prop.test(c(WAG_b1_sb,WAG_b2_sb),c(3000,3000)) # no diff p-value = 0.3143
prop.test(c(WAG_b1_sb,WAG_b3_sb),c(3000,3000)) # no diff p-value = 0.169
prop.test(c(WAG_b1_sb,WAG_b10_sb),c(3000,3000)) # no diff p-value = 0.5404
prop.test(c(WAG_b1_sb,WAG_u_sb),c(3000,3000)) # no diff p-value = 0.6521
prop.test(c(WAG_b1_sb,WAG_cu_sb),c(3000,3000)) # no diff p-value = 0.4267
prop.test(c(WAG_b1_sb,WAG_cud_sb),c(3000,3000)) # no diff p-value = 0.542
prop.test(c(WAG_b1_sb,WAG_NJ_sb),c(3000,3000)) # sign diff p-value = 3.917e-11
prop.test(c(WAG_b1_sb,WAG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_b2_sb,WAG_b3_sb),c(3000,3000)) # no diff p-value = 0.7317
prop.test(c(WAG_b2_sb,WAG_b10_sb),c(3000,3000)) # no diff p-value = 0.1
prop.test(c(WAG_b2_sb,WAG_u_sb),c(3000,3000)) # no diff p-value = 0.5969
prop.test(c(WAG_b2_sb,WAG_cu_sb),c(3000,3000)) # no diff p-value = 0.8533
prop.test(c(WAG_b2_sb,WAG_cud_sb),c(3000,3000)) # no diff p-value = 0.7114
prop.test(c(WAG_b2_sb,WAG_NJ_sb),c(3000,3000)) # sign diff p-value = 2.432e-08
prop.test(c(WAG_b2_sb,WAG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_b3_sb,WAG_b10_sb),c(3000,3000)) # sign diff p-value = 0.04401
prop.test(c(WAG_b3_sb,WAG_u_sb),c(3000,3000)) # no diff p-value = 0.3691
prop.test(c(WAG_b3_sb,WAG_cu_sb),c(3000,3000)) # no diff p-value = 0.5795
prop.test(c(WAG_b3_sb,WAG_cud_sb),c(3000,3000)) # no diff p-value = 0.4597
prop.test(c(WAG_b3_sb,WAG_NJ_sb),c(3000,3000)) # sign diff p-value = 1.89e-07
prop.test(c(WAG_b3_sb,WAG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_b10_sb,WAG_u_sb),c(3000,3000)) # no diff p-value = 0.2759
prop.test(c(WAG_b10_sb,WAG_cu_sb),c(3000,3000)) # no diff p-value = 0.1517
prop.test(c(WAG_b10_sb,WAG_cud_sb),c(3000,3000)) # no diff p-value = 0.2119
prop.test(c(WAG_b10_sb,WAG_NJ_sb),c(3000,3000)) # sign diff p-value = 4.386e-13
prop.test(c(WAG_b10_sb,WAG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_u_sb,WAG_cu_sb),c(3000,3000)) # no diff p-value = 0.7509
prop.test(c(WAG_u_sb,WAG_cud_sb),c(3000,3000)) # no diff p-value = 0.8947
prop.test(c(WAG_u_sb,WAG_NJ_sb),c(3000,3000)) # sign diff p-value = 8.701e-10
prop.test(c(WAG_u_sb,WAG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_cu_sb,WAG_cud_sb),c(3000,3000)) # no diff p-value = 0.8739
prop.test(c(WAG_cu_sb,WAG_NJ_sb),c(3000,3000)) # sign diff p-value = 7.091e-09
prop.test(c(WAG_cu_sb,WAG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_cud_sb,WAG_NJ_sb),c(3000,3000)) # sign diff p-value = 2.325e-09
prop.test(c(WAG_cud_sb,WAG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_NJ_sb,WAG_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

# Number of suc. DAY_sb 
DAY_b1_sb = 406 + 387 + 409 
DAY_b1_sb # 1202
DAY_b2_sb = 453 + 389 + 363 
DAY_b2_sb # 1205
DAY_b3_sb = 436 + 431 + 337 
DAY_b3_sb # 1204
DAY_b10_sb = 369 + 336 + 378 
DAY_b10_sb # 1083
DAY_u_sb = 464 + 361 + 362 
DAY_u_sb # 1187
DAY_cu_sb = 450 + 366 + 390 
DAY_cu_sb # 1206
DAY_cud_sb = 450 + 320 + 407 
DAY_cud_sb # 1177
DAY_NJ_sb = 512 + 473 + 435 
DAY_NJ_sb # 1420
DAY_ML_sb = 596 + 594 + 608 
DAY_ML_sb #   1798

DAY_b1_sb/30 # 40.07%
DAY_b2_sb/30 # 40.17%
DAY_b3_sb/30 # 40.13% 
DAY_b10_sb/30 # 36.1%
DAY_u_sb/30 # 39.57%
DAY_cu_sb/30 # 40.2%
DAY_cud_sb/30 # 39.23%
DAY_NJ_sb/30 # 47.33%
DAY_ML_sb/30 # 59.93%

# DAY_sb binomial test 

prop.test(c(DAY_b1_sb,DAY_b2_sb),c(3000,3000)) # no diff  p-value = 0.958
prop.test(c(DAY_b1_sb,DAY_b3_sb),c(3000,3000)) # no diff p-value = 0.979
prop.test(c(DAY_b1_sb,DAY_b10_sb),c(3000,3000)) # sign diff p-value = 0.001706
prop.test(c(DAY_b1_sb,DAY_u_sb),c(3000,3000)) # no diff p-value = 0.712
prop.test(c(DAY_b1_sb,DAY_cu_sb),c(3000,3000)) # no diff p-value = 0.937
prop.test(c(DAY_b1_sb,DAY_cud_sb),c(3000,3000)) # no diff p-value = 0.5265
prop.test(c(DAY_b1_sb,DAY_NJ_sb),c(3000,3000)) # sign diff p-value = 1.624e-08
prop.test(c(DAY_b1_sb,DAY_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_b1_sb,DAY_b3_sb),c(3000,3000)) # no diff p-value = 0.979
prop.test(c(DAY_b2_sb,DAY_b10_sb),c(3000,3000)) # sign diff p-value = 0.001299
prop.test(c(DAY_b2_sb,DAY_u_sb),c(3000,3000)) # no diff p-value = 0.654
prop.test(c(DAY_b2_sb,DAY_cu_sb),c(3000,3000)) # no diff p-value = 1
prop.test(c(DAY_b2_sb,DAY_cud_sb),c(3000,3000)) # no diff p-value = 0.4762
prop.test(c(DAY_b2_sb,DAY_NJ_sb),c(3000,3000)) # sign diff p-value = 2.56e-08
prop.test(c(DAY_b2_sb,DAY_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_b3_sb,DAY_b10_sb),c(3000,3000)) # sign diff p-value = 0.001424
prop.test(c(DAY_b3_sb,DAY_u_sb),c(3000,3000)) # no diff p-value = 0.6731
prop.test(c(DAY_b3_sb,DAY_cu_sb),c(3000,3000)) # no diff p-value = 0.979
prop.test(c(DAY_b3_sb,DAY_cud_sb),c(3000,3000)) # no diff p-value = 0.4927
prop.test(c(DAY_b3_sb,DAY_NJ_sb),c(3000,3000)) # sign diff p-value = 2.201e-08
prop.test(c(DAY_b3_sb,DAY_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_b10_sb,DAY_u_sb),c(3000,3000)) # sign diff p-value = 0.006109
prop.test(c(DAY_b10_sb,DAY_cu_sb),c(3000,3000)) # sign diff p-value = 0.001185
prop.test(c(DAY_b10_sb,DAY_cud_sb),c(3000,3000)) # sign diff p-value = 0.01322
prop.test(c(DAY_b10_sb,DAY_NJ_sb),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_b10_sb,DAY_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_u_sb,DAY_cu_sb),c(3000,3000)) # no diff p-value = 0.6351
prop.test(c(DAY_u_sb,DAY_cud_sb),c(3000,3000)) # no diff p-value = 0.812
prop.test(c(DAY_u_sb,DAY_NJ_sb),c(3000,3000)) # sign diff p-value = 1.519e-09
prop.test(c(DAY_u_sb,DAY_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_cu_sb,DAY_cud_sb),c(3000,3000)) # no diff p-value = 0.4601
prop.test(c(DAY_cu_sb,DAY_NJ_sb),c(3000,3000)) # sign diff p-value = 2.975e-08
prop.test(c(DAY_cu_sb,DAY_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_cud_sb,DAY_NJ_sb),c(3000,3000)) # sign diff p-value = 2.871e-10
prop.test(c(DAY_cud_sb,DAY_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_NJ_sb,DAY_ML_sb),c(3000,3000)) # sign diff p-value < 2.2e-16
