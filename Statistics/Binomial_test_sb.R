# Number of succ. JC 
JC_sb_b1 = 481 + 448 + 470
JC_sb_b1 # 1399
JC_sb_b2 = 485 + 437 + 494 
JC_sb_b2 # 1416
JC_sb_b3 = 424 + 477 + 514 
JC_sb_b3 # 1415
JC_sb_b10 = 455 + 520 + 440 
JC_sb_b10 # 1415
JC_sb_u = 471 + 471 + 447 
JC_sb_u # 1389
JC_sb_cu = 454 + 508 + 455 
JC_sb_cu # 1417
JC_sb_cud = 499 + 454 + 474 
JC_sb_cud # 1427
JC_sb_NJ = 556 + 292 + 305 
JC_sb_NJ # 1153
JC_sb_ML = 475 + 470 + 495 
JC_sb_ML # 1440

JC_sb_b1/30 # 46.63%
JC_sb_b2/30 # 47.2%
JC_sb_b3/30 # 47.17%
JC_sb_b10/30 # 47.17%
JC_sb_u/30 # 46.3%
JC_sb_cu/30 # 47.23%
JC_sb_cud/30 # 47.57% 
JC_sb_NJ/30 # 38.43%
JC_sb_ML/30 # 48%
  
# Binomal test JC

prop.test(c(JC_sb_b1,JC_sb_b2),c(3000,3000)) # no diff p-value = 0.6789
prop.test(c(JC_sb_b1,JC_sb_b3),c(3000,3000)) # no diff p-value = 0.698
prop.test(c(JC_sb_b1,JC_sb_b10),c(3000,3000)) # no diff p-value = 0.698
prop.test(c(JC_sb_b1,JC_sb_u),c(3000,3000)) # no diff p-value = 0.8158
prop.test(c(JC_sb_b1,JC_sb_cu),c(3000,3000)) # no diff p-value = 0.6601
prop.test(c(JC_sb_b1,JC_sb_cud),c(3000,3000)) # no diff p-value = 0.485
prop.test(c(JC_sb_b1,JC_sb_NJ),c(3000,3000)) # sign diff p-value = 1.578e-10 
prop.test(c(JC_sb_b1,JC_sb_ML),c(3000,3000)) # no diff  p-value = 0.301

prop.test(c(JC_sb_b2,JC_sb_b3),c(3000,3000)) # no diff p-value = 1
prop.test(c(JC_sb_b2,JC_sb_b10),c(3000,3000)) # no diff p-value = 1
prop.test(c(JC_sb_b2,JC_sb_u),c(3000,3000)) # no diff p-value = 0.5011
prop.test(c(JC_sb_b2,JC_sb_cu),c(3000,3000)) # no diff  p-value = 1
prop.test(c(JC_sb_b2,JC_sb_cud),c(3000,3000)) # no diff p-value = 0.796
prop.test(c(JC_sb_b2,JC_sb_NJ),c(3000,3000)) # sign diff p-value = 8.159e-12
prop.test(c(JC_sb_b2,JC_sb_ML),c(3000,3000)) # no diff p-value = 0.5521

prop.test(c(JC_sb_b3,JC_sb_b10),c(3000,3000)) # no diff p-value = 1
prop.test(c(JC_sb_b3,JC_sb_u),c(3000,3000)) # no diff p-value = 0.5177
prop.test(c(JC_sb_b3,JC_sb_cu),c(3000,3000)) # no diff p-value = 0.9794
prop.test(c(JC_sb_b3,JC_sb_cud),c(3000,3000)) # no diff p-value = 0.7761
prop.test(c(JC_sb_b3,JC_sb_NJ),c(3000,3000)) # sign diff p-value = 9.762e-12
prop.test(c(JC_sb_b3,JC_sb_ML),c(3000,3000)) # no diff p-value = 0.535

prop.test(c(JC_sb_b10,JC_sb_u),c(3000,3000)) # no diff p-value = 0.5177
prop.test(c(JC_sb_b10,JC_sb_cu),c(3000,3000)) # no diff p-value = 0.9794
prop.test(c(JC_sb_b10,JC_sb_cud),c(3000,3000)) # no diff p-value = 0.7761
prop.test(c(JC_sb_b10,JC_sb_NJ),c(3000,3000)) # sign diff p-value = 9.762e-12
prop.test(c(JC_sb_b10,JC_sb_ML),c(3000,3000)) # no diff p-value = 0.535

prop.test(c(JC_sb_u,JC_sb_cu),c(3000,3000)) # no diff p-value = 0.4848
prop.test(c(JC_sb_u,JC_sb_cud),c(3000,3000)) # no diff p-value = 0.3385
prop.test(c(JC_sb_u,JC_sb_NJ),c(3000,3000)) # sign diff p-value = 8.271e-10
prop.test(c(JC_sb_u,JC_sb_ML),c(3000,3000)) # no diff p-value = 0.196

prop.test(c(JC_sb_cu,JC_sb_cud),c(3000,3000)) # no diff p-value = 0.816
prop.test(c(JC_sb_cu,JC_sb_NJ),c(3000,3000)) # sign diff p-value = 6.815e-12
prop.test(c(JC_sb_cu,JC_sb_ML),c(3000,3000)) # no diff p-value = 0.5696

prop.test(c(JC_sb_cud,JC_sb_NJ),c(3000,3000)) # sign diff p-value = 1.088e-12
prop.test(c(JC_sb_cud,JC_sb_ML),c(3000,3000)) # no diff p-value = 0.7565

prop.test(c(JC_sb_NJ,JC_sb_ML),c(3000,3000)) # sign diff p-value = 9.097e-14

# Number of succ. HKY 
HKY_sb_b1 = 461 + 472 + 509
HKY_sb_b1 # 1442
HKY_sb_b2 = 482 + 469 + 477
HKY_sb_b2 # 1428
HKY_sb_b3 = 444 + 471 + 505
HKY_sb_b3 # 1420
HKY_sb_b10 = 419 + 484 + 531
HKY_sb_b10 # 1434
HKY_sb_u = 443 + 471 + 471
HKY_sb_u # 1385
HKY_sb_cu = 479 + 439 + 521
HKY_sb_cu # 1439
HKY_sb_cud = 478 + 453 + 485
HKY_sb_cud # 1416 
HKY_sb_ML = 465 + 473 + 500
HKY_sb_ML # 1438 
HKY_sb_NJ = 461 + 434 + 403 
HKY_sb_NJ # 1298

HKY_sb_b1/30 # 48.07%
HKY_sb_b2/30 # 47.6% 
HKY_sb_b3/30 # 47.33%
HKY_sb_b10/30 # 47.8%
HKY_sb_u/30 # 46.17%
HKY_sb_cu/30 # 47.97%
HKY_sb_cud/30 # 47.2%
HKY_sb_NJ/30 # 43.27%
HKY_sb_ML/30 # 47.93%

# Binomal test HKY 

prop.test(c(HKY_sb_b1,HKY_sb_b2),c(3000,3000)) # no diff, p-value = 0.7369
prop.test(c(HKY_sb_b1,HKY_sb_b3),c(3000,3000)) # no diff, p-value = 0.5873
prop.test(c(HKY_sb_b1,HKY_sb_b10),c(3000,3000)) # no diff, p-value = 0.8565
prop.test(c(HKY_sb_b1,HKY_sb_u),c(3000,3000)) # no diff, p-value = 0.1475
prop.test(c(HKY_sb_b1,HKY_sb_cu),c(3000,3000)) # no diff, p-value = 0.9588
prop.test(c(HKY_sb_b1,HKY_sb_cud),c(3000,3000)) # no diff, p-value = 0.5181
prop.test(c(HKY_sb_b1,HKY_ML),c(3000,3000)) # no diff, p-value = 0.9382
prop.test(c(HKY_sb_b1,HKY_NJ),c(3000,3000)) # sign diff p-value = 0.0002104

prop.test(c(HKY_sb_b2,HKY_sb_b3),c(3000,3000)) # no diff,  p-value = 0.8564
prop.test(c(HKY_sb_b2,HKY_sb_b10),c(3000,3000)) # no diff, p-value = 0.8972
prop.test(c(HKY_sb_b2,HKY_sb_u),c(3000,3000)) # no diff,  p-value = 0.2772
prop.test(c(HKY_sb_b2,HKY_sb_cu),c(3000,3000)) # no diff, p-value = 0.7961
prop.test(c(HKY_sb_b2,HKY_sb_cud),c(3000,3000)) # no diff,  p-value = 0.7761
prop.test(c(HKY_sb_b2,HKY_ML),c(3000,3000)) # no diff, p-value = 0.8161
prop.test(c(HKY_sb_b2,HKY_NJ),c(3000,3000)) # sign diff p-value = 0.0008236
          
prop.test(c(HKY_sb_b3,HKY_sb_b10),c(3000,3000)) # no diff, p-value = 0.7368
prop.test(c(HKY_sb_b3,HKY_sb_u),c(3000,3000)) # no diff, p-value = 0.379
prop.test(c(HKY_sb_b3,HKY_sb_cu),c(3000,3000)) # no diff, p-value = 0.6417
prop.test(c(HKY_sb_b3,HKY_sb_cud),c(3000,3000)) # no diff, p-value = 0.9382
prop.test(c(HKY_sb_b3,HKY_ML),c(3000,3000)) # no diff, p-value = 0.6603
prop.test(c(HKY_sb_b3,HKY_NJ),c(3000,3000)) # sign diff p-value = 0.0017

prop.test(c(HKY_sb_b10,HKY_sb_u),c(3000,3000)) # no diff, p-value = 0.2144
prop.test(c(HKY_sb_b10,HKY_sb_cu),c(3000,3000)) # no diff, p-value = 0.9177
prop.test(c(HKY_sb_b10,HKY_sb_cud),c(3000,3000)) # no diff, p-value = 0.6603
prop.test(c(HKY_sb_b10,HKY_ML),c(3000,3000)) # no diff,  p-value = 0.9382
prop.test(c(HKY_sb_b10,HKY_NJ),c(3000,3000)) # p-value = 0.0004658

prop.test(c(HKY_sb_u,HKY_sb_cu),c(3000,3000)) # no diff, p-value = 0.1704
prop.test(c(HKY_sb_u,HKY_sb_cud),c(3000,3000)) # no diff, p-value = 0.4376
prop.test(c(HKY_sb_u,HKY_ML),c(3000,3000)) # no diff, p-value = 0.1786
prop.test(c(HKY_sb_u,HKY_NJ),c(3000,3000)) # sign diff p-value = 0.02555

prop.test(c(HKY_sb_cu,HKY_sb_cud),c(3000,3000)) # no diff, p-value = 0.5696
prop.test(c(HKY_sb_cu,HKY_ML),c(3000,3000)) # no diff,  p-value = 1
prop.test(c(HKY_sb_cu,HKY_NJ),c(3000,3000)) # sign diff p-value = 0.0002848

prop.test(c(HKY_sb_cud,HKY_ML),c(3000,3000)) # no diff,  p-value = 0.5872
prop.test(c(HKY_sb_cud,HKY_NJ),c(3000,3000)) # sign diff p-value = 0.002407

prop.test(c(HKY_ML,HKY_NJ),c(3000,3000)) # sign diff p-value = 0.0003146

# Num of succ. F81 
F81_sb_b1 = 480 + 486 + 468 
F81_sb_b1 # 1434 
F81_sb_b2 = 463 + 449 + 502 
F81_sb_b2 # 1414
F81_sb_b3 = 455 + 497 + 465
F81_sb_b3 # 1417 
F81_sb_b10 = 428 + 480 + 539 
F81_sb_b10 # 1447 
F81_sb_u = 470 + 469 + 443 
F81_sb_u # 1382 
F81_sb_cu = 481 + 486 + 450 
F81_sb_cu # 1417 
F81_sb_cud = 494 + 465 + 459 
F81_sb_cud # 1418 
F81_sb_ML = 488 + 470 + 475 
F81_sb_ML # 1433 
F81_sb_NJ = 559 + 315 + 303
F81_sb_NJ #  1177

F81_sb_b1/30 # 47.8%
F81_sb_b2/30 # 47.13%
F81_sb_b3/30 # 47.23%
F81_sb_b10/30 # 48.23%
F81_sb_u/30 # 46.07%
F81_sb_cu/30 # 47.23%
F81_sb_cud/30 # 47.27%
F81_sb_NJ/30 # 39.23%
F81_sb_ML/30 # 47.77%

# Binomial test F81 

prop.test(c(F81_sb_b1,F81_sb_b2),c(3000,3000)) # no diff, p-value = 0.6233
prop.test(c(F81_sb_b1,F81_sb_b3),c(3000,3000)) # no diff, p-value = 0.6791
prop.test(c(F81_sb_b1,F81_sb_b10),c(3000,3000)) # no diff, p-value = 0.7565
prop.test(c(F81_sb_b1,F81_sb_u),c(3000,3000)) # no diff, p-value = 0.1871
prop.test(c(F81_sb_b1,F81_sb_cu),c(3000,3000)) # no diff, p-value = 0.6791
prop.test(c(F81_sb_b1,F81_sb_cud),c(3000,3000)) # no diff, p-value = 0.6982
prop.test(c(F81_sb_b1,F81_ML_sb),c(3000,3000)) # no diff, p-value = 1
prop.test(c(F81_sb_b1,F81_NJ_sb),c(3000,3000)) # sign diff p-value = 2.626e-11

prop.test(c(F81_sb_b2,F81_sb_b3),c(3000,3000)) # no diff,  p-value = 0.9588
prop.test(c(F81_sb_b2,F81_sb_b10),c(3000,3000)) # no diff, p-value = 0.4082
prop.test(c(F81_sb_b2,F81_sb_u),c(3000,3000)) # no diff,  p-value = 0.4224
prop.test(c(F81_sb_b2,F81_sb_cu),c(3000,3000)) # no diff, p-value = 0.9588
prop.test(c(F81_sb_b2,F81_sb_cud),c(3000,3000)) # no diff, p-value = 0.9382
prop.test(c(F81_sb_b2,F81_ML_sb),c(3000,3000)) # no diff, p-value = 0.6417
prop.test(c(F81_sb_b2,F81_NJ_sb),c(3000,3000)) # sign diff p-value = 7.703e-10

prop.test(c(F81_sb_b3,F81_sb_b10),c(3000,3000)) # no diff, p-value = 0.4535
prop.test(c(F81_sb_b3,F81_sb_u),c(3000,3000)) # no diff, p-value = 0.3789
prop.test(c(F81_sb_b3,F81_sb_cu),c(3000,3000)) # no diff, p-value = 1
prop.test(c(F81_sb_b3,F81_sb_cud),c(3000,3000)) # no diff, p-value = 1
prop.test(c(F81_sb_b3,F81_ML_sb),c(3000,3000)) # no diff, p-value = 0.6982
prop.test(c(F81_sb_b3,F81_NJ_sb),c(3000,3000)) # sign diff p-value = 4.717e-10

prop.test(c(F81_sb_b10,F81_sb_u),c(3000,3000)) # no diff, pretty close,  p-value = 0.09789
prop.test(c(F81_sb_b10,F81_sb_cu),c(3000,3000)) # no diff, p-value = 0.4535
prop.test(c(F81_sb_b10,F81_sb_cud),c(3000,3000)) # no diff, p-value = 0.4693
prop.test(c(F81_sb_b10,F81_ML_sb),c(3000,3000)) # no diff,   p-value = 0.7369
prop.test(c(F81_sb_b10,F81_NJ_sb),c(3000,3000)) # sign diff p-value = 2.546e-12

prop.test(c(F81_sb_u,F81_sb_cu),c(3000,3000)) # no diff, p-value = 0.3789
prop.test(c(F81_sb_u,F81_sb_cud),c(3000,3000)) # no diff, p-value = 0.3651
prop.test(c(F81_sb_u,F81_ML_sb),c(3000,3000)) # no diff, p-value = 0.1959
prop.test(c(F81_sb_u,F81_NJ_sb),c(3000,3000)) # sign diff p-value = 1.009e-07

prop.test(c(F81_sb_cu,F81_sb_cud),c(3000,3000)) # no diff, p-value = 1 
prop.test(c(F81_sb_cu,F81_ML_sb),c(3000,3000)) # no diff, p-value = 0.6982
prop.test(c(F81_sb_cu,F81_NJ_sb),c(3000,3000)) # sign diff p-value = 4.717e-10

prop.test(c(F81_sb_cud,F81_ML_sb),c(3000,3000)) # no diff, p-value = 0.7174
prop.test(c(F81_sb_cud,F81_NJ_sb),c(3000,3000)) # sign diff p-value = 4e-10

prop.test(c(F81_ML_sb,F81_NJ_sb),c(3000,3000)) # sign diff p-value = 3.128e-11

# Number of succ. F84 
F84_sb_b1 = 443 + 464 + 514 
F84_sb_b1  # 1421
F84_sb_b2 = 460 + 459 + 493 
F84_sb_b2  # 1412
F84_sb_b3 = 469 + 470 + 499
F84_sb_b3 # 1438  
F84_sb_b10 = 379 + 544 + 499 
F84_sb_b10  # 1422
F84_sb_u = 430 + 478 + 458 
F84_sb_u # 1366
F84_sb_cu = 468 + 477 + 480 
F84_sb_cu  # 1425
F84_sb_cud = 434 + 514 + 467 
F84_sb_cud  # 1415 
F84_ML_sb = 477 + 487 + 486 
F84_ML_sb  # 1450 
F84_NJ_sb = 484 + 412 + 403 
F84_NJ_sb #  1299

F84_sb_b1/30 # 47.37%
F84_sb_b2/30 # 47.07%
F84_sb_b3/30 # 47.93%
F84_sb_b10/30 # 47.4%
F84_sb_u/30 # 45.53%
F84_sb_cu/30 # 47.5%
F84_sb_cud/30 # 47.17%
F84_NJ_sb/30 # 43.3%
F84_ML_sb/30 # 48.33%

# Binomial test F84 

prop.test(c(F84_sb_b1,F84_sb_b2),c(3000,3000)) # no diff, p-value = 0.8361
prop.test(c(F84_sb_b1,F84_sb_b3),c(3000,3000)) # no diff, p-value = 0.6792
prop.test(c(F84_sb_b1,F84_sb_b10),c(3000,3000)) # no diff, p-value = 1
prop.test(c(F84_sb_b1,F84_sb_u),c(3000,3000)) # no diff, p-value = 0.1622
prop.test(c(F84_sb_b1,F84_sb_cu),c(3000,3000)) # no diff, p-value = 0.9382
prop.test(c(F84_sb_b1,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.8971
prop.test(c(F84_sb_b1,F84_ML_sb),c(3000,3000)) # no diff, p-value = 0.4693
prop.test(c(F84_sb_b1,F84_NJ_sb),c(3000,3000)) # sign diff p-value = 0.001702

prop.test(c(F84_sb_b2,F84_sb_b3),c(3000,3000)) # no diff,  p-value = 0.5181
prop.test(c(F84_sb_b2,F84_sb_b10),c(3000,3000)) # no diff, p-value = 0.816
prop.test(c(F84_sb_b2,F84_sb_u),c(3000,3000)) # no diff,  p-value = 0.244
prop.test(c(F84_sb_b2,F84_sb_cu),c(3000,3000)) # no diff, p-value = 0.7563
prop.test(c(F84_sb_b2,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.9587
prop.test(c(F84_sb_b2,F84_ML_sb),c(3000,3000)) # no diff, p-value = 0.3389
prop.test(c(F84_sb_b2,F84_NJ_sb),c(3000,3000)) #sign diff p-value = 0.003669

prop.test(c(F84_sb_b3,F84_sb_b10),c(3000,3000)) # no diff, p-value = 0.6982
prop.test(c(F84_sb_b3,F84_sb_u),c(3000,3000)) # no diff, p-value = 0.06619 (very close)
prop.test(c(F84_sb_b3,F84_sb_cu),c(3000,3000)) # no diff, p-value = 0.7564
prop.test(c(F84_sb_b3,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.5695
prop.test(c(F84_sb_b3,F84_ML_sb),c(3000,3000)) # no diff, p-value = 0.7762
prop.test(c(F84_sb_b3,F84_NJ_sb),c(3000,3000)) # sign diff p-value = 0.0003477

prop.test(c(F84_sb_b10,F84_sb_u),c(3000,3000)) # no diff, p-value = 0.1545
prop.test(c(F84_sb_b10,F84_sb_cu),c(3000,3000)) # no diff, p-value = 0.9588
prop.test(c(F84_sb_b10,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.8767
prop.test(c(F84_sb_b10,F84_ML_sb),c(3000,3000)) # no diff,   p-value = 0.4853
prop.test(c(F84_sb_b10,F84_NJ_sb),c(3000,3000)) # sign diff p-value = 0.001558

prop.test(c(F84_sb_u,F84_sb_cu),c(3000,3000)) # no diff, p-value = 0.1333
prop.test(c(F84_sb_u,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.214
prop.test(c(F84_sb_u,F84_ML_sb),c(3000,3000)) # sign diff, p-value = 0.03179
prop.test(c(F84_sb_u,F84_NJ_sb),c(3000,3000)) # no diff p-value = 0.08637

prop.test(c(F84_sb_cu,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.816
prop.test(c(F84_sb_cu,F84_ML_sb),c(3000,3000)) # no diff, p-value = 0.5351
prop.test(c(F84_sb_cu,F84_NJ_sb),c(3000,3000)) # sign diff p-value = 0.00119

prop.test(c(F84_sb_cud,F84_ML_sb),c(3000,3000)) # no diff, p-value = 0.3795
prop.test(c(F84_sb_cud,F84_NJ_sb),c(3000,3000)) # sign diff p-value = 0.002856

prop.test(c(F84_ML_sb,F84_NJ_sb),c(3000,3000)) # sign diff p-value = 0.0001017

# Number of succ. K2P 
K2P_sb_b1 = 477 + 456 + 488 
K2P_sb_b1 # 1421 
K2P_sb_b2 = 472 + 492 + 484 
K2P_sb_b2 #   1448
K2P_sb_b3 = 462 + 487 + 500 
K2P_sb_b3 # 1449 
K2P_sb_b10 = 515 + 470 + 426 
K2P_sb_b10 # 1411
K2P_sb_u = 429 + 494 + 492 
K2P_sb_u # 1415 
K2P_sb_cu = 451 + 502 + 481 
K2P_sb_cu # 1434
K2P_sb_cud = 428 + 509 + 496 
K2P_sb_cud # 1433
K2P_NJ_sb = 466 + 426 + 389
K2P_NJ_sb #  1281 
K2P_ML_sb = 475 + 470 + 495 
K2P_ML_sb # 1440

K2P_sb_b1/30 # 47.37%
K2P_sb_b2/30 # 48.27%
K2P_sb_b3/30 # 48.3%
K2P_sb_b10/30 # 47.03%
K2P_sb_u/30 # 47.17%
K2P_sb_cu/30 # 47.8%
K2P_sb_cud/30 # 47.77%
K2P_NJ_sb/30 # 42.7%
K2P_ML_sb/30 # 48%

# Binomial test K2P 

prop.test(c(K2P_sb_b1,K2P_sb_b2),c(3000,3000)) # no diff p-value = 0.5016
prop.test(c(K2P_sb_b1,K2P_sb_b3),c(3000,3000)) # no diff,  p-value = 0.4853
prop.test(c(K2P_sb_b1,K2P_sb_b10),c(3000,3000)) # no diff, p-value = 0.816
prop.test(c(K2P_sb_b1,K2P_sb_u),c(3000,3000)) # no diff,   p-value = 0.8971
prop.test(c(K2P_sb_b1,K2P_sb_cu),c(3000,3000)) # no diff, p-value = 0.7564
prop.test(c(K2P_sb_b1,K2P_sb_cud),c(3000,3000)) # no diff, p-value = 0.7761
prop.test(c(K2P_sb_b1,K2P_NJ_sb),c(3000,3000)) # sign diff p-value = 0.00031
prop.test(c(K2P_sb_b1,K2P_ML_sb),c(3000,3000)) # no diff, p-value = 0.6417

prop.test(c(K2P_sb_b2,K2P_sb_b3),c(3000,3000)) # no diff p-value = 1
prop.test(c(K2P_sb_b2,K2P_sb_b10),c(3000,3000)) # no diff p-value = 0.3521
prop.test(c(K2P_sb_b2,K2P_sb_u),c(3000,3000)) # no diff p-value = 0.4082
prop.test(c(K2P_sb_b2,K2P_sb_cu),c(3000,3000)) # no diff p-value = 0.7369
prop.test(c(K2P_sb_b2,K2P_sb_cud),c(3000,3000)) # no diff p-value = 0.7175
prop.test(c(K2P_sb_b2,K2P_NJ_sb),c(3000,3000)) # sign diff p-value = 1.68e-05
prop.test(c(K2P_sb_b2,K2P_ML_sb),c(3000,3000)) # no diff p-value = 0.8565

prop.test(c(K2P_sb_b3,K2P_sb_b10),c(3000,3000)) # no diff, p-value = 0.3389
prop.test(c(K2P_sb_b3,K2P_sb_u),c(3000,3000)) # no diff, p-value = 0.3937
prop.test(c(K2P_sb_b3,K2P_sb_cu),c(3000,3000)) # no diff, p-value = 0.7175
prop.test(c(K2P_sb_b3,K2P_sb_cud),c(3000,3000)) # no diff, p-value = 0.6983
prop.test(c(K2P_sb_b3,K2P_NJ_sb),c(3000,3000)) # sign diff p-value = 1.495e-05
prop.test(c(K2P_sb_b3,K2P_ML_sb),c(3000,3000)) # no diff, p-value = 0.8362

prop.test(c(K2P_sb_b10,K2P_sb_u),c(3000,3000)) # no diff, p-value = 0.9382
prop.test(c(K2P_sb_b10,K2P_sb_cu),c(3000,3000)) # no diff, p-value = 0.5695
prop.test(c(K2P_sb_b10,K2P_sb_cud),c(3000,3000)) # no diff, p-value = 0.5872
prop.test(c(K2P_sb_b10,K2P_NJ_sb),c(3000,3000)) # sign diff p-value = 0.0008126
prop.test(c(K2P_sb_b10,K2P_ML_sb),c(3000,3000)) # no diff,  p-value = 0.4692

prop.test(c(K2P_sb_u,K2P_sb_cu),c(3000,3000)) # no diff, p-value = 0.6417
prop.test(c(K2P_sb_u,K2P_sb_cud),c(3000,3000)) # no diff, p-value = 0.6603
prop.test(c(K2P_sb_u,K2P_NJ_sb),c(3000,3000)) # sign diff p-value = 0.0005568
prop.test(c(K2P_sb_u,K2P_ML_sb),c(3000,3000)) # no diff, p-value = 0.535

prop.test(c(K2P_sb_cu,K2P_sb_cud),c(3000,3000)) # no diff, p-value = 1
prop.test(c(K2P_sb_cu,K2P_NJ_sb),c(3000,3000)) # sign diff p-value = 8.065e-05
prop.test(c(K2P_sb_cu,K2P_ML_sb),c(3000,3000)) # no diff, p-value = 0.8972

prop.test(c(K2P_sb_cud,K2P_NJ_sb),c(3000,3000)) # sign diff p-value = 8.979e-05
prop.test(c(K2P_sb_cud,K2P_ML_sb),c(3000,3000)) # no diff, p-value = 0.8768

prop.test(c(K2P_NJ_sb,K2P_ML_sb),c(3000,3000)) # sign diff p-value = 4.18e-05

# Numb of suc GTR 
GTR_sb_b1 = 458 + 494 + 447 
GTR_sb_b1 # 1399 
GTR_sb_b2 = 460 + 497 + 452 
GTR_sb_b2 # 1409
GTR_sb_b3 = 467 + 480 + 491 
GTR_sb_b3 # 1438
GTR_sb_b10 = 501 + 467 + 456 
GTR_sb_b10 # 1424
GTR_sb_u = 432 + 494 + 499 
GTR_sb_u # 1425
GTR_sb_cu = 463 + 499 + 455 
GTR_sb_cu # 1417
GTR_sb_cud = 409 + 507 + 491 
GTR_sb_cud # 1407
GTR_sb_ML = 466 + 490 + 495 
GTR_sb_ML # 1451
GTR_sb_NJ = 483 + 412 + 403
GTR_sb_NJ # 1298

GTR_sb_b1/30 # 46.63%
GTR_sb_b2/30 # 46.97%
GTR_sb_b3/30 # 47.93%
GTR_sb_b10/30 # 47.47%
GTR_sb_u/30 # 47.5%
GTR_sb_cu/30 # 47.23%
GTR_sb_cud/30 # 46.9%
GTR_sb_NJ/30 # 43.27%
GTR_sb_ML/30 # 48.37%

# Binomial test GTR_sb

prop.test(c(GTR_sb_b1,GTR_sb_b2),c(3000,3000)) # no diff p-value = 0.8159
prop.test(c(GTR_sb_b1,GTR_sb_b3),c(3000,3000)) # no diff p-value = 0.3258
prop.test(c(GTR_sb_b1,GTR_sb_b10),c(3000,3000)) # no diff p-value = 0.5348
prop.test(c(GTR_sb_b1,GTR_sb_u),c(3000,3000)) # no diff p-value = 0.5179
prop.test(c(GTR_sb_b1,GTR_sb_cu),c(3000,3000)) # no diff p-value = 0.6601
prop.test(c(GTR_sb_b1,GTR_sb_cud),c(3000,3000)) # no diff p-value = 0.8563
prop.test(c(GTR_sb_b1,GTR_sb_ML),c(3000,3000)) # no diff p-value = 0.1873
prop.test(c(GTR_sb_b1,GTR_sb_NJ),c(3000,3000)) # sign diff p-value = 0.009452

prop.test(c(GTR_sb_b2,GTR_sb_b3),c(3000,3000)) # no diff p-value = 0.4691
prop.test(c(GTR_sb_b2,GTR_sb_b10),c(3000,3000)) # no diff p-value = 0.7173
prop.test(c(GTR_sb_b2,GTR_sb_u),c(3000,3000)) # no diff p-value = 0.6981
prop.test(c(GTR_sb_b2,GTR_sb_cu),c(3000,3000)) # no diff p-value = 0.8563
prop.test(c(GTR_sb_b2,GTR_sb_cud),c(3000,3000)) # no diff p-value = 0.9794
prop.test(c(GTR_sb_b2,GTR_sb_ML),c(3000,3000)) # no diff p-value = 0.2892
prop.test(c(GTR_sb_b2,GTR_sb_NJ),c(3000,3000)) # sign diff p-value = 0.00432

prop.test(c(GTR_sb_b3,GTR_sb_b10),c(3000,3000)) # no diff p-value = 0.7369
prop.test(c(GTR_sb_b3,GTR_sb_u),c(3000,3000)) # no diff  p-value = 0.7564
prop.test(c(GTR_sb_b3,GTR_sb_cu),c(3000,3000)) # no diff  p-value = 0.6052
prop.test(c(GTR_sb_b3,GTR_sb_cud),c(3000,3000)) # no diff p-value = 0.438
prop.test(c(GTR_sb_b3,GTR_sb_ML),c(3000,3000)) # no diff p-value = 0.7565
prop.test(c(GTR_sb_b3,GTR_sb_NJ),c(3000,3000)) # sign diff p-value = 0.0003146

prop.test(c(GTR_sb_b10,GTR_sb_u),c(3000,3000)) # no diff p-value = 1
prop.test(c(GTR_sb_b10,GTR_sb_cu),c(3000,3000)) #no diff p-value = 0.8767
prop.test(c(GTR_sb_b10,GTR_sb_cud),c(3000,3000)) # no diff p-value = 0.679
prop.test(c(GTR_sb_b10,GTR_sb_ML),c(3000,3000)) # no diff p-value = 0.5016
prop.test(c(GTR_sb_b10,GTR_sb_NJ),c(3000,3000)) # sign diff p-value = 0.001189

prop.test(c(GTR_sb_u,GTR_sb_cu),c(3000,3000)) # no diff  p-value = 0.8564
prop.test(c(GTR_sb_u,GTR_sb_cud),c(3000,3000)) # no diff p-value = 0.6602
prop.test(c(GTR_sb_u,GTR_sb_ML),c(3000,3000)) # no diff p-value = 0.5182
prop.test(c(GTR_sb_u,GTR_sb_NJ),c(3000,3000)) # sign diff p-value = 0.001086

prop.test(c(GTR_sb_cu,GTR_sb_cud),c(3000,3000)) # no diff p-value = 0.8159
prop.test(c(GTR_sb_cu,GTR_sb_ML),c(3000,3000)) # no diff p-value = 0.3937
prop.test(c(GTR_sb_cu,GTR_sb_NJ),c(3000,3000)) # sign diff p-value = 0.002209

prop.test(c(GTR_sb_cud,GTR_sb_ML),c(3000,3000)) # no diff p-value = 0.2664
prop.test(c(GTR_sb_cud,GTR_sb_NJ),c(3000,3000)) # sign diff p-value = 0.005077

prop.test(c(GTR_sb_ML,GTR_sb_NJ),c(3000,3000)) # sign diff p-value = 8.202e-05

