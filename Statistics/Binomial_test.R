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
HKY_ML = 465 + 473 + 500
HKY_ML # 1438 

# Binomal test HKY 

prop.test(c(HKY_sb_b1,HKY_sb_b2),c(3000,3000)) # no diff, p-value = 0.7369
prop.test(c(HKY_sb_b1,HKY_sb_b3),c(3000,3000)) # no diff, p-value = 0.5873
prop.test(c(HKY_sb_b1,HKY_sb_b10),c(3000,3000)) # no diff, p-value = 0.8565
prop.test(c(HKY_sb_b1,HKY_sb_u),c(3000,3000)) # no diff, p-value = 0.1475
prop.test(c(HKY_sb_b1,HKY_sb_cu),c(3000,3000)) # no diff, p-value = 0.9588
prop.test(c(HKY_sb_b1,HKY_sb_cud),c(3000,3000)) # no diff, p-value = 0.5181
prop.test(c(HKY_sb_b1,HKY_ML),c(3000,3000)) # no diff, p-value = 0.9382

prop.test(c(HKY_sb_b2,HKY_sb_b3),c(3000,3000)) # no diff,  p-value = 0.8564
prop.test(c(HKY_sb_b2,HKY_sb_b10),c(3000,3000)) # no diff, p-value = 0.8972
prop.test(c(HKY_sb_b2,HKY_sb_u),c(3000,3000)) # no diff,  p-value = 0.2772
prop.test(c(HKY_sb_b2,HKY_sb_cu),c(3000,3000)) # no diff, p-value = 0.7961
prop.test(c(HKY_sb_b2,HKY_sb_cud),c(3000,3000)) # no diff,  p-value = 0.7761
prop.test(c(HKY_sb_b2,HKY_ML),c(3000,3000)) # no diff, p-value = 0.8161
          
prop.test(c(HKY_sb_b3,HKY_sb_b10),c(3000,3000)) # no diff, p-value = 0.7368
prop.test(c(HKY_sb_b3,HKY_sb_u),c(3000,3000)) # no diff, p-value = 0.379
prop.test(c(HKY_sb_b3,HKY_sb_cu),c(3000,3000)) # no diff, p-value = 0.6417
prop.test(c(HKY_sb_b3,HKY_sb_cud),c(3000,3000)) # no diff, p-value = 0.9382
prop.test(c(HKY_sb_b3,HKY_ML),c(3000,3000)) # no diff, p-value = 0.6603

prop.test(c(HKY_sb_b10,HKY_sb_u),c(3000,3000)) # no diff, p-value = 0.2144
prop.test(c(HKY_sb_b10,HKY_sb_cu),c(3000,3000)) # no diff, p-value = 0.9177
prop.test(c(HKY_sb_b10,HKY_sb_cud),c(3000,3000)) # no diff, p-value = 0.6603
prop.test(c(HKY_sb_b10,HKY_ML),c(3000,3000)) # no diff,  p-value = 0.9382

prop.test(c(HKY_sb_u,HKY_sb_cu),c(3000,3000)) # no diff, p-value = 0.1704
prop.test(c(HKY_sb_u,HKY_sb_cud),c(3000,3000)) # no diff, p-value = 0.4376
prop.test(c(HKY_sb_u,HKY_ML),c(3000,3000)) # no diff, p-value = 0.1786

prop.test(c(HKY_sb_cu,HKY_sb_cud),c(3000,3000)) # no diff, p-value = 0.5696
prop.test(c(HKY_sb_cu,HKY_ML),c(3000,3000)) # no diff,  p-value = 1

prop.test(c(HKY_sb_cud,HKY_ML),c(3000,3000)) # no diff,  p-value = 0.5872

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
F81_ML_sb = 488 + 470 + 475 
F81_ML_sb # 1433 

# Binomial test F81 

prop.test(c(F81_sb_b1,F81_sb_b2),c(3000,3000)) # no diff, p-value = 0.6233
prop.test(c(F81_sb_b1,F81_sb_b3),c(3000,3000)) # no diff, p-value = 0.6791
prop.test(c(F81_sb_b1,F81_sb_b10),c(3000,3000)) # no diff, p-value = 0.7565
prop.test(c(F81_sb_b1,F81_sb_u),c(3000,3000)) # no diff, p-value = 0.1871
prop.test(c(F81_sb_b1,F81_sb_cu),c(3000,3000)) # no diff, p-value = 0.6791
prop.test(c(F81_sb_b1,F81_sb_cud),c(3000,3000)) # no diff, p-value = 0.6982
prop.test(c(F81_sb_b1,F81_ML_sb),c(3000,3000)) # no diff, p-value = 1


prop.test(c(F81_sb_b2,F81_sb_b3),c(3000,3000)) # no diff,  p-value = 0.9588
prop.test(c(F81_sb_b2,F81_sb_b10),c(3000,3000)) # no diff, p-value = 0.4082
prop.test(c(F81_sb_b2,F81_sb_u),c(3000,3000)) # no diff,  p-value = 0.4224
prop.test(c(F81_sb_b2,F81_sb_cu),c(3000,3000)) # no diff, p-value = 0.9588
prop.test(c(F81_sb_b2,F81_sb_cud),c(3000,3000)) # no diff, p-value = 0.9382
prop.test(c(F81_sb_b2,F81_ML_sb),c(3000,3000)) # no diff, p-value = 0.6417

prop.test(c(F81_sb_b3,F81_sb_b10),c(3000,3000)) # no diff, p-value = 0.4535
prop.test(c(F81_sb_b3,F81_sb_u),c(3000,3000)) # no diff, p-value = 0.3789
prop.test(c(F81_sb_b3,F81_sb_cu),c(3000,3000)) # no diff, p-value = 1
prop.test(c(F81_sb_b3,F81_sb_cud),c(3000,3000)) # no diff, p-value = 1
prop.test(c(F81_sb_b3,F81_ML_sb),c(3000,3000)) # no diff, p-value = 0.6982

prop.test(c(F81_sb_b10,F81_sb_u),c(3000,3000)) # no diff, pretty close,  p-value = 0.09789
prop.test(c(F81_sb_b10,F81_sb_cu),c(3000,3000)) # no diff, p-value = 0.4535
prop.test(c(F81_sb_b10,F81_sb_cud),c(3000,3000)) # no diff, p-value = 0.4693
prop.test(c(F81_sb_b10,F81_ML_sb),c(3000,3000)) # no diff,   p-value = 0.7369

prop.test(c(F81_sb_u,F81_sb_cu),c(3000,3000)) # no diff, p-value = 0.3789
prop.test(c(F81_sb_u,F81_sb_cud),c(3000,3000)) # no diff, p-value = 0.3651
prop.test(c(F81_sb_u,F81_ML_sb),c(3000,3000)) # no diff, p-value = 0.1959

prop.test(c(F81_sb_cu,F81_sb_cud),c(3000,3000)) # no diff, p-value = 1 
prop.test(c(F81_sb_cu,F81_ML_sb),c(3000,3000)) # no diff, p-value = 0.6982

prop.test(c(F81_sb_cud,F81_ML_sb),c(3000,3000)) # no diff, p-value = 0.7174

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

# Binomial test F84 

prop.test(c(F84_sb_b1,F84_sb_b2),c(3000,3000)) # no diff, p-value = 0.8361
prop.test(c(F84_sb_b1,F84_sb_b3),c(3000,3000)) # no diff, p-value = 0.6792
prop.test(c(F84_sb_b1,F84_sb_b10),c(3000,3000)) # no diff, p-value = 1
prop.test(c(F84_sb_b1,F84_sb_u),c(3000,3000)) # no diff, p-value = 0.1622
prop.test(c(F84_sb_b1,F84_sb_cu),c(3000,3000)) # no diff, p-value = 0.9382
prop.test(c(F84_sb_b1,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.8971
prop.test(c(F84_sb_b1,F84_ML_sb),c(3000,3000)) # no diff, p-value = 0.4693

prop.test(c(F84_sb_b2,F84_sb_b3),c(3000,3000)) # no diff,  p-value = 0.5181
prop.test(c(F84_sb_b2,F84_sb_b10),c(3000,3000)) # no diff, p-value = 0.816
prop.test(c(F84_sb_b2,F84_sb_u),c(3000,3000)) # no diff,  p-value = 0.244
prop.test(c(F84_sb_b2,F84_sb_cu),c(3000,3000)) # no diff, p-value = 0.7563
prop.test(c(F84_sb_b2,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.9587
prop.test(c(F84_sb_b2,F84_ML_sb),c(3000,3000)) # no diff, p-value = 0.3389

prop.test(c(F84_sb_b3,F84_sb_b10),c(3000,3000)) # no diff, p-value = 0.6982
prop.test(c(F84_sb_b3,F84_sb_u),c(3000,3000)) # no diff, p-value = 0.06619 (very close)
prop.test(c(F84_sb_b3,F84_sb_cu),c(3000,3000)) # no diff, p-value = 0.7564
prop.test(c(F84_sb_b3,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.5695
prop.test(c(F84_sb_b3,F84_ML_sb),c(3000,3000)) # no diff, p-value = 0.7762

prop.test(c(F84_sb_b10,F84_sb_u),c(3000,3000)) # no diff, p-value = 0.1545
prop.test(c(F84_sb_b10,F84_sb_cu),c(3000,3000)) # no diff, p-value = 0.9588
prop.test(c(F84_sb_b10,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.8767
prop.test(c(F84_sb_b10,F84_ML_sb),c(3000,3000)) # no diff,   p-value = 0.4853

prop.test(c(F84_sb_u,F84_sb_cu),c(3000,3000)) # no diff, p-value = 0.1333
prop.test(c(F84_sb_u,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.214
prop.test(c(F84_sb_u,F84_ML_sb),c(3000,3000)) # sign diff, p-value = 0.03179

prop.test(c(F84_sb_cu,F84_sb_cud),c(3000,3000)) # no diff, p-value = 0.816
prop.test(c(F84_sb_cu,F84_ML_sb),c(3000,3000)) # no diff, p-value = 0.5351

prop.test(c(F84_sb_cud,F84_ML_sb),c(3000,3000)) # no diff, p-value = 0.3795

# Number of succ. K2P 
K2P_sb_b1 = 477 + 456 + 488 
K2P_sb_b1 # 1421 
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
K2P_ML_sb = 475 + 470 + 495 
K2P_ML_sb # 1440

# Binomial test K2P 

prop.test(c(K2P_sb_b1,K2P_sb_b3),c(3000,3000)) # no diff,  p-value = 0.4853
prop.test(c(K2P_sb_b1,K2P_sb_b10),c(3000,3000)) # no diff, p-value = 0.816
prop.test(c(K2P_sb_b1,K2P_sb_u),c(3000,3000)) # no diff,   p-value = 0.8971
prop.test(c(K2P_sb_b1,K2P_sb_cu),c(3000,3000)) # no diff, p-value = 0.7564
prop.test(c(K2P_sb_b1,K2P_sb_cud),c(3000,3000)) # no diff, p-value = 0.7761
prop.test(c(K2P_sb_b1,K2P_ML_sb),c(3000,3000)) # no diff, p-value = 0.6417

prop.test(c(K2P_sb_b3,K2P_sb_b10),c(3000,3000)) # no diff, p-value = 0.3389
prop.test(c(K2P_sb_b3,K2P_sb_u),c(3000,3000)) # no diff, p-value = 0.3937
prop.test(c(K2P_sb_b3,K2P_sb_cu),c(3000,3000)) # no diff, p-value = 0.7175
prop.test(c(K2P_sb_b3,K2P_sb_cud),c(3000,3000)) # no diff, p-value = 0.6983
prop.test(c(K2P_sb_b3,K2P_ML_sb),c(3000,3000)) # no diff, p-value = 0.8362

prop.test(c(K2P_sb_b10,K2P_sb_u),c(3000,3000)) # no diff, p-value = 0.9382
prop.test(c(K2P_sb_b10,K2P_sb_cu),c(3000,3000)) # no diff, p-value = 0.5695
prop.test(c(K2P_sb_b10,K2P_sb_cud),c(3000,3000)) # no diff, p-value = 0.5872
prop.test(c(K2P_sb_b10,K2P_ML_sb),c(3000,3000)) # no diff,  p-value = 0.4692

prop.test(c(K2P_sb_u,K2P_sb_cu),c(3000,3000)) # no diff, p-value = 0.6417
prop.test(c(K2P_sb_u,K2P_sb_cud),c(3000,3000)) # no diff, p-value = 0.6603
prop.test(c(K2P_sb_u,K2P_ML_sb),c(3000,3000)) # no diff, p-value = 0.535

prop.test(c(K2P_sb_cu,K2P_sb_cud),c(3000,3000)) # no diff, p-value = 1
prop.test(c(K2P_sb_cu,K2P_ML_sb),c(3000,3000)) # no diff, p-value = 0.8972

prop.test(c(K2P_sb_cud,K2P_ML_sb),c(3000,3000)) # no diff, p-value = 0.8768
