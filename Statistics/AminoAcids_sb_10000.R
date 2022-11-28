# WAG sb 10000 succ. rate
WAG_sb_10000_b1 = 429 + 529 + 449 
WAG_sb_10000_b1 # 1407
WAG_sb_10000_b2 = 471 + 546 + 470 
WAG_sb_10000_b2 # 1487
WAG_sb_10000_b3 = 500 + 466 + 534
WAG_sb_10000_b3 # 1500
WAG_sb_10000_b10 = 421 + 502 + 411
WAG_sb_10000_b10 # 1334
WAG_sb_10000_u = 472 + 504 + 463 
WAG_sb_10000_u # 1439
WAG_sb_10000_cu = 512 + 450 + 422 
WAG_sb_10000_cu # 1384
WAG_sb_10000_cud = 485 + 448 + 454 
WAG_sb_10000_cud # 1387
WAG_sb_10000_NJ = 619 + 603 + 577 
WAG_sb_10000_NJ # 1799
WAG_sb_10000_ML = 865 + 875 + 870 
WAG_sb_10000_ML #  2610

WAG_sb_10000_b1/30 # 46.9%
WAG_sb_10000_b2/30 # 49.57%
WAG_sb_10000_b3/30 # 50.0%
WAG_sb_10000_b10/30 # 44.47%
WAG_sb_10000_u/30 # 47.97%
WAG_sb_10000_cu/30 # 46.13%
WAG_sb_10000_cud/30 # 46.23%
WAG_sb_10000_NJ/30 # 59.97%
WAG_sb_10000_ML/30 # 87.0%

# WAG_sb 10000 binomial test 

prop.test(c(WAG_sb_10000_b1 ,WAG_sb_10000_b2),c(3000,3000)) # sign diff p-value = 0.04125
prop.test(c(WAG_sb_10000_b1 ,WAG_sb_10000_b3),c(3000,3000)) # sign diff p-value = 0.01747
prop.test(c(WAG_sb_10000_b1 ,WAG_sb_10000_b10),c(3000,3000)) # almost diff p-value = 0.06204
prop.test(c(WAG_sb_10000_b1 ,WAG_sb_10000_u),c(3000,3000)) # no diff  p-value = 0.4229
prop.test(c(WAG_sb_10000_b1 ,WAG_sb_10000_cu),c(3000,3000)) # no diff p-value = 0.5691
prop.test(c(WAG_sb_10000_b1 ,WAG_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.6229
prop.test(c(WAG_sb_10000_b1 ,WAG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_sb_10000_b1 ,WAG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_sb_10000_b2 ,WAG_sb_10000_b3),c(3000,3000)) # no diff p-value = 0.7567
prop.test(c(WAG_sb_10000_b2 ,WAG_sb_10000_b10),c(3000,3000)) # sign diff p-value = 8.437e-05
prop.test(c(WAG_sb_10000_b2 ,WAG_sb_10000_u),c(3000,3000)) # no diff  p-value = 0.2248
prop.test(c(WAG_sb_10000_b2 ,WAG_sb_10000_cu),c(3000,3000)) # sign diff p-value = 0.008387
prop.test(c(WAG_sb_10000_b2 ,WAG_sb_10000_cud),c(3000,3000)) # sign diff p-value = 0.01051
prop.test(c(WAG_sb_10000_b2 ,WAG_sb_10000_NJ),c(3000,3000)) # sign diff p-value = 7.221e-16
prop.test(c(WAG_sb_10000_b2 ,WAG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_sb_10000_b3 ,WAG_sb_10000_b10),c(3000,3000)) # sign diff p-value = 1.983e-05
prop.test(c(WAG_sb_10000_b3 ,WAG_sb_10000_u),c(3000,3000)) # no diff p-value = 0.1213
prop.test(c(WAG_sb_10000_b3 ,WAG_sb_10000_cu),c(3000,3000)) # sign diff p-value = 0.002963
prop.test(c(WAG_sb_10000_b3 ,WAG_sb_10000_cud),c(3000,3000)) # sign diff p-value = 0.003805
prop.test(c(WAG_sb_10000_b3 ,WAG_sb_10000_NJ),c(3000,3000)) # sign diff p-value = 1.052e-14
prop.test(c(WAG_sb_10000_b3 ,WAG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16


prop.test(c(WAG_sb_10000_b10 ,WAG_sb_10000_u),c(3000,3000)) # sign diff p-value = 0.007081
prop.test(c(WAG_sb_10000_b10 ,WAG_sb_10000_cu),c(3000,3000)) # no diff p-value = 0.2038
prop.test(c(WAG_sb_10000_b10 ,WAG_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.1775
prop.test(c(WAG_sb_10000_b10 ,WAG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_sb_10000_b10 ,WAG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_sb_10000_u ,WAG_sb_10000_cu),c(3000,3000)) # no diff p-value = 0.1625
prop.test(c(WAG_sb_10000_u ,WAG_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.1872
prop.test(c(WAG_sb_10000_u ,WAG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_sb_10000_u ,WAG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_sb_10000_cu ,WAG_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.9587
prop.test(c(WAG_sb_10000_cu ,WAG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_sb_10000_cu ,WAG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(WAG_sb_10000_cud ,WAG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(WAG_sb_10000_cud ,WAG_sb_10000_ML),c(3000,3000)) #sign diff p-value < 2.2e-16

prop.test(c(WAG_sb_10000_NJ ,WAG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# DAY sb 10000 succ. rate

DAY_sb_10000_b1 = 488 + 436 + 512 
DAY_sb_10000_b1 # 1436
DAY_sb_10000_b2 = 536 + 527 + 419 
DAY_sb_10000_b2 # 1482
DAY_sb_10000_b3 = 467 + 580 + 410 
DAY_sb_10000_b3 # 1457
DAY_sb_10000_b10 = 458 + 470 + 451 
DAY_sb_10000_b10 # 1379
DAY_sb_10000_u = 546 + 436 + 478
DAY_sb_10000_u # 1460
DAY_sb_10000_cu = 449 + 456 + 439 
DAY_sb_10000_cu # 1344
DAY_sb_10000_cud = 519 + 392 + 477
DAY_sb_10000_cud #  1388
DAY_sb_10000_NJ = 608 + 623 + 588
DAY_sb_10000_NJ # 1819
DAY_sb_10000_ML = 843 + 884 + 871 
DAY_sb_10000_ML #  2598

DAY_sb_10000_b1/30 # 47.87%
DAY_sb_10000_b2/30 # 49.4%
DAY_sb_10000_b3/30 # 48.57%
DAY_sb_10000_b10/30 # 45.97%
DAY_sb_10000_u/30 # 48.67%
DAY_sb_10000_cu/30 # 44.8%
DAY_sb_10000_cud/30 # 46.27%
DAY_sb_10000_NJ/30 # 60.63%
DAY_sb_10000_ML/30 # 86.6%

# DAY_sb 10000 binomial test 
  
prop.test(c(DAY_sb_10000_b1 ,DAY_sb_10000_b2),c(3000,3000)) # no diff p-value = 0.2451
prop.test(c(DAY_sb_10000_b1 ,DAY_sb_10000_b3),c(3000,3000)) # no diff p-value = 0.6053
prop.test(c(DAY_sb_10000_b1 ,DAY_sb_10000_b10),c(3000,3000)) # no diff p-value = 0.1474
prop.test(c(DAY_sb_10000_b1 ,DAY_sb_10000_u),c(3000,3000)) # no diff p-value = 0.5524
prop.test(c(DAY_sb_10000_b1 ,DAY_sb_10000_cu),c(3000,3000)) # sign diff p-value = 0.01848
prop.test(c(DAY_sb_10000_b1 ,DAY_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.2241
prop.test(c(DAY_sb_10000_b1 ,DAY_sb_10000_NJ),c(3000,3000)) # no diff p-value < 2.2e-16
prop.test(c(DAY_sb_10000_b1 ,DAY_sb_10000_ML),c(3000,3000)) # no diff p-value < 2.2e-16

prop.test(c(DAY_sb_10000_b2 ,DAY_sb_10000_b3),c(3000,3000)) # no diff p-value = 0.5354
prop.test(c(DAY_sb_10000_b2 ,DAY_sb_10000_b10),c(3000,3000)) # sign diff p-value = 0.008378
prop.test(c(DAY_sb_10000_b2 ,DAY_sb_10000_u),c(3000,3000)) # no diff p-value = 0.5876
prop.test(c(DAY_sb_10000_b2 ,DAY_sb_10000_cu),c(3000,3000)) # sign diff p-value = 0.0003952
prop.test(c(DAY_sb_10000_b2 ,DAY_sb_10000_cud),c(3000,3000)) # sign diff p-value = 0.01624
prop.test(c(DAY_sb_10000_b2 ,DAY_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_sb_10000_b2 ,DAY_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_sb_10000_b3 ,DAY_sb_10000_b10),c(3000,3000)) # sign diff p-value = 0.04647
prop.test(c(DAY_sb_10000_b3 ,DAY_sb_10000_u),c(3000,3000)) # no diff p-value = 0.9588
prop.test(c(DAY_sb_10000_b3 ,DAY_sb_10000_cu),c(3000,3000)) # sign diff p-value = 0.003753
prop.test(c(DAY_sb_10000_b3 ,DAY_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.07873
prop.test(c(DAY_sb_10000_b3 ,DAY_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_sb_10000_b3 ,DAY_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_sb_10000_b10 ,DAY_sb_10000_u),c(3000,3000)) # sign diff p-value = 0.03859
prop.test(c(DAY_sb_10000_b10 ,DAY_sb_10000_cu),c(3000,3000)) # no diff p-value = 0.378
prop.test(c(DAY_sb_10000_b10 ,DAY_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.8359
prop.test(c(DAY_sb_10000_b10 ,DAY_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_sb_10000_b10 ,DAY_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_sb_10000_u ,DAY_sb_10000_cu),c(3000,3000)) # sign diff p-value = 0.002924
prop.test(c(DAY_sb_10000_u ,DAY_sb_10000_cud),c(3000,3000)) # almost diff p-value = 0.06642
prop.test(c(DAY_sb_10000_u ,DAY_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_sb_10000_u ,DAY_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_sb_10000_cu ,DAY_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.265
prop.test(c(DAY_sb_10000_cu ,DAY_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_sb_10000_cu ,DAY_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_sb_10000_cud ,DAY_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(DAY_sb_10000_cud ,DAY_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(DAY_sb_10000_NJ ,DAY_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16


# JTT sb 10000 succ. rate

JTT_sb_10000_b1 = 541 + 456 + 473 
JTT_sb_10000_b1 # 1470
JTT_sb_10000_b2 = 491 + 482 + 539
JTT_sb_10000_b2 # 1512
JTT_sb_10000_b3 = 562 + 489 + 442 
JTT_sb_10000_b3 # 1493
JTT_sb_10000_b10 = 473 + 494 + 375
JTT_sb_10000_b10 # 1342
JTT_sb_10000_u = 517 + 443 + 492 
JTT_sb_10000_u # 1452
JTT_sb_10000_cu = 429 + 525 + 473 
JTT_sb_10000_cu # 1427
JTT_sb_10000_cud = 465 + 463 + 452 
JTT_sb_10000_cud # 1380
JTT_sb_10000_NJ = 628 + 610 + 597 
JTT_sb_10000_NJ # 1835
JTT_sb_10000_ML = 859 + 891 + 882 
JTT_sb_10000_ML #  2632

JTT_sb_10000_b1/30 # 49.0%
JTT_sb_10000_b2/30 # 50.4%
JTT_sb_10000_b3/30 # 49.77%
JTT_sb_10000_b10/30 # 44.73%
JTT_sb_10000_u/30 # 48.4%
JTT_sb_10000_cu/30 # 47.57%
JTT_sb_10000_cud/30 # 46.0%
JTT_sb_10000_NJ/30 # 61.17%
JTT_sb_10000_ML/30 # 87.73%

# JTT_sb 10000 binomial test 

prop.test(c(JTT_sb_10000_b1 ,JTT_sb_10000_b2),c(3000,3000)) # no diff p-value = 0.2898
prop.test(c(JTT_sb_10000_b1 ,JTT_sb_10000_b3),c(3000,3000)) # no diff p-value = 0.57
prop.test(c(JTT_sb_10000_b1 ,JTT_sb_10000_b10),c(3000,3000)) # sign diff p-value = 0.001018
prop.test(c(JTT_sb_10000_b1 ,JTT_sb_10000_u),c(3000,3000)) # no diff p-value = 0.6606
prop.test(c(JTT_sb_10000_b1 ,JTT_sb_10000_cu),c(3000,3000)) # no diff p-value = 0.2779
prop.test(c(JTT_sb_10000_b1 ,JTT_sb_10000_cud),c(3000,3000)) # sign diff p-value = 0.0214
prop.test(c(JTT_sb_10000_b1 ,JTT_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_sb_10000_b1 ,JTT_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_sb_10000_b2 ,JTT_sb_10000_b3),c(3000,3000)) # no diff p-value = 0.6421
prop.test(c(JTT_sb_10000_b2 ,JTT_sb_10000_b10),c(3000,3000)) # sign diff p-value = 1.25e-05
prop.test(c(JTT_sb_10000_b2 ,JTT_sb_10000_u),c(3000,3000)) # no diff p-value = 0.1276
prop.test(c(JTT_sb_10000_b2 ,JTT_sb_10000_cu),c(3000,3000)) # sign diff p-value = 0.03006
prop.test(c(JTT_sb_10000_b2 ,JTT_sb_10000_cud),c(3000,3000)) # sign diff p-value = 0.0007128
prop.test(c(JTT_sb_10000_b2 ,JTT_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_sb_10000_b2 ,JTT_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_sb_10000_b3 ,JTT_sb_10000_b10),c(3000,3000)) # sign diff p-value = 0.0001049
prop.test(c(JTT_sb_10000_b3 ,JTT_sb_10000_u),c(3000,3000)) # sign diff p-value = 0.0001049
prop.test(c(JTT_sb_10000_b3 ,JTT_sb_10000_cu),c(3000,3000)) # almost diff p-value = 0.09317
prop.test(c(JTT_sb_10000_b3 ,JTT_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.003799
prop.test(c(JTT_sb_10000_b3 ,JTT_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_sb_10000_b3 ,JTT_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_sb_10000_b10 ,JTT_sb_10000_u),c(3000,3000)) # sign diff p-value = 0.004787
prop.test(c(JTT_sb_10000_b10 ,JTT_sb_10000_cu),c(3000,3000)) # sign diff  p-value = 0.02961
prop.test(c(JTT_sb_10000_b10 ,JTT_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.3373
prop.test(c(JTT_sb_10000_b10 ,JTT_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_sb_10000_b10 ,JTT_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_sb_10000_u ,JTT_sb_10000_cu),c(3000,3000)) # no diff p-value = 0.5351
prop.test(c(JTT_sb_10000_u ,JTT_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.06634
prop.test(c(JTT_sb_10000_u ,JTT_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_sb_10000_u ,JTT_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_sb_10000_cu ,JTT_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.234
prop.test(c(JTT_sb_10000_cu ,JTT_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_sb_10000_cu ,JTT_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_sb_10000_cud ,JTT_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(JTT_sb_10000_cud ,JTT_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(JTT_sb_10000_NJ ,JTT_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

# LG sb 10000 succ. rate

LG_sb_10000_b1 = 503 + 507 + 446 
LG_sb_10000_b1 # 1456
LG_sb_10000_b2 = 451 + 491 + 530 
LG_sb_10000_b2 # 1472
LG_sb_10000_b3 = 509 + 502 + 460 
LG_sb_10000_b3 # 1471
LG_sb_10000_b10 = 460 + 416 + 458
LG_sb_10000_b10 # 1334
LG_sb_10000_u = 430 + 513 + 488
LG_sb_10000_u # 1431
LG_sb_10000_cu = 369 + 524 + 469
LG_sb_10000_cu # 1362
LG_sb_10000_cud = 467 + 498 + 415 
LG_sb_10000_cud # 1380
LG_sb_10000_NJ = 631 + 639 + 592 
LG_sb_10000_NJ # 1862
LG_sb_10000_ML = 854 + 883 + 883 
LG_sb_10000_ML #  2620

LG_sb_10000_b1/30 #  48.53%
LG_sb_10000_b2/30 # 49.07%
LG_sb_10000_b3/30 # 49.03%
LG_sb_10000_b10/30 # 44.47%
LG_sb_10000_u/30 # 47.7%
LG_sb_10000_cu/30 # 45.4%
LG_sb_10000_cud/30 # 46.0%
LG_sb_10000_NJ/30 # 62.07%
LG_sb_10000_ML/30 # 87.33%

# LG_sb 10000 binomial test 

prop.test(c(LG_sb_10000_b1 ,LG_sb_10000_b2),c(3000,3000)) # no diff p-value = 0.6985
prop.test(c(LG_sb_10000_b1 ,LG_sb_10000_b3),c(3000,3000)) # no diff p-value = 0.7177
prop.test(c(LG_sb_10000_b1 ,LG_sb_10000_b10),c(3000,3000)) # sign diff p-value = 0.001737
prop.test(c(LG_sb_10000_b1 ,LG_sb_10000_u),c(3000,3000)) # no diff p-value = 0.5352
prop.test(c(LG_sb_10000_b1 ,LG_sb_10000_cu),c(3000,3000)) # sign diff  p-value = 0.01614
prop.test(c(LG_sb_10000_b1 ,LG_sb_10000_cud),c(3000,3000)) # almost diff p-value = 0.05245
prop.test(c(LG_sb_10000_b1 ,LG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_sb_10000_b1 ,LG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_sb_10000_b2 ,LG_sb_10000_b3),c(3000,3000)) # no diff p-value = 1
prop.test(c(LG_sb_10000_b2 ,LG_sb_10000_b10),c(3000,3000)) # sign diff p-value = 0.000393
prop.test(c(LG_sb_10000_b2 ,LG_sb_10000_u),c(3000,3000)) # no diff p-value = 0.3014
prop.test(c(LG_sb_10000_b2 ,LG_sb_10000_cu),c(3000,3000)) # sign diff p-value = 0.004822
prop.test(c(LG_sb_10000_b2 ,LG_sb_10000_cud),c(3000,3000)) # sign diff p-value = 0.01865
prop.test(c(LG_sb_10000_b2 ,LG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_sb_10000_b2 ,LG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_sb_10000_b3 ,LG_sb_10000_b10),c(3000,3000)) # sign diff p-value = 0.0004333
prop.test(c(LG_sb_10000_b3 ,LG_sb_10000_u),c(3000,3000)) # no diff p-value = 0.3137
prop.test(c(LG_sb_10000_b3 ,LG_sb_10000_cu),c(3000,3000)) # sign diff p-value = 0.005224
prop.test(c(LG_sb_10000_b3 ,LG_sb_10000_cud),c(3000,3000)) # sign diff p-value = 0.01998
prop.test(c(LG_sb_10000_b3 ,LG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_sb_10000_b3 ,LG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_sb_10000_b10 ,LG_sb_10000_u),c(3000,3000)) # sign diff p-value = 0.01291
prop.test(c(LG_sb_10000_b10 ,LG_sb_10000_cu),c(3000,3000)) # no diff p-value = 0.4835
prop.test(c(LG_sb_10000_b10 ,LG_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.2431
prop.test(c(LG_sb_10000_b10 ,LG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_sb_10000_b10 ,LG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_sb_10000_u ,LG_sb_10000_cu),c(3000,3000)) # almost diff p-value = 0.07842
prop.test(c(LG_sb_10000_u ,LG_sb_10000_cud),c(3000,3000)) # no diff p-value = 0.1958
prop.test(c(LG_sb_10000_u ,LG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_sb_10000_u ,LG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_sb_10000_cu ,LG_sb_10000_cud),c(3000,3000)) # no diff  p-value = 0.6595
prop.test(c(LG_sb_10000_cu ,LG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_sb_10000_cu ,LG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_sb_10000_cud ,LG_sb_10000_NJ),c(3000,3000)) # sign diff p-value < 2.2e-16
prop.test(c(LG_sb_10000_cud ,LG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16

prop.test(c(LG_sb_10000_NJ ,LG_sb_10000_ML),c(3000,3000)) # sign diff p-value < 2.2e-16
