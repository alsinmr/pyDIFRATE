try:
	from chimerax.core.commands import run as rc

	import os
	import numpy as np

	di=np.array([0,4,17,18,19,20,21,24,26,30,\
		38,39,40,41,42,57,58,59,60,61,\
		64,65,66,67,68,74,75,76,77,78,\
		89,90,91,92,93,99,100,101,102,103,\
		111,112,113,114,115,126,127,128,129,130,\
		140,141,142,143,144,157,158,159,160,161,\
		171,172,173,174,175,187,188,189,190,191,\
		202,203,204,205,206,217,218,219,220,221,\
		234,235,236,237,238,256,257,258,259,260,\
		272,273,274,275,276,287,288,289,290,291,\
		304,305,306,307,308,325,326,327,328,329,\
		332,334,338,346,347,348,349,350,363,364,\
		365,366,367,377,378,379,380,381,391,393,\
		397,405,406,407,408,409,429,430,431,432,\
		433,436,437,438,439,440,451,452,453,454,\
		455,470,471,472,473,474,486,488,492,500,\
		501,502,503,504,512,514,518,526,527,528,\
		529,530,541,543,547,555,556,557,558,559,\
		570,571,572,573,574,589,590,591,592,593,\
		608,609,610,611,612,620,621,622,623,624,\
		631,632,633,634,635,645,646,647,648,649,\
		667,668,669,670,671,686,687,688,689,690,\
		705,706,707,708,709,720,721,722,723,724,\
		736,737,738,739,740,753,754,755,756,757,\
		769,770,771,772,773,785,786,787,788,789,\
		804,805,806,807,808,823,824,825,826,827,\
		842,843,844,845,846,852,853,854,855,856,\
		873,874,875,876,877,883,884,885,886,887,\
		894,895,896,897,898,913,914,915,916,917,\
		932,933,934,935,936,951,952,953,954,955,\
		970,971,972,973,974,977,978,979,980,981,\
		993,994,995,996,997,1012,1013,1014,1015,1016,\
		1019,1020,1021,1022,1023,1033,1034,1035,1036,1037,\
		1044,1045,1046,1047,1048,1063,1064,1065,1066,1067,\
		1079,1080,1081,1082,1083,1098,1099,1100,1101,1102,\
		1115,1116,1117,1118,1119,1131,1132,1133,1134,1135,\
		1147,1148,1149,1150,1151,1166,1167,1168,1169,1170,\
		1188,1189,1190,1191,1192,1208,1209,1210,1211,1212,\
		1230,1231,1232,1233,1234,1241,1242,1243,1244,1245,\
		1258,1259,1260,1261,1262,1282,1283,1284,1285,1286,\
		1296,1297,1298,1299,1300,1312,1313,1314,1315,1316,\
		1326,1327,1328,1329,1330,1340,1341,1342,1343,1344,\
		1360,1361,1362,1363,1364,1380,1381,1382,1383,1384,\
		1399,1400,1401,1402,1403,1409,1410,1411,1412,1413,\
		1423,1424,1425,1426,1427,1442,1443,1444,1445,1446,\
		1452,1453,1454,1455,1456,1468,1469,1470,1471,1472,\
		1478,1479,1480,1481,1482,1490,1491,1492,1493,1494,\
		1509,1510,1511,1512,1513,1528,1529,1530,1531,1532,\
		1544,1545,1546,1547,1548,1558,1559,1560,1561,1562,\
		1572,1573,1574,1575,1576,1591,1592,1593,1594,1595,\
		1602,1603,1604,1605,1606,1621,1623,1627,1635,1636,\
		1637,1638,1639,1655,1656,1657,1658,1659,1669,1670,\
		1671,1672,1673,1688,1689,1690,1691,1692,1702,1703,\
		1704,1705,1706,1723,1724,1725,1726,1727,1737,1738,\
		1739,1740,1741,1756,1757,1758,1759,1760,1773,1774,\
		1775,1776,1777,1780,1781,1782,1783,1784,1795,1796,\
		1797,1798,1799,1819,1820,1821,1822,1823,1841,1842,\
		1843,1844,1845,1858,1859,1860,1861,1862,1865,1867,\
		1871,1879,1880,1881,1882,1883,1895,1896,1897,1898,\
		1899,1914,1915,1916,1917,1918,1924,1925,1926,1927,\
		1928,1941,1942,1943,1944,1945,1960,1961,1962,1963,\
		1964,1976,1978,1982,1990,1991,1992,1993,1994,2011,\
		2012,2013,2014,2015,2021,2022,2023,2024,2025,2038,\
		2039,2040,2041,2042,2045,2046,2047,2048,2049,2064,\
		2065,2066,2067,2068,2074,2075,2076,2077,2078,2090,\
		2091,2092,2093,2094,2107,2108,2109,2110,2111,2123,\
		2124,2125,2126,2127,2134,2135,2136,2137,2138,2148,\
		2149,2150,2151,2152,2167,2168,2169,2170,2171,2181,\
		2182,2183,2184,2185,2200,2201,2202,2203,2204,2214,\
		2215,2216,2217,2218,2230,2231,2232,2233,2234,2249,\
		2250,2251,2252,2253,2259,2260,2261,2262,2263,2278,\
		2279,2280,2281,2282,2290,2291,2292,2293,2294,2314,\
		2315,2316,2317,2318,2331,2332,2333,2334,2335,2355,\
		2356,2357,2358,2359,2366,2367,2368,2369,2370,2385,\
		2386,2387,2388,2389,2401,2402,2403,2404,2405,2422,\
		2423,2424,2425,2426,2439,2440,2441,2442,2443,2458,\
		2459,2460,2461,2462,2473,2474,2475,2476,2477,2484,\
		2485,2486,2487,2488,2506,2507,2508,2509,2510,2525,\
		2526,2527,2528,2529,2536,2537,2538,2539,2540,2558,\
		2559,2560,2561,2562,2582,2583,2584,2585,2586,2601,\
		2602,2603,2604,2605,2612,2613,2614,2615,2616,2632,\
		2633,2634,2635,2636,2651,2652,2653,2654,2655,2670,\
		2671,2672,2673,2674,2689,2690,2691,2692,2693,2696,\
		2697,2698,2699,2700,2715,2716,2717,2718,2719,2725,\
		2726,2727,2728,2729,2749,2750,2751,2752,2753,2756,\
		2757,2758,2759,2760,2775,2776,2777,2778,2779,2786,\
		2787,2788,2789,2790,2796,2797,2798,2799,2800,2815,\
		2816,2817,2818,2819,2834,2835,2836,2837,2838,2844,\
		2845,2846,2847,2848,2855,2857,2861,2869,2870,2871,\
		2872,2873,2888,2889,2890,2891,2892,2898,2899,2900,\
		2901,2902,2917,2918,2919,2920,2921,2937,2938,2939,\
		2940,2941,2961,2962,2963,2964,2965,2976,2977,2978,\
		2979,2980,2997,2998,2999,3000,3001,3008,3009,3010,\
		3011,3012,3027,3028,3029,3030,3031,3046,3047,3048,\
		3049,3050,3061,3062,3063,3064,3065,3080,3081,3082,\
		3083,3084,3099,3101,3105,3113,3114,3115,3116,3117,\
		3125,3126,3127,3128,3129,3145,3146,3147,3148,3149,\
		3160,3161,3162,3163,3164,3179,3180,3181,3182,3183,\
		3195,3196,3197,3198,3199,3205,3206,3207,3208,3209,\
		3215,3216,3217,3218,3219,3229,3230,3231,3232,3233,\
		3244,3245,3246,3247,3248,3266,3267,3268,3269,3270,\
		3290,3292,3296,3304,3305,3306,3307,3308,3311,3312,\
		3313,3314,3315,3326,3327,3328,3329,3330,3341,3342,\
		3343,3344,3345,3363,3364,3365,3366,3367,3374,3375,\
		3376,3377,3378,3393,3394,3395,3396,3397,3414,3415,\
		3416,3417,3418,3421,3422,3423,3424,3425,3435,3436,\
		3437,3438,3439,3451,3452,3453,3454,3455,3472,3473,\
		3474,3475,3476,3483,3484,3485,3486,3487,3502,3503,\
		3504,3505,3506,3513,3514,3515,3516,3517,3524,3525,\
		3526,3527,3528,3543,3544,3545,3546,3547,3562,3563,\
		3564,3565,3566,3581,3582,3583,3584,3585,3600,3601,\
		3602,3603,3604,3621,3622,3623,3624,3625,3637,3638,\
		3639,3640,3641,3656,3658,3662,3670,3671,3672,3673,\
		3674,3689,3690,3691,3692,3693,3696,3697,3698,3699,\
		3700,3715,3716,3717,3718,3719,3734,3735,3736,3737,\
		3738,3745,3746,3747,3748,3749,3765,3766,3767,3768,\
		3769,3776,3777,3778,3779,3780,3797,3798,3799,3800,\
		3801,3811,3812,3813,3814,3815,3835,3836,3837,3838,\
		3839,3854,3855,3856,3857,3858,3878,3879,3880,3881,\
		3882,3889,3890,3891,3892,3893,3911,3912,3913,3914,\
		3915,3930,3931,3932,3933,3934,3952,3953,3954,3955,\
		3956,3966,3967,3968,3969,3970,3983,3984,3985,3986,\
		3987,3999,4000,4001,4002,4003,4010,4012,4016,4024,\
		4025,4026,4027,4028,4031,4032,4033,4034,4035,4041,\
		4042,4043,4044,4045,4051,4052,4053,4054,4055,4065,\
		4066,4067,4068,4069,4077,4078,4079,4080,4081,4094,\
		4095,4096,4097,4098,4115,4116,4117,4118,4119,4132,\
		4133,4134,4135,4136,4149,4150,4151,4152,4153,4173,\
		4174,4175,4176,4177,4197,4198,4199,4200,4201,4214,\
		4215,4216,4217,4218,4236,4237,4238,4239,4240,4250,\
		4251,4252,4253,4254,4264,4265,4266,4267,4268,4286,\
		4287,4288,4289,4290,4303,4304,4305,4306,4307,4322,\
		4323,4324,4325,4326,4338,4339,4340,4341,4342,4348,\
		4349,4350,4351,4352,4364,4365,4366,4367,4368,4380,\
		4381,4382,4383,4384,4396,4397,4398,4399,4400,4412,\
		4413,4414,4415,4416,4432,4433,4434,4435,4436,4442,\
		4443,4444,4445,4446,4458,4459,4460,4461,4462,4469,\
		4470,4471,4472,4473,4493,4494,4495,4496,4497,4512,\
		4514,4518,4526,4527,4528,4529,4530,4545,4546,4547,\
		4548,4549,4562,4563,4564,4565,4566,4572,4573,4574,\
		4575,4576,4592,4593,4594,4595,4596,4609,4610,4611,\
		4612,4613,4628,4629,4630,4631,4632,4638,4639,4640,\
		4641,4642,4654,4655,4656,4657,4658,4666,4667,4668,\
		4669,4670,4685,4686,4687,4688,4689,4697,4698,4699,\
		4700,4701,4708,4709,4710,4711,4712,4725,4726,4727,\
		4728,4729,4741,4742,4743,4744,4745,4760,4761,4762,\
		4763,4764,4772,4773,4774,4775,4776,4791,4792,4793,\
		4794,4795,4813,4814,4815,4816,4817,4828,4829,4830,\
		4831,4832,4849,4850,4851,4852,4853,4871,4872,4873,\
		4874,4875,4890,4891,4892,4893,4894,4909,4910,4911,\
		4912,4913,4929,4930,4931,4932,4933,4943,4944,4945,\
		4946,4947,4959,4960,4961,4962,4963,4979,4980,4981,\
		4982,4983,4996,4997,4998,4999,5000,5015,5016,5017,\
		5018,5019,5034,5035,5036,5037,5038,5044,5045,5046,\
		5047,5048,5061,5062,5063,5064,5065,5071,5072,5073,\
		5074,5075,5082,5083,5084,5085,5086,5096,5097,5098,\
		5099,5100,5116,5117,5118,5119,5120,5126,5127,5128,\
		5129,5130,5140,5142,5146,5154,5155,5156,5157,5158,\
		5173,5174,5175,5176,5177,5192,5193,5194,5195,5196,\
		5213,5214,5215,5216,5217,5220,5221,5222,5223,5224,\
		5244,5245,5246,5247,5248,5261,5262,5263,5264,5265,\
		5275,5276,5277,5278,5279,5286,5287,5288,5289,5290,\
		5300,5301,5302,5303,5304,5321,5322,5323,5324,5325,\
		5345,5346,5347,5348,5349,5367,5368,5369,5370,5371,\
		5377,5378,5379,5380,5381,5397,5398,5399,5400,5401,\
		5416,5417,5418,5419,5420,5427,5428,5429,5430,5431,\
		5437,5438,5439,5440,5441,5457,5458,5459,5460,5461,\
		5481,5482,5483,5484,5485,5491,5492,5493,5494,5495,\
		5506,5507,5508,5509,5510,5523,5524,5525,5526,5527,\
		5547,5548,5549,5550,5551,5566,5567,5568,5569,5570,\
		5578,5579,5580,5581,5582,5588,5589,5590,5591,5592,\
		5607,5608,5609,5610,5611,5624,5625,5626,5627,5628,\
		5635,5636,5637,5638,5639,5650,5651,5652,5653,5654,\
		5666,5667,5668,5669,5670,5677,5678,5679,5680,5681,\
		5693,5694,5695,5696,5697,5707,5708,5709,5710,5711,\
		5727,5728,5729,5730,5731,5749,5750,5751,5752,5753,\
		5759,5760,5761,5762,5763,5781,5782,5783,5784,5785,\
		5803,5804,5805,5806,5807,5817,5818,5819,5820,5821,\
		5836,5837,5838,5839,5840,5851,5852,5853,5854,5855,\
		5867,5868,5869,5870,5871,5891,5892,5893,5894,5895,\
		5913,5914,5915,5916,5917,5927,5928,5929,5930,5931,\
		5938,5939,5940,5941,5942,5945,5947,5951,5959,5960,\
		5961,5962,5963,5973,5974,5975,5976,5977,5985,5986,\
		5987,5988,5989,5996,5997,5998,5999,6000,6016,6017,\
		6018,6019,6020,6030,6031,6032,6033,6034,6045,6046,\
		6047,6048,6049,6055,6056,6057,6058,6059,6069,6070,\
		6071,6072,6073,6083,6084,6085,6086,6087,6099,6100,\
		6101,6102,6103,6118,6119,6120,6121,6122,6133,6134,\
		6135,6136,6137,6150,6151,6152,6153,6154,6167,6168,\
		6169,6170,6171,6184,6185,6186,6187,6188,6201,6202,\
		6203,6204,6205,6218,6219,6220,6221,6222,6235,6236,\
		6237,6238,6239,6252,6253,6254,6255,6256,6269]).astype("uint32")