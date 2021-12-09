try:
	from chimerax.core.commands import run as rc

	import os
	import numpy as np

	di=np.array([6,8,24,26,28,35,37,39,50,52,\
		54,72,74,76,86,88,90,96,98,100,\
		113,115,117,132,134,136,151,153,155,161,\
		163,165,181,183,185,196,198,200,210,212,\
		214,222,224,226,234,236,238,244,246,248,\
		261,263,265,280,282,286,294,296,298,313,\
		314,315,316,317,323,324,325,326,327,340,\
		341,342,343,344,359,360,361,362,363,379,\
		380,381,382,383,393,394,395,396,397,412,\
		413,414,415,416,422,423,424,425,426,441,\
		442,443,444,445,451,452,453,454,455,472,\
		474,476,479,480,481,482,483,489,490,491,\
		492,493,505,506,507,508,509,524,525,526,\
		527,528,543,544,545,546,547,562,563,564,\
		565,566,569,570,571,572,573,585,586,587,\
		588,589,596,597,598,599,600,603,604,605,\
		606,607,617,618,619,620,621,636,637,638,\
		639,640,646,647,648,649,650,665,666,667,\
		668,669,684,685,686,687,688,703,704,705,\
		706,707,722,723,724,725,726,741,742,743,\
		744,745,760,761,762,763,764,782,784,786,\
		799,801,803,821,823,825,836,838,840,853,\
		855,857,877,879,881,891,893,895,907,908,\
		909,910,911,921,922,923,924,925,935,937,\
		939,954,955,956,957,958,973,974,975,976,\
		977,992,993,994,995,996,1008,1009,1010,1011,\
		1012,1022,1023,1024,1025,1026,1041,1042,1043,1044,\
		1045,1052,1053,1054,1055,1056,1072,1073,1074,1075,\
		1076,1083,1084,1085,1086,1087,1095,1096,1097,1098,\
		1099,1114,1115,1116,1117,1118,1133,1134,1135,1136,\
		1137,1149,1150,1151,1152,1153,1159,1160,1161,1162,\
		1163,1178,1180,1182,1195,1197,1199,1206,1208,1210,\
		1225,1227,1231,1239,1240,1241,1242,1243,1259,1260,\
		1261,1262,1263,1273,1274,1275,1276,1277,1293,1295,\
		1297,1309,1310,1311,1312,1313,1330,1331,1332,1333,\
		1334,1344,1345,1346,1347,1348,1363,1364,1365,1366,\
		1367,1380,1382,1384,1392,1394,1396,1409,1411,1413,\
		1433,1435,1437,1449,1451,1453,1469,1471,1473,1476,\
		1478,1480,1491,1492,1493,1494,1495,1501,1502,1503,\
		1504,1505,1518,1519,1520,1521,1522,1528,1529,1530,\
		1531,1532,1550,1551,1552,1553,1554,1569,1570,1571,\
		1572,1573,1583,1585,1589,1597,1599,1601,1617,1618,\
		1619,1620,1621,1633,1634,1635,1636,1637,1650,1651,\
		1652,1653,1654,1661,1662,1663,1664,1665,1677,1678,\
		1679,1680,1681,1688,1689,1690,1691,1692,1707,1708,\
		1709,1710,1711,1721,1722,1723,1724,1725,1737,1739,\
		1741,1748,1749,1750,1751,1752,1767,1768,1769,1770,\
		1771,1787,1788,1789,1790,1791,1798,1799,1800,1801,\
		1802,1817,1819,1821,1833,1834,1835,1836,1837,1852,\
		1853,1854,1855,1856,1871,1872,1873,1874,1875,1881,\
		1882,1883,1884,1885,1897,1898,1899,1900,1901,1912,\
		1913,1914,1915,1916,1936,1938,1940,1953,1955,1957,\
		1970,1972,1974,1989,1991,1993,2008,2010,2012,2027,\
		2029,2031,2041,2043,2047,2055,2057,2059,2079,2081,\
		2083,2086,2088,2090,2110,2112,2114,2134,2136,2140,\
		2148,2150,2152,2162,2164,2166,2176,2177,2178,2179,\
		2180,2200,2202,2204,2217,2218,2219,2220,2221,2227,\
		2228,2229,2230,2231,2248,2249,2250,2251,2252,2264,\
		2265,2266,2267,2268,2271,2272,2273,2274,2275,2290,\
		2291,2292,2293,2294,2300,2301,2302,2303,2304,2316,\
		2317,2318,2319,2320,2335,2336,2337,2338,2339,2359,\
		2360,2361,2362,2363,2375,2377,2379,2394,2395,2396,\
		2397,2398,2404,2405,2406,2407,2408,2420,2421,2422,\
		2423,2424,2430,2431,2432,2433,2434,2441,2442,2443,\
		2444,2445,2452,2454,2456,2471,2473,2477,2485,2486,\
		2487,2488,2489,2505,2506,2507,2508,2509,2524,2525,\
		2526,2527,2528,2543,2544,2545,2546,2547,2564,2566,\
		2568,2581,2583,2585,2597,2599,2601,2614,2616,2618,\
		2628,2630,2632,2640,2642,2644,2655,2657,2661,2669,\
		2671,2673,2689,2691,2693,2706,2708,2710,2720,2722,\
		2724,2736,2738,2740,2750,2752,2754,2769,2771,2773,\
		2781,2783,2785,2791,2793,2795,2812,2814,2816,2834,\
		2836,2838,2846,2848,2850,2868,2870,2872,2889,2891,\
		2893,2905,2907,2909,2915,2917,2919,2935,2937,2939,\
		2947,2949,2951,2964,2966,2968,2984,2986,2990,2998,\
		3000,3002,3009,3011,3013,3021,3022,3023,3024,3025,\
		3032,3033,3034,3035,3036,3049,3050,3051,3052,3053,\
		3073,3074,3075,3076,3077,3092,3093,3094,3095,3096,\
		3103,3104,3105,3106,3107,3124,3125,3126,3127,3128,\
		3138,3139,3140,3141,3142,3152,3153,3154,3155,3156,\
		3171,3172,3173,3174,3175,3190,3191,3192,3193,3194,\
		3209,3210,3211,3212,3213,3225,3226,3227,3228,3229,\
		3244,3246,3248,3261,3263,3265,3282,3284,3286,3302,\
		3304,3306,3309,3311,3315,3323,3324,3325,3326,3327,\
		3342,3343,3344,3345,3346,3353,3354,3355,3356,3357,\
		3373,3374,3375,3376,3377,3392,3393,3394,3395,3396,\
		3412,3413,3414,3415,3416,3431,3432,3433,3434,3435,\
		3442,3443,3444,3445,3446,3463,3464,3465,3466,3467,\
		3483,3484,3485,3486,3487,3505,3507,3509,3524,3526,\
		3528,3545,3547,3549,3564,3566,3568,3588,3590,3592,\
		3607,3609,3611,3629,3631,3633,3653,3655,3657,3677,\
		3679,3681,3691,3693,3695,3705,3707,3709,3722,3724,\
		3726,3739,3741,3743,3751,3753,3755,3773,3775,3777,\
		3790,3792,3794,3814,3816,3818,3826,3828,3830,3840,\
		3842,3844,3862,3864,3866,3883,3885,3887,3907,3909,\
		3911,3918,3920,3922,3929,3931,3933,3944,3946,3948,\
		3958,3960,3962,3980,3982,3984,4004,4006,4008,4023,\
		4025,4027,4037,4039,4041,4056,4057,4058,4059,4060,\
		4073,4074,4075,4076,4077,4092,4093,4094,4095,4096,\
		4111,4112,4113,4114,4115,4122,4123,4124,4125,4126,\
		4141,4142,4143,4144,4145,4157,4158,4159,4160,4161,\
		4173,4174,4175,4176,4177,4183,4184,4185,4186,4187,\
		4203,4204,4205,4206,4207,4213,4214,4215,4216,4217,\
		4229,4230,4231,4232,4233,4240,4241,4242,4243,4244,\
		4264,4265,4266,4267,4268,4283,4285,4289,4297,4298,\
		4299,4300,4301,4316,4317,4318,4319,4320,4330,4331,\
		4332,4333,4334,4349,4350,4351,4352,4353,4369,4370,\
		4371,4372,4373,4383,4384,4385,4386,4387,4397,4398,\
		4399,4400,4401,4413,4414,4415,4416,4417,4433,4434,\
		4435,4436,4437,4445,4446,4447,4448,4449,4469,4471,\
		4473,4483,4485,4487,4500,4502,4504,4517,4519,4521,\
		4536,4538,4540,4555,4557,4559,4565,4567,4569,4579,\
		4581,4583,4589,4590,4591,4592,4593,4603,4604,4605,\
		4606,4607,4620,4622,4624,4634,4635,4636,4637,4638,\
		4653,4654,4655,4656,4657,4672,4673,4674,4675,4676,\
		4692,4693,4694,4695,4696,4711,4712,4713,4714,4715,\
		4730,4731,4732,4733,4734,4741,4742,4743,4744,4745,\
		4758,4759,4760,4761,4762,4777,4778,4779,4780,4781,\
		4791,4792,4793,4794,4795,4801,4802,4803,4804,4805,\
		4818,4819,4820,4821,4822,4837,4838,4839,4840,4841,\
		4848,4849,4850,4851,4852,4862,4863,4864,4865,4866,\
		4873,4875,4877,4889,4890,4891,4892,4893,4903,4905,\
		4909,4917,4918,4919,4920,4921,4936,4937,4938,4939,\
		4940,4956,4957,4958,4959,4960,4977,4978,4979,4980,\
		4981,4984,4985,4986,4987,4988,5004,5005,5006,5007,\
		5008,5023,5025,5027,5037,5039,5041,5059,5061,5063,\
		5073,5075,5077,5093,5095,5097,5110,5112,5114,5134,\
		5136,5138,5146,5148,5150,5165,5167,5169,5182,5184,\
		5186,5202,5204,5206,5222,5224,5226,5242,5244,5246,\
		5256,5258,5260,5276,5278,5280,5288]).astype("uint32")

	ids=np.array([313,314,315,316,323,324,325,326,340,341,\
		342,343,359,360,361,362,379,380,381,382,\
		393,394,395,396,412,413,414,415,422,423,\
		424,425,441,442,443,444,451,452,453,454,\
		479,480,481,482,489,490,491,492,505,506,\
		507,508,524,525,526,527,543,544,545,546,\
		562,563,564,565,569,570,571,572,585,586,\
		587,588,596,597,598,599,603,604,605,606,\
		617,618,619,620,636,637,638,639,646,647,\
		648,649,665,666,667,668,684,685,686,687,\
		703,704,705,706,722,723,724,725,741,742,\
		743,744,760,761,762,763,907,908,909,910,\
		921,922,923,924,954,955,956,957,973,974,\
		975,976,992,993,994,995,1008,1009,1010,1011,\
		1022,1023,1024,1025,1041,1042,1043,1044,1052,1053,\
		1054,1055,1072,1073,1074,1075,1083,1084,1085,1086,\
		1095,1096,1097,1098,1114,1115,1116,1117,1133,1134,\
		1135,1136,1149,1150,1151,1152,1159,1160,1161,1162,\
		1239,1240,1241,1242,1259,1260,1261,1262,1273,1274,\
		1275,1276,1309,1310,1311,1312,1330,1331,1332,1333,\
		1344,1345,1346,1347,1363,1364,1365,1366,1491,1492,\
		1493,1494,1501,1502,1503,1504,1518,1519,1520,1521,\
		1528,1529,1530,1531,1550,1551,1552,1553,1569,1570,\
		1571,1572,1617,1618,1619,1620,1633,1634,1635,1636,\
		1650,1651,1652,1653,1661,1662,1663,1664,1677,1678,\
		1679,1680,1688,1689,1690,1691,1707,1708,1709,1710,\
		1721,1722,1723,1724,1748,1749,1750,1751,1767,1768,\
		1769,1770,1787,1788,1789,1790,1798,1799,1800,1801,\
		1833,1834,1835,1836,1852,1853,1854,1855,1871,1872,\
		1873,1874,1881,1882,1883,1884,1897,1898,1899,1900,\
		1912,1913,1914,1915,2176,2177,2178,2179,2217,2218,\
		2219,2220,2227,2228,2229,2230,2248,2249,2250,2251,\
		2264,2265,2266,2267,2271,2272,2273,2274,2290,2291,\
		2292,2293,2300,2301,2302,2303,2316,2317,2318,2319,\
		2335,2336,2337,2338,2359,2360,2361,2362,2394,2395,\
		2396,2397,2404,2405,2406,2407,2420,2421,2422,2423,\
		2430,2431,2432,2433,2441,2442,2443,2444,2485,2486,\
		2487,2488,2505,2506,2507,2508,2524,2525,2526,2527,\
		2543,2544,2545,2546,3021,3022,3023,3024,3032,3033,\
		3034,3035,3049,3050,3051,3052,3073,3074,3075,3076,\
		3092,3093,3094,3095,3103,3104,3105,3106,3124,3125,\
		3126,3127,3138,3139,3140,3141,3152,3153,3154,3155,\
		3171,3172,3173,3174,3190,3191,3192,3193,3209,3210,\
		3211,3212,3225,3226,3227,3228,3323,3324,3325,3326,\
		3342,3343,3344,3345,3353,3354,3355,3356,3373,3374,\
		3375,3376,3392,3393,3394,3395,3412,3413,3414,3415,\
		3431,3432,3433,3434,3442,3443,3444,3445,3463,3464,\
		3465,3466,3483,3484,3485,3486,4056,4057,4058,4059,\
		4073,4074,4075,4076,4092,4093,4094,4095,4111,4112,\
		4113,4114,4122,4123,4124,4125,4141,4142,4143,4144,\
		4157,4158,4159,4160,4173,4174,4175,4176,4183,4184,\
		4185,4186,4203,4204,4205,4206,4213,4214,4215,4216,\
		4229,4230,4231,4232,4240,4241,4242,4243,4264,4265,\
		4266,4267,4297,4298,4299,4300,4316,4317,4318,4319,\
		4330,4331,4332,4333,4349,4350,4351,4352,4369,4370,\
		4371,4372,4383,4384,4385,4386,4397,4398,4399,4400,\
		4413,4414,4415,4416,4433,4434,4435,4436,4445,4446,\
		4447,4448,4589,4590,4591,4592,4603,4604,4605,4606,\
		4634,4635,4636,4637,4653,4654,4655,4656,4672,4673,\
		4674,4675,4692,4693,4694,4695,4711,4712,4713,4714,\
		4730,4731,4732,4733,4741,4742,4743,4744,4758,4759,\
		4760,4761,4777,4778,4779,4780,4791,4792,4793,4794,\
		4801,4802,4803,4804,4818,4819,4820,4821,4837,4838,\
		4839,4840,4848,4849,4850,4851,4862,4863,4864,4865,\
		4889,4890,4891,4892,4917,4918,4919,4920,4936,4937,\
		4938,4939,4956,4957,4958,4959,4977,4978,4979,4980,\
		4984,4985,4986,4987,5004,5005,5006,5007]).astype("uint32")

	r=np.array([1.052190,1.052190,1.052190,1.052190,1.226754,1.226754,1.226754,1.226754,1.171208,1.171208,\
		1.171208,1.171208,1.353635,1.353635,1.353635,1.353635,1.371079,1.371079,1.371079,1.371079,\
		1.192812,1.192812,1.192812,1.192812,1.232747,1.232747,1.232747,1.232747,1.058128,1.058128,\
		1.058128,1.058128,1.333054,1.333054,1.333054,1.333054,1.143108,1.143108,1.143108,1.143108,\
		1.062685,1.062685,1.062685,1.062685,1.286806,1.286806,1.286806,1.286806,1.218035,1.218035,\
		1.218035,1.218035,1.069501,1.069501,1.069501,1.069501,1.391118,1.391118,1.391118,1.391118,\
		1.204427,1.204427,1.204427,1.204427,1.293042,1.293042,1.293042,1.293042,1.069334,1.069334,\
		1.069334,1.069334,1.346911,1.346911,1.346911,1.346911,1.162438,1.162438,1.162438,1.162438,\
		1.061606,1.061606,1.061606,1.061606,1.064623,1.064623,1.064623,1.064623,1.299665,1.299665,\
		1.299665,1.299665,1.213491,1.213491,1.213491,1.213491,1.064062,1.064062,1.064062,1.064062,\
		1.405841,1.405841,1.405841,1.405841,1.209760,1.209760,1.209760,1.209760,1.312921,1.312921,\
		1.312921,1.312921,1.071130,1.071130,1.071130,1.071130,1.024713,1.024713,1.024713,1.024713,\
		1.084251,1.084251,1.084251,1.084251,1.125015,1.125015,1.125015,1.125015,0.941138,0.941138,\
		0.941138,0.941138,1.144340,1.144340,1.144340,1.144340,1.180328,1.180328,1.180328,1.180328,\
		0.982311,0.982311,0.982311,0.982311,1.002548,1.002548,1.002548,1.002548,1.194957,1.194957,\
		1.194957,1.194957,1.145944,1.145944,1.145944,1.145944,0.928880,0.928880,0.928880,0.928880,\
		1.102544,1.102544,1.102544,1.102544,1.070549,1.070549,1.070549,1.070549,1.174211,1.174211,\
		1.174211,1.174211,1.018986,1.018986,1.018986,1.018986,1.069017,1.069017,1.069017,1.069017,\
		1.168155,1.168155,1.168155,1.168155,0.963701,0.963701,0.963701,0.963701,1.118913,1.118913,\
		1.118913,1.118913,1.073177,1.073177,1.073177,1.073177,0.968971,0.968971,0.968971,0.968971,\
		1.223263,1.223263,1.223263,1.223263,1.289329,1.289329,1.289329,1.289329,1.524298,1.524298,\
		1.524298,1.524298,1.068357,1.068357,1.068357,1.068357,1.307967,1.307967,1.307967,1.307967,\
		1.540097,1.540097,1.540097,1.540097,1.093222,1.093222,1.093222,1.093222,1.069601,1.069601,\
		1.069601,1.069601,1.198654,1.198654,1.198654,1.198654,0.936704,0.936704,0.936704,0.936704,\
		1.100557,1.100557,1.100557,1.100557,1.293699,1.293699,1.293699,1.293699,1.088567,1.088567,\
		1.088567,1.088567,0.921190,0.921190,0.921190,0.921190,1.179415,1.179415,1.179415,1.179415,\
		1.069559,1.069559,1.069559,1.069559,0.939845,0.939845,0.939845,0.939845,1.186225,1.186225,\
		1.186225,1.186225,1.115175,1.115175,1.115175,1.115175,0.924606,0.924606,0.924606,0.924606,\
		1.241548,1.241548,1.241548,1.241548,1.016811,1.016811,1.016811,1.016811,0.916609,0.916609,\
		0.916609,0.916609,1.168463,1.168463,1.168463,1.168463,1.172956,1.172956,1.172956,1.172956,\
		0.944262,0.944262,0.944262,0.944262,1.560579,1.560579,1.560579,1.560579,0.900000,0.900000,\
		0.900000,0.900000,1.284011,1.284011,1.284011,1.284011,0.964327,0.964327,0.964327,0.964327,\
		0.979763,0.979763,0.979763,0.979763,1.193521,1.193521,1.193521,1.193521,1.119517,1.119517,\
		1.119517,1.119517,0.911570,0.911570,0.911570,0.911570,1.066106,1.066106,1.066106,1.066106,\
		1.208235,1.208235,1.208235,1.208235,1.015712,1.015712,1.015712,1.015712,1.197722,1.197722,\
		1.197722,1.197722,1.216472,1.216472,1.216472,1.216472,0.984188,0.984188,0.984188,0.984188,\
		1.115163,1.115163,1.115163,1.115163,1.362953,1.362953,1.362953,1.362953,0.909270,0.909270,\
		0.909270,0.909270,1.175255,1.175255,1.175255,1.175255,1.194047,1.194047,1.194047,1.194047,\
		0.926949,0.926949,0.926949,0.926949,1.579459,1.579459,1.579459,1.579459,1.461178,1.461178,\
		1.461178,1.461178,0.919336,0.919336,0.919336,0.919336,1.040578,1.040578,1.040578,1.040578,\
		1.032111,1.032111,1.032111,1.032111,0.933315,0.933315,0.933315,0.933315,0.959561,0.959561,\
		0.959561,0.959561,1.090474,1.090474,1.090474,1.090474,1.036574,1.036574,1.036574,1.036574,\
		0.900257,0.900257,0.900257,0.900257,1.033866,1.033866,1.033866,1.033866,1.135500,1.135500,\
		1.135500,1.135500,0.977064,0.977064,0.977064,0.977064,0.934492,0.934492,0.934492,0.934492,\
		1.162880,1.162880,1.162880,1.162880,1.289806,1.289806,1.289806,1.289806,1.006958,1.006958,\
		1.006958,1.006958,0.928492,0.928492,0.928492,0.928492,1.191590,1.191590,1.191590,1.191590,\
		1.160928,1.160928,1.160928,1.160928,0.936137,0.936137,0.936137,0.936137,1.032687,1.032687,\
		1.032687,1.032687,1.238967,1.238967,1.238967,1.238967,0.962985,0.962985,0.962985,0.962985,\
		0.939073,0.939073,0.939073,0.939073,0.977581,0.977581,0.977581,0.977581,0.970448,0.970448,\
		0.970448,0.970448,1.005117,1.005117,1.005117,1.005117,0.952309,0.952309,0.952309,0.952309,\
		0.937934,0.937934,0.937934,0.937934,0.988830,0.988830,0.988830,0.988830,0.988068,0.988068,\
		0.988068,0.988068,0.959186,0.959186,0.959186,0.959186,0.984960,0.984960,0.984960,0.984960,\
		1.171568,1.171568,1.171568,1.171568,1.069086,1.069086,1.069086,1.069086,0.941382,0.941382,\
		0.941382,0.941382,1.134403,1.134403,1.134403,1.134403,0.929559,0.929559,0.929559,0.929559,\
		1.127018,1.127018,1.127018,1.127018,1.212695,1.212695,1.212695,1.212695,1.084545,1.084545,\
		1.084545,1.084545,0.964365,0.964365,0.964365,0.964365,1.190800,1.190800,1.190800,1.190800,\
		1.193340,1.193340,1.193340,1.193340,0.985248,0.985248,0.985248,0.985248,1.048564,1.048564,\
		1.048564,1.048564,1.116514,1.116514,1.116514,1.116514,0.908927,0.908927,0.908927,0.908927,\
		1.176945,1.176945,1.176945,1.176945,1.005654,1.005654,1.005654,1.005654,0.955642,0.955642,\
		0.955642,0.955642,1.150358,1.150358,1.150358,1.150358,1.141592,1.141592,1.141592,1.141592,\
		0.912113,0.912113,0.912113,0.912113,1.060346,1.060346,1.060346,1.060346,1.206823,1.206823,\
		1.206823,1.206823,1.118318,1.118318,1.118318,1.118318,0.914656,0.914656,0.914656,0.914656,\
		1.162172,1.162172,1.162172,1.162172,1.260670,1.260670,1.260670,1.260670,1.057328,1.057328,\
		1.057328,1.057328,1.018547,1.018547,1.018547,1.018547,1.317649,1.317649,1.317649,1.317649,\
		0.907617,0.907617,0.907617,0.907617,1.000485,1.000485,1.000485,1.000485,0.903974,0.903974,\
		0.903974,0.903974,1.080632,1.080632,1.080632,1.080632,1.108696,1.108696,1.108696,1.108696,\
		0.926796,0.926796,0.926796,0.926796,0.979487,0.979487,0.979487,0.979487]).astype("float")

	clr=np.array([[210,174,136,255],
		[210,174,136,255],
		[210,174,136,255],
		[210,174,136,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,164,128,255],
		[210,164,128,255],
		[210,164,128,255],
		[210,164,128,255],
		[210,163,128,254],
		[210,163,128,254],
		[210,163,128,254],
		[210,163,128,254],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,174,136,255],
		[210,174,136,255],
		[210,174,136,255],
		[210,174,136,255],
		[210,164,129,255],
		[210,164,129,255],
		[210,164,129,255],
		[210,164,129,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,166,130,254],
		[210,166,130,254],
		[210,166,130,254],
		[210,166,130,254],
		[210,168,132,255],
		[210,168,132,255],
		[210,168,132,255],
		[210,168,132,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,162,127,255],
		[210,162,127,255],
		[210,162,127,255],
		[210,162,127,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,166,130,255],
		[210,166,130,255],
		[210,166,130,255],
		[210,166,130,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,164,128,254],
		[210,164,128,254],
		[210,164,128,254],
		[210,164,128,254],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,165,130,255],
		[210,165,130,255],
		[210,165,130,255],
		[210,165,130,255],
		[210,168,132,255],
		[210,168,132,255],
		[210,168,132,255],
		[210,168,132,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,162,127,255],
		[210,162,127,255],
		[210,162,127,255],
		[210,162,127,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,165,129,255],
		[210,165,129,255],
		[210,165,129,255],
		[210,165,129,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,170,132,255],
		[210,170,132,255],
		[210,170,132,255],
		[210,170,132,255],
		[210,177,137,254],
		[210,177,137,254],
		[210,177,137,254],
		[210,177,137,254],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,171,133,254],
		[210,171,133,254],
		[210,171,133,254],
		[210,171,133,254],
		[210,178,139,254],
		[210,178,139,254],
		[210,178,139,254],
		[210,178,139,254],
		[210,172,134,254],
		[210,172,134,254],
		[210,172,134,254],
		[210,172,134,254],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,175,137,255],
		[210,175,137,255],
		[210,175,137,255],
		[210,175,137,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,166,130,254],
		[210,166,130,254],
		[210,166,130,254],
		[210,166,130,254],
		[210,157,124,255],
		[210,157,124,255],
		[210,157,124,255],
		[210,157,124,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,165,129,255],
		[210,165,129,255],
		[210,165,129,255],
		[210,165,129,255],
		[210,157,123,255],
		[210,157,123,255],
		[210,157,123,255],
		[210,157,123,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,166,130,255],
		[210,166,130,255],
		[210,166,130,255],
		[210,166,130,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,170,133,254],
		[210,170,133,254],
		[210,170,133,254],
		[210,170,133,254],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,178,139,254],
		[210,178,139,254],
		[210,178,139,254],
		[210,178,139,254],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,167,131,254],
		[210,167,131,254],
		[210,167,131,254],
		[210,167,131,254],
		[210,175,137,255],
		[210,175,137,255],
		[210,175,137,255],
		[210,175,137,255],
		[210,179,139,254],
		[210,179,139,254],
		[210,179,139,254],
		[210,179,139,254],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,156,123,255],
		[210,156,123,255],
		[210,156,123,255],
		[210,156,123,255],
		[210,180,140,255],
		[210,180,140,255],
		[210,180,140,255],
		[210,180,140,255],
		[210,166,130,255],
		[210,166,130,255],
		[210,166,130,255],
		[210,166,130,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,179,139,254],
		[210,179,139,254],
		[210,179,139,254],
		[210,179,139,254],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,175,137,255],
		[210,175,137,255],
		[210,175,137,255],
		[210,175,137,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,168,132,255],
		[210,168,132,255],
		[210,168,132,255],
		[210,168,132,255],
		[210,177,137,255],
		[210,177,137,255],
		[210,177,137,255],
		[210,177,137,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,163,128,255],
		[210,163,128,255],
		[210,163,128,255],
		[210,163,128,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,169,132,254],
		[210,169,132,254],
		[210,169,132,254],
		[210,169,132,254],
		[210,179,139,254],
		[210,179,139,254],
		[210,179,139,254],
		[210,179,139,254],
		[210,156,123,255],
		[210,156,123,255],
		[210,156,123,255],
		[210,156,123,255],
		[210,160,125,255],
		[210,160,125,255],
		[210,160,125,255],
		[210,160,125,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,175,136,254],
		[210,175,136,254],
		[210,175,136,254],
		[210,175,136,254],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,171,134,254],
		[210,171,134,254],
		[210,171,134,254],
		[210,171,134,254],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,166,130,254],
		[210,166,130,254],
		[210,166,130,254],
		[210,166,130,254],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,170,133,254],
		[210,170,133,254],
		[210,170,133,254],
		[210,170,133,254],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,175,136,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,168,131,255],
		[210,177,138,254],
		[210,177,138,254],
		[210,177,138,254],
		[210,177,138,254],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,137,255],
		[210,177,137,255],
		[210,177,137,255],
		[210,177,137,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,171,134,255],
		[210,171,134,255],
		[210,171,134,255],
		[210,171,134,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,178,139,255],
		[210,171,134,255],
		[210,171,134,255],
		[210,171,134,255],
		[210,171,134,255],
		[210,168,132,255],
		[210,168,132,255],
		[210,168,132,255],
		[210,168,132,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,174,136,255],
		[210,174,136,255],
		[210,174,136,255],
		[210,174,136,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,170,133,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,178,138,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,171,133,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,174,135,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,169,132,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,170,133,254],
		[210,170,133,254],
		[210,170,133,254],
		[210,170,133,254],
		[210,167,130,255],
		[210,167,130,255],
		[210,167,130,255],
		[210,167,130,255],
		[210,174,136,254],
		[210,174,136,254],
		[210,174,136,254],
		[210,174,136,254],
		[210,175,137,255],
		[210,175,137,255],
		[210,175,137,255],
		[210,175,137,255],
		[210,165,129,255],
		[210,165,129,255],
		[210,165,129,255],
		[210,165,129,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,176,137,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,173,135,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,172,134,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,179,139,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255],
		[210,177,138,255]]).astype("uint8")
	mdl=session.open_command.open_data("/Users/albertsmith/Documents/Dynamics/MF_MD_theory/Figures/backbone3D/scene_Y1.cxs")[0]
	session.models.add(mdl)
	rc(session,"~display")
	rc(session,"~ribbon")
	from chimerax.atomic.pbgroup import PseudobondGroup
	pbg_list=list()
	for mdl in session.models:
		if isinstance(mdl,PseudobondGroup):pbg_list.append(mdl)
	for mdl in pbg_list:mdl.delete()
	if len(session.models)>1:
		atoms=session.models[1].atoms
		rc(session,"display #1.1")
	else:
		atoms=session.models[0].atoms
		rc(session,"display #1")
	hide=getattr(atoms,"hides")
	hide[:]=1
	hide[di]=0
	setattr(atoms,"hides",hide)
	rc(session,"style ball")
	rc(session,"size stickRadius 0.2")
	rc(session,"color all tan")
	r0=getattr(atoms,"radii").copy()
	clr0=getattr(atoms,"colors").copy()
	r0[:]=.8
	r0[ids]=r
	clr0[ids]=clr
	setattr(atoms,"radii",r0)
	setattr(atoms,"colors",clr0)
	rc(session,"window 350 500")
	rc(session,"save /Users/albertsmith/Documents/Dynamics/MF_MD_theory/Figures/backbone3D/fr3_rho3.png ")
except:
	pass
finally:
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script415745.py")
	rc(session,"exit")
