try:
	from chimerax.core.commands import run as rc

	import os
	import numpy as np

	di=np.array([0,1,2,16,17,18,38,39,40,57,\
		58,59,69,70,71,79,80,81,98,99,\
		100,114,115,116,121,122,123,145,146,147,\
		159,160,161,170,171,172,173,174,177,178,\
		179,180,181,182,202,203,204,214,215,216,\
		217,218,219,220,221,230,231,232,233,234,\
		235,257,258,259,271,272,273,286,287,288,\
		301,302,303,325,326,327,328,329,332,333,\
		334,335,336,337,359,360,361,362,363,364,\
		365,369,370,371,375,376,377,392,393,394,\
		395,396,397,398,399,405,406,407,411,412,\
		413,418,419,420,432,433,434,435,436,437,\
		438,442,443,444,448,449,450,464,465,466,\
		478,479,480,481,482,485,486,487,488,489,\
		490,491,492,495,496,497,498,499,500,508,\
		509,510,527,528,529,544,545,546,551,552,\
		553,558,559,560,577,578,579,601,602,603,\
		620,621,622,631,632,633,643,644,645,660,\
		661,662,674,675,676,688,689,690,702,703,\
		704,713,714,715,716,717,718,719,723,724,\
		725,729,730,731,744,745,746,758,759,760,\
		774,775,776,790,791,792,797,798,799,819,\
		820,821,826,827,828,841,842,843,852,853,\
		854,876,877,878,879,880,881,882,886,887,\
		888,892,893,894,895,896,897,898,899,905,\
		906,907,911,912,913,914,915,916,917,918,\
		927,928,929,930,931,932,937,938,939,951,\
		952,953,966,967,968,987,988,989,994,995,\
		996,1001,1002,1003,1023,1024,1025,1030,1031,1032,\
		1050,1051,1052,1074,1075,1076,1086,1087,1088,1100,\
		1101,1102,1117,1118,1119,1134,1135,1136,1151,1152,\
		1153,1168,1169,1170,1185,1186,1187]).astype("uint32")