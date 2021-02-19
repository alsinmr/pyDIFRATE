try:
	from chimerax.core.commands import run as rc

	import os
	import numpy as np

	di=np.array([0,1,2,3,19,20,21,22,28,36,\
		37,38,39,44,55,56,57,58,66,75,\
		76,77,78,82,91,92,93,94,100,113,\
		114,115,116,120,127,128,129,130,135,146,\
		147,148,149,153,160,161,162,163,164,167,\
		168,169,170,176,189,190,191,192,196,203,\
		204,205,206,211,222,223,224,225,229,236,\
		237,238,239,244,255,256,257,258,264,270,\
		271,272,273,277,286,287,288,295,301,302,\
		303,304,315,316,317,318,321,326,327,328,\
		329,334,338,339,340,341,345,352,353,354,\
		355,360,371,372,373,374,380,386,387,388,\
		389,394,400,401,402,403,407,416,417,418,\
		419,425,438,439,440,441,443,448,449,450,\
		451,457,470,471,472,473,478,489,490,491,\
		492,498,506,507,508,509,514,518,519,520,\
		521,527,540,541,542,543,549,555,556,557,\
		558,559,562,563,564,570,581,582,583,595,\
		596,597,598,609,610,611,612,617,621,622,\
		623,624,630,638,639,640,641,647,655,656,\
		657,658,666,679,680,681,682,687,698,699,\
		700,701,706,717,718,719,720,728,737,738,\
		739,740,742,747,748,749,750,751,754,755,\
		756,757,763,776,777,778,779,785,793,794,\
		795,796,801,812,813,814,815,821,827,828,\
		829,830,835,839,840,841,842,843,846,847,\
		848,849,857,870,871,872,873,877,884,885,\
		886,887,892,903,904,905,906,909,914,915,\
		916,917,922,926,927,928,929,938,947,948,\
		949,950,955,961,962,963,964,969,980,981,\
		982,983,989,997,998,999,1000,1006,1019,1020,\
		1021,1022,1028,1034,1035,1036,1037,1040,1045,1046,\
		1047,1048,1052,1059,1060,1061,1062,1067,1078,1079,\
		1080,1081,1088,1095,1096,1097,1098,1103,1114,1115,\
		1116,1117,1121,1130,1131,1132,1133,1138,1149,1150,\
		1151,1152,1160,1173,1174,1175,1176,1181,1192,1193,\
		1194,1195,1203,1216,1217,1218,1219,1220,1223,1224,\
		1225,1228]).astype("uint32")
	rc(session,"open /Users/albertsmith/Documents/GitHub/pyDIFRATE/Struct/1d3z.pdb_5.pdb")
	rc(session,"~display")
	rc(session,"~ribbon")
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
except:
	pass
finally:
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script408196.py")
