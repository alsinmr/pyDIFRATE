try:
	from chimerax.core.commands import run as rc

	import os
	import numpy as np

	di=np.array([1065,1066,1067,1068,1069,1070,1071,1072,1073,1074,\
		1075,1076,1077,1078,1079,1080,1081,1082,1083,1084,\
		1085,1086,1087,1088,1089,1090,1091,1092,1093,1094,\
		1095,1096,1097,1098,1099,1100,1101,1102,1103,1104,\
		1105,1106,1107,1108,1109,1110,1111,1112,1113,1114,\
		1115,1116,1117,1118,1119,1120,1121,1122,1123,1124,\
		1125,1126,1127,1128,1129,1130,1131,1132,1133,1134,\
		1135,1136,1137,1138,1139,1140,1141,1142,1143,1144,\
		1145,1146,1147,1148,1149,1150,1151,1152,1153,1154,\
		1155,1156,1157,1158,1159,1160,1161,1162,1163,1164,\
		1165,1166,1167,1168,1169,1170,1171,1172,1173,1174,\
		1175,1176,1177,1178,1179,1180,1181,1182,1183,1184,\
		1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,\
		1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,\
		1205,1206,1207,1208,1209,1210,1211,1212,1213,1214,\
		1215,1216,1217,1218,1219,1220,1221,1222,1223,1224,\
		1225,1226,1227,1228,1229,1230,1231,1232,1233,1234,\
		1235,1236,1237,1238,1239,1240,1241,1242,1243,1244,\
		1245,1246,1247,1248,1249,1250,1251,1252,1253,1254,\
		1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,\
		1265,1266,1267,1268,1269,1270,1271,1272,1273,1274,\
		1275,1276,1277,1278,1279,1280,1281,1282,1283,1284,\
		1285,1286,1287,1288,1289,1290,1291,1292,1293,1294,\
		1295,1296,1297,1298,1299,1300,1301,1302,1303,1304,\
		1305,1306,1307,1308,1309,1310,1311,1312,1313,1314,\
		1315,1316,1317,1318,1319,1320,1321,1322,1323,1324,\
		1325,1326,1327,1328,1329,1330,1331,1332,1333,1334,\
		1335,1336,1337,1338,1339,1340,1341,1342,1343,1344,\
		1345,1346,1347,1348,1349,1350,1351,1352,1353,1354,\
		1355,1356,1357,1358,1359,1360,1361,1362,1363,1364,\
		1365,1366,1367,1368,1369,1370,1371,1372,1373,1374,\
		1375,1376,1377,1378,1379,1380,1381,1382,1383,1384,\
		1385,1386,1387,1388,1389,1390,1391,1392,1393,1394,\
		1395,1396,1397,1398,1399,1400,1401,1402,1403,1404,\
		1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,\
		1415,1416,1417,1418,1419,1420,1421,1422,1423,1424,\
		1425,1426,1427,1428,1429,1430,1431,1432,1433,1434,\
		1435,1436,1437,1438,1439,1440,1441,1442,1443,1444,\
		1445,1446,1447,1448,1449,1450,1451,1452,1453,1454,\
		1455,1456,1457,1458,1459,1460,1461,1462,1463,1464,\
		1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,\
		1475,1476,1477,1478,1479,1480,1481,1482,1483,1484,\
		1485,1486,1487,1488,1489,1490,1491,1492,1493,1494,\
		1495,1496,1497,1498,1499,1500,1501,1502,1503,1504,\
		1505,1506,1507,1508,1509,1510,1511,1512,1513,1514,\
		1515,1516,1517,1518,1519,1520,1521,1522,1523,1524,\
		1525,1526,1527,1528,1529,1530,1531,1532,1533,1534,\
		1535,1536,1537,1538,1539,1540,1541,1542,1543,1544,\
		1545,1546,1547,1548,1549,1550,1551,1552,1553,1554,\
		1555,1556,1557,1558,1559,1560,1561,1562,1563,1564,\
		1565,1566,1567,1568,1569,1570,1571,1572,1573,1574,\
		1575,1576,1577,1578,1579,1580,1581,1582,1583,1584,\
		1585,1586,1587,1588,1589,1590,1591,1592,1593,1594,\
		1595,1596,1597,1598,1599,1600,1601,1602,1603,1604,\
		1605,1606,1607,1608,1609,1610,1611,1612,1613,1614,\
		1615,1616,1617,1618,1619,1620,1621,1622,1623,1624,\
		1625,1626,1627,1628,1629,1630,1631,1632,1633,1634,\
		1635,1636,1637,1638,1639,1640,1641,1642,1643,1644,\
		1645,1646,1647,1648,1649,1650,1651,1652,1653,1654,\
		1655,1656,1657,1658,1659,1660,1661,1662,1663,1664,\
		1665,1666,1667,1668,1669,1670,1671,1672,1673,1674,\
		1675,1676,1677,1678,1679,1680,1681,1682,1683,1684,\
		1685,1686,1687,1688,1689,1690,1691,1692,1693,1694,\
		1695,1696,1697,1698,1699,1700,1701,1702,1703,1704,\
		1705,1706,1707,1708,1709,1710,1711,1712,1713,1714,\
		1715,1716,1717,1718,1719,1720,1721,1722,1723,1724,\
		1725,1726,1727,1728,1729,1730,1731,1732,1733,1734,\
		1735,1736,1737,1738,1739,1740,1741,1742,1743,1744,\
		1745,1746,1747,1748,1749,1750,1751,1752,1753,1754,\
		1755,1756,1757,1758,1759,1760,1761,1762,1763,1764,\
		1765,1766,1767,1768,1769,1770,1771,1772,1773,1774,\
		1775,1776,1777,1778,1779,1780,1781,1782,1783,1784,\
		1785,1786,1787,1788,1789,1790,1791,1792,1793,1794,\
		1795,1796,1797,1798,1799,1800,1801,1802,1803,1804,\
		1805,1806,1807,1808,1809,1810,1811,1812,1813,1814,\
		1815,1816,1817,1818,1819,1820,1821,1822,1823,1824,\
		1825,1826,1827,1828,1829,1830,1831,1832,1833,1834,\
		1835,1836,1837,1838,1839,1840,1841,1842,1843,1844,\
		1845,1846,1847,1848,1849,1850,1851,1852,1853,1854,\
		1855,1856,1857,1858,1859,1860,1861,1862,1863,1864,\
		1865,1866,1867,1868,1869,1870,1871,1872,1873,1874,\
		1875,1876,1877,1878,1879,1880,1881,1882,1883,1884,\
		1885,1886,1887,1888,1889,1890,1891,1892,1893,1894,\
		1895,1896,1897,1898,1899,1900,1901,1902,1903,1904,\
		1905,1906,1907,1908,1909,1910,1911,1912,1913,1914,\
		1915,1916,1917,1918,1919,1920,1921,1922,1923,1924,\
		1925,1926,1927,1928,1929,1930,1931,1932,1933,1934,\
		1935,1936,1937,1938,1939,1940,1941,1942,1943,1944,\
		1945,1946,1947,1948,1949,1950,1951,1952,1953,1954,\
		1955,1956,1957,1958,1959,1960,1961,1962,1963,1964,\
		1965,1966,1967,1968,1969,1970,1971,1972,1973,1974,\
		1975,1976,1977,1978,1979,1980,1981,1982,1983,1984,\
		1985,1986,1987,1988,1989,1990,1991,1992,1993,1994,\
		1995,1996,1997,1998,1999,2000,2001,2002,2003,2004,\
		2005,2006,2007,2008,2009,2010,2011,2012,2013,2014,\
		2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,\
		2025,2026,2027,2028,2029,2030,2031,2032,2033,2034,\
		2035,2036,2037,2038,2039,2040,2041,2042,2043,2044,\
		2045,2046,2047,2048,2049,2050,2051,2052,2053,2054,\
		2055,2056,2057,2058,2059,2060,2061,2062,2063,2064,\
		2065,2066,2067,2068,2069,2070,2071,2072,2073,2074,\
		2075,2076,2077,2078,2079,2080,2081,2082,2083,2084,\
		2085,2086,2087,2088,2089,2090,2091,2092,2093,2094,\
		2095,2096,2097,2098,2099,2100,2101,2102,2103,2104,\
		2105,2106,2107,2108,2109,2110,2111,2112,2113,2114,\
		2115,2116,2117,2118,2119,2120,2121,2122,2123,2124,\
		2125,2126,2127,2128,2129]).astype("uint32")

	ids=np.array([1073,1074,1075,1076,1080,1081,1082,1083,1102,1103,\
		1104,1105,1114,1115,1116,1117,1121,1122,1123,1124,\
		1133,1134,1135,1136,1137,1138,1139,1140,1203,1204,\
		1205,1206,1249,1250,1251,1252,1256,1257,1258,1259,\
		1358,1359,1360,1361,1394,1395,1396,1397,1398,1399,\
		1400,1401,1430,1431,1432,1433,1434,1435,1436,1437,\
		1467,1468,1469,1470,1471,1472,1473,1474,1483,1484,\
		1485,1486,1487,1488,1489,1490,1511,1512,1513,1514,\
		1521,1522,1523,1524,1531,1532,1533,1534,1546,1547,\
		1548,1549,1550,1551,1552,1553,1593,1594,1595,1596,\
		1600,1601,1602,1603,1636,1637,1638,1639,1643,1644,\
		1645,1646,1748,1749,1750,1751,1752,1753,1754,1755,\
		1793,1794,1795,1796,1797,1798,1799,1800,1809,1810,\
		1811,1812,1813,1814,1815,1816,1911,1912,1913,1914,\
		1915,1916,1917,1918,1930,1931,1932,1933,1934,1935,\
		1936,1937,1946,1947,1948,1949,1953,1954,1955,1956]).astype("uint32")

	r=np.array([1.409468,1.411378,1.405033,1.411993,1.233628,1.235019,1.232883,1.232983,2.289579,2.291232,\
		2.288982,2.288524,1.767868,1.765805,1.771496,1.766303,1.127194,1.127279,1.126127,1.128175,\
		1.460224,1.457449,1.456615,1.466607,1.226928,1.234177,1.225127,1.221480,1.769664,1.770852,\
		1.768382,1.769757,1.093329,1.092764,1.092078,1.095146,1.321810,1.322129,1.318407,1.324894,\
		2.462949,2.460618,2.464393,2.463836,1.275611,1.274964,1.277532,1.274338,1.339677,1.335304,\
		1.344177,1.339550,1.621273,1.622010,1.623789,1.618021,2.111307,2.113791,2.112903,2.107226,\
		1.759136,1.744509,1.784843,1.748057,1.702190,1.702129,1.710986,1.693456,1.574062,1.570816,\
		1.582367,1.569002,1.565004,1.570134,1.571586,1.553291,2.002662,1.975911,2.007927,2.024147,\
		2.021644,2.014041,2.017589,2.033303,2.065122,2.058412,2.068511,2.068443,1.732327,1.721613,\
		1.734760,1.740609,1.696237,1.684123,1.722589,1.682000,1.508046,1.507066,1.506194,1.510877,\
		1.147772,1.147408,1.146150,1.149757,1.483796,1.481010,1.489943,1.480433,1.128270,1.131771,\
		1.125109,1.127931,1.635046,1.636755,1.608269,1.660112,2.163296,2.177458,2.145228,2.167202,\
		1.343590,1.340213,1.349512,1.341045,1.588352,1.577244,1.603634,1.584177,1.494594,1.490408,\
		1.504704,1.488669,1.560278,1.560219,1.567043,1.553572,1.197899,1.196701,1.197998,1.198998,\
		1.196039,1.198047,1.195094,1.194975,1.575416,1.568862,1.588512,1.568873,2.072821,2.075538,\
		2.075605,2.067319,1.224165,1.222684,1.224720,1.225089,1.301350,1.304350,1.305782,1.293917]).astype("float")

	clr=np.array([[215,173,123,255],
		[215,173,123,255],
		[215,173,124,255],
		[215,173,123,255],
		[213,175,129,255],
		[213,175,129,255],
		[213,175,129,255],
		[213,175,129,255],
		[225,161,96,255],
		[225,161,96,255],
		[225,161,96,255],
		[225,161,96,255],
		[219,168,112,255],
		[219,168,112,255],
		[219,168,112,255],
		[219,168,112,255],
		[212,176,132,255],
		[212,176,132,255],
		[212,177,132,255],
		[212,176,132,255],
		[216,172,122,255],
		[216,172,122,255],
		[216,172,122,255],
		[216,172,122,255],
		[213,175,129,254],
		[213,175,129,255],
		[213,175,129,255],
		[213,175,129,255],
		[219,168,112,255],
		[219,168,112,255],
		[219,168,112,255],
		[219,168,112,255],
		[212,177,133,255],
		[212,177,133,255],
		[212,177,133,255],
		[212,177,133,255],
		[214,174,126,255],
		[214,174,126,254],
		[214,174,126,255],
		[214,174,126,255],
		[227,159,90,255],
		[227,159,90,255],
		[227,159,90,255],
		[227,159,90,255],
		[214,175,128,255],
		[214,175,128,255],
		[214,174,128,255],
		[214,175,128,254],
		[214,174,126,255],
		[214,174,126,255],
		[214,174,126,255],
		[214,174,126,255],
		[218,170,117,255],
		[218,170,117,255],
		[218,170,117,255],
		[218,170,117,255],
		[223,163,101,255],
		[223,163,101,255],
		[223,163,101,255],
		[223,164,101,255],
		[219,168,112,255],
		[219,168,113,255],
		[219,168,112,255],
		[219,168,113,255],
		[219,169,114,255],
		[219,169,114,255],
		[219,169,114,255],
		[218,169,115,255],
		[217,171,118,255],
		[217,171,118,255],
		[217,170,118,255],
		[217,171,118,255],
		[217,171,119,255],
		[217,171,118,254],
		[217,171,118,255],
		[217,171,119,255],
		[222,165,105,255],
		[222,165,106,255],
		[222,165,105,255],
		[222,165,104,255],
		[222,165,104,255],
		[222,165,104,255],
		[222,165,104,255],
		[222,164,104,255],
		[223,164,103,255],
		[223,164,103,255],
		[223,164,103,255],
		[223,164,103,255],
		[219,168,113,255],
		[219,169,114,255],
		[219,168,113,255],
		[219,168,113,255],
		[218,169,114,255],
		[218,169,115,255],
		[219,169,114,255],
		[218,169,115,255],
		[216,171,120,255],
		[216,171,120,255],
		[216,171,120,255],
		[216,171,120,255],
		[212,176,132,255],
		[212,176,132,255],
		[212,176,132,255],
		[212,176,132,255],
		[216,172,121,255],
		[216,172,121,255],
		[216,172,121,255],
		[216,172,121,255],
		[212,176,132,254],
		[212,176,132,254],
		[212,177,132,255],
		[212,176,132,255],
		[218,170,116,255],
		[218,170,116,255],
		[217,170,117,255],
		[218,169,116,255],
		[224,163,100,255],
		[224,163,99,254],
		[224,163,100,255],
		[224,163,100,255],
		[214,174,126,255],
		[214,174,126,255],
		[215,174,125,255],
		[214,174,126,254],
		[217,170,118,255],
		[217,171,118,254],
		[217,170,117,255],
		[217,170,118,254],
		[216,172,121,255],
		[216,172,121,255],
		[216,171,120,255],
		[216,172,121,255],
		[217,171,119,255],
		[217,171,119,255],
		[217,171,118,255],
		[217,171,119,255],
		[213,176,130,255],
		[213,176,130,255],
		[213,176,130,255],
		[213,176,130,255],
		[213,176,130,255],
		[213,176,130,255],
		[213,176,130,255],
		[213,176,130,254],
		[217,171,118,255],
		[217,171,118,255],
		[217,170,118,255],
		[217,171,118,255],
		[223,164,103,254],
		[223,164,102,255],
		[223,164,102,255],
		[223,164,103,255],
		[213,175,129,255],
		[213,175,129,255],
		[213,175,129,255],
		[213,175,129,255],
		[214,174,127,255],
		[214,174,127,255],
		[214,174,127,255],
		[214,174,127,255]]).astype("uint8")
	mdl=session.open_command.open_data("/Users/albertsmith/Documents/GitHub/pyDIFRATE/Struct/HETs_MET_4pw.xtc_102001.pdb")[0]
	session.models.add(mdl)
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
except:
	pass
finally:
	os.remove("/Users/albertsmith/Documents/GitHub/pyDIFRATE/chimera/chimera_script880959.py")
