////////////////////////////////////////////////////////////////////////
// $Id: Ac228Lines.cc,v 1.1 2011/07/06 21:15:29 kevmc Exp $
//  File:        Ac228Lines.hh                                         //
//  Description: Generate special Ac228 gammas for CDMS testing       //
//                                                                    //
//  Author:      Kevin McCarthy (MIT)                                 //
//  Date:        20110706                                          //
//                                                                    //
//  20110706  K. McCarthy -- modeled from Am241Lines.cc               //
//////////////////////////////////////////////////////////////////////// 

#include "CDMSsources/Ac228Lines.hh"


// Constructor fills spectrum

Ac228Lines::Ac228Lines(const G4ThreeVector& dir, G4int verbose)
  : CDMSGammaLines("Ac228Lines", dir, verbose) {
    	AddLine(42.46*keV, 0.009);
	AddLine(56.96*keV, 0.019);
	AddLine(57.766*keV, 0.47);
	AddLine(77.34*keV, 0.026);
	AddLine(99.509*keV, 1.26);
	AddLine(100.41*keV, 0.093);
	AddLine(114.54*keV, 0.0098);
	AddLine(129.065*keV, 2.42);
	AddLine(135.51*keV, 0.018);
	AddLine(137.91*keV, 0.024);
	AddLine(141.01*keV, 0.050);
	AddLine(145.849*keV, 0.158);
	AddLine(153.977*keV, 0.722);
	AddLine(168.42*keV, 0.0030);
	AddLine(168.65*keV, 0.010);
	AddLine(173.964*keV, 0.035);
	AddLine(184.54*keV, 0.070);
	AddLine(191.353*keV, 0.123);
	AddLine(199.407*keV, 0.315);
	AddLine(204.026*keV, 0.112);
	AddLine(209.253*keV, 3.89);
	AddLine(214.85*keV, 0.029);
	AddLine(223.80*keV, 0.054);
	AddLine(231.42*keV, 0.025);
	AddLine(257.49*keV, 0.030);
	AddLine(263.62*keV, 0.040);
	AddLine(270.245*keV, 3.46);
	AddLine(278.70*keV, 0.031);
	AddLine(282.00*keV, 0.072);
	AddLine(321.646*keV, 0.226);
	AddLine(326.04*keV, 0.033);
	AddLine(327.45*keV, 0.12);
	AddLine(328.000*keV, 2.95);
	AddLine(332.370*keV, 0.40);
	AddLine(338.320*keV, 11.27);
	AddLine(340.98*keV, 0.369);
	AddLine(356.94*keV, 0.0170);
	AddLine(372.57*keV, 0.0067);
	AddLine(377.99*keV, 0.025);
	AddLine(384.63*keV, 0.0067);
	AddLine(389.12*keV, 0.0103);
	AddLine(397.94*keV, 0.027);
	AddLine(399.62*keV, 0.029);
	AddLine(409.462*keV, 1.92);
	AddLine(416.30*keV, 0.0132);
	AddLine(419.42*keV, 0.021);
	AddLine(440.44*keV, 0.121);
	AddLine(449.21*keV, 0.048);
	AddLine(452.51*keV, 0.015);
	AddLine(457.35*keV, 0.0150);
	AddLine(463.004*keV, 4.40);
	AddLine(466.40*keV, 0.029);
	AddLine(470.20*keV, 0.013);
	AddLine(471.76*keV, 0.033);
	AddLine(474.79*keV, 0.022);
	AddLine(478.40*keV, 0.209);
	AddLine(480.94*keV, 0.023);
	AddLine(490.33*keV, 0.0111);
	AddLine(492.30*keV, 0.0235);
	AddLine(497.49*keV, 0.0059);
	AddLine(503.823*keV, 0.182);
	AddLine(508.959*keV, 0.45);
	AddLine(515.06*keV, 0.049);
	AddLine(520.151*keV, 0.067);
	AddLine(523.131*keV, 0.103);
	AddLine(540.68*keV, 0.026);
	AddLine(546.45*keV, 0.201);
	AddLine(548.73*keV, 0.023);
	AddLine(555.07*keV, 0.046);
	AddLine(562.500*keV, 0.87);
	AddLine(570.88*keV, 0.182);
	AddLine(572.29*keV, 0.150);
	AddLine(583.41*keV, 0.111);
	AddLine(590.65*keV, 0.017);
	AddLine(610.64*keV, 0.023);
	AddLine(616.20*keV, 0.080);
	AddLine(620.33*keV, 0.080);
	AddLine(623.27*keV, 0.011);
	AddLine(627.23*keV, 0.014);
	AddLine(629.40*keV, 0.045);
	AddLine(634.18*keV, 0.0106);
	AddLine(640.34*keV, 0.054);
	AddLine(649.03*keV, 0.040);
	AddLine(651.48*keV, 0.090);
	AddLine(660.1*keV, 0.005);
	AddLine(663.88*keV, 0.028);
	AddLine(666.47*keV, 0.005);
	AddLine(672.00*keV, 0.026);
	AddLine(674.16*keV, 0.109);
	AddLine(677.07*keV, 0.062);
	AddLine(684.0*keV, 0.019);
	AddLine(688.11*keV, 0.067);
	AddLine(692.47*keV, 0.0056);
	AddLine(699.08*keV, 0.037);
	AddLine(701.747*keV, 0.173);
	AddLine(707.41*keV, 0.155);
	AddLine(718.31*keV, 0.019);
	AddLine(726.863*keV, 0.62);
	AddLine(737.72*keV, 0.037);
	AddLine(755.315*keV, 1.00);
	AddLine(770.2*keV, 0.0063);
	AddLine(772.291*keV, 1.49);
	AddLine(776.52*keV, 0.019);
	AddLine(778.1*keV, 0.022);
	AddLine(782.142*keV, 0.485);
	AddLine(791.44*keV, 0.010);
	AddLine(794.947*keV, 4.25);
	AddLine(813.77*keV, 0.0070);
	AddLine(816.62*keV, 0.030);
	AddLine(824.934*keV, 0.050);
	AddLine(830.486*keV, 0.540);
	AddLine(835.710*keV, 1.61);
	AddLine(840.377*keV, 0.91);
	AddLine(853.17*keV, 0.0088);
	AddLine(853.97*keV, 0.0031);
	AddLine(870.45*keV, 0.044);
	AddLine(873.11*keV, 0.031);
	AddLine(874.45*keV, 0.047);
	AddLine(877.39*keV, 0.014);
	AddLine(880.76*keV, 0.0062);
	AddLine(887.33*keV, 0.027);
	AddLine(901.26*keV, 0.016);
	AddLine(904.19*keV, 0.77);
	AddLine(911.204*keV, 25.8);
	AddLine(919.01*keV, 0.027);
	AddLine(921.98*keV, 0.0147);
	AddLine(922.5*keV, 0.0147);
	AddLine(924.3*keV, 0.0075);
	AddLine(930.93*keV, 0.0124);
	AddLine(939.87*keV, 0.009);
	AddLine(944.196*keV, 0.095);
	AddLine(947.982*keV, 0.106);
	AddLine(958.61*keV, 0.28);
	AddLine(964.766*keV, 4.99);
	AddLine(968.971*keV, 15.8);
	AddLine(975.98*keV, 0.050);
	AddLine(979.48*keV, 0.026);
	AddLine(987.88*keV, 0.077);
	AddLine(988.63*keV, 0.077);
	AddLine(1013.58*keV, 0.0046);
	AddLine(1016.44*keV, 0.019);
	AddLine(1017.92*keV, 0.0057);
	AddLine(1019.86*keV, 0.021);
	AddLine(1033.248*keV, 0.201);
	AddLine(1039.84*keV, 0.044);
	AddLine(1040.92*keV, 0.044);
	AddLine(1053.09*keV, 0.013);
	AddLine(1054.22*keV, 0.018);
	AddLine(1062.55*keV, 0.010);
	AddLine(1065.19*keV, 0.132);
	AddLine(1074.71*keV, 0.010);
	AddLine(1088.18*keV, 0.0059);
	AddLine(1095.679*keV, 0.129);
	AddLine(1103.41*keV, 0.0150);
	AddLine(1110.610*keV, 0.019);
	AddLine(1117.63*keV, 0.054);
	AddLine(1135.24*keV, 0.0098);
	AddLine(1142.85*keV, 0.0103);
	AddLine(1148.16*keV, 0.0059);
	AddLine(1153.52*keV, 0.139);
	AddLine(1157.14*keV, 0.0070);
	AddLine(1164.55*keV, 0.065);
	AddLine(1175.31*keV, 0.024);
	AddLine(1190.81*keV, 0.0062);
	AddLine(1217.03*keV, 0.021);
	AddLine(1229.40*keV, 0.0075);
	AddLine(1245.16*keV, 0.095);
	AddLine(1247.08*keV, 0.50);
	AddLine(1249.97*keV, 0.062);
	AddLine(1276.69*keV, 0.014);
	AddLine(1286.27*keV, 0.050);
	AddLine(1287.78*keV, 0.080);
	AddLine(1309.71*keV, 0.019);
	AddLine(1315.31*keV, 0.015);
	AddLine(1337.33*keV, 0.0049);
	AddLine(1344.59*keV, 0.0090);
	AddLine(1347.50*keV, 0.015);
	AddLine(1357.78*keV, 0.020);
	AddLine(1365.71*keV, 0.014);
	AddLine(1374.24*keV, 0.014);
	AddLine(1378.23*keV, 0.0059);
	AddLine(1385.39*keV, 0.0106);
	AddLine(1401.49*keV, 0.012);
	AddLine(1415.66*keV, 0.021);
	AddLine(1430.95*keV, 0.035);
	AddLine(1434.22*keV, 0.0080);
	AddLine(1438.01*keV, 0.0059);
	AddLine(1451.40*keV, 0.0106);
	AddLine(1459.138*keV, 0.83);
	AddLine(1469.71*keV, 0.020);
	AddLine(1480.37*keV, 0.016);
	AddLine(1495.93*keV, 0.86);
	AddLine(1501.57*keV, 0.46);
	AddLine(1529.05*keV, 0.057);
	AddLine(1537.87*keV, 0.047);
	AddLine(1548.65*keV, 0.038);
	AddLine(1557.10*keV, 0.178);
	AddLine(1559.78*keV, 0.020);
	AddLine(1571.52*keV, 0.0057);
	AddLine(1573.26*keV, 0.033);
	AddLine(1580.53*keV, 0.60);
	AddLine(1588.19*keV, 3.22);
	AddLine(1609.41*keV, 0.0077);
	AddLine(1625.06*keV, 0.255);
	AddLine(1630.627*keV, 1.51);
	AddLine(1638.281*keV, 0.47);
	AddLine(1666.523*keV, 0.178);
	AddLine(1671.64*keV, 0.0041);
	AddLine(1677.67*keV, 0.054);
	AddLine(1684.01*keV, 0.015);
	AddLine(1686.12*keV, 0.095);
	AddLine(1700.59*keV, 0.0101);
	AddLine(1702.44*keV, 0.048);
	AddLine(1706.17*keV, 0.0085);
	AddLine(1713.47*keV, 0.0054);
	AddLine(1721.4*keV, 0.0057);
	AddLine(1724.20*keV, 0.029);
	AddLine(1738.22*keV, 0.018);
	AddLine(1740.4*keV, 0.011);
	AddLine(1741.74*keV, 0.0080);
	AddLine(1745.28*keV, 0.0065);
	AddLine(1750.54*keV, 0.0080);
	AddLine(1758.11*keV, 0.035);
	AddLine(1772.2*keV, 0.0018);
	AddLine(1784.4*keV, 0.0059);
	AddLine(1787.3*keV, 0.0013);
	AddLine(1795.1*keV, 0.0021);
	AddLine(1797.5*keV, 0.0021);
	AddLine(1800.86*keV, 0.0044);
	AddLine(1823.21*keV, 0.044);
	AddLine(1826.7*keV, 0.0021);
	AddLine(1835.29*keV, 0.038);
	AddLine(1842.14*keV, 0.042);
	AddLine(1850.13*keV, 0.0044);
	AddLine(1870.81*keV, 0.0243);
	AddLine(1879.6*keV, 0.0013);
	AddLine(1887.12*keV, 0.090);
	AddLine(1900.14*keV, 0.0028);
	AddLine(1907.18*keV, 0.0119);
	AddLine(1915.9*keV, 0.0008);
	AddLine(1919.5*keV, 0.0021);
	AddLine(1929.78*keV, 0.0199);
	AddLine(1936.3*keV, 0.0021);
	AddLine(1944.20*keV, 0.0021);
	AddLine(1952.37*keV, 0.059);
	AddLine(1955.9*keV, 0.0008);
	AddLine(1958.4*keV, 0.0015);
	AddLine(1965.22*keV, 0.0204);
	AddLine(1971.9*keV, 0.0036);
	AddLine(1979.3*keV, 0.0018);
	AddLine(2000.9*keV, 0.0010);
	AddLine(2029.4*keV, 0.0018);

  NormalizeLines();
}

Ac228Lines::~Ac228Lines() {}
