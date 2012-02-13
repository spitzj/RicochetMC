// $Id: CosmicMuon.cc,v 1.3 2011/07/21 21:22:18 kelsey Exp $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        CosmicMuon.cc                                                //
//  Description: generator of cosmic ray muons according to the distributions //
//               in Cassiday (PRD 7, 2022 (1973) ) which uses the Groom       //
//               parameterization                                             //
//  Author:      unknown, adpated for use in CDMS framework by                //
//               D.H. Wright (SLAC)                                           //
//  Date:        24 May 2011                                                  //
//                                                                            //
//  20110716  M. Kelsey -- Drop unnecessary functions, improve const-ness     //
////////////////////////////////////////////////////////////////////////////////

#include "SurfaceCosmicMuon.hh"
#include "G4ThreeVector.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "Randomize.hh"
#include <math.h>


// To describe vertical flux of cosmic rays
// dN/dx = A * [exp(-B*h) + C*exp(-D*h)]
// where h is slant depth
//
// VerticalFlux_A = I_1 = 12.8e-6 (Cassiday) 
// VerticalFlux_C = I_2/I_1 = 1.13e-6/12.8e-6 (Cassiday)
// VerticalFlux_B = 1/lambda_1 = 1/366
// VerticalFlux_D = 1/lambda_2 = 1/794

SurfaceCosmicMuon::SurfaceCosmicMuon()
  : CDMSVParticleGenerator("CosmicMuon"),
    VerticalFlux_A(8.6e-6 /s/cm2), VerticalFlux_B(1./450./m),
    VerticalFlux_C(0.44/8.6), VerticalFlux_D(1./870./m),
// : VerticalFlux_A(12.8e-6 /s/cm2), VerticalFlux_B(1./366./m),
//   VerticalFlux_C(1.13/12.8), VerticalFlux_D(1./794./m),

   // Parameters for theta and energy distribution of cosmic rays
   ThetaEnergyDist_A(twopi*VerticalFlux_A*VerticalFlux_B),
   ThetaEnergyDist_B(VerticalFlux_B),
   ThetaEnergyDist_C(VerticalFlux_C*VerticalFlux_D/VerticalFlux_B),
   ThetaEnergyDist_D(VerticalFlux_D),

    // Parameters for theta distribution of cosmic rays
   ThetaDist_A(twopi*VerticalFlux_A),
   ThetaDist_C(VerticalFlux_C),

   // Parameters for energy distribution of cosmic rays
   EnergyDist_B(VerticalFlux_B),
   EnergyDist_D(VerticalFlux_D),

   // Loss rate is modeled piecewise linearly:  dE/dx = a(E + c)
   // For 1 < E < 10	(and E < 1)
   EnergyLossRate_1a(0.004255/m), EnergyLossRate_1c(0.1820/0.004255*GeV),

   // For 10 < E < 100
   EnergyLossRate_10a(0.0005956/m), EnergyLossRate_10c(0.2160/0.0005956*GeV),

   // For 100 < E < 1000
   EnergyLossRate_100a(0.0004263/m), EnergyLossRate_100c(0.2280/0.0004263*GeV),

   // For 100 < E < 10000+
   EnergyLossRate_1000a(0.0004425/m), EnergyLossRate_1000c(0.2020/0.0004425*GeV),

   // Additive constants for piecewise defined stopping distance
   StoppingDistance_0GeV(0.*m),
   StoppingDistance_D1(StoppingDistance_0GeV - log(1 + 0*GeV/EnergyLossRate_1c)
			/EnergyLossRate_1a),
   StoppingDistance_10GeV(StoppingDistance_D1 + log(1 + 10*GeV/EnergyLossRate_1c)
			   /EnergyLossRate_1a),
   StoppingDistance_D10(StoppingDistance_10GeV - log(1 + 10*GeV /EnergyLossRate_10c)
			 /EnergyLossRate_10a),
   StoppingDistance_100GeV(StoppingDistance_D10 + log(1 + 100*GeV/EnergyLossRate_10c)
			    /EnergyLossRate_10a),
   StoppingDistance_D100(StoppingDistance_100GeV
                  - log(1 + 100*GeV/EnergyLossRate_100c)/EnergyLossRate_100a),
   StoppingDistance_1000GeV(StoppingDistance_D100
                  + log(1 + 1000*GeV/EnergyLossRate_100c)/EnergyLossRate_100a),
   StoppingDistance_D1000(StoppingDistance_1000GeV
                  - log(1 + 1000*GeV/EnergyLossRate_1000c)/EnergyLossRate_1000a),

   muThetaMin(0.), muThetaMax(halfpi), muThetaBins(1000),
   muEnergyMin(1*GeV), muEnergyMax(3001*GeV), muEnergyBins(1000),
   MWEDepth(1000*m)
{
  Initialize();
}


void SurfaceCosmicMuon::Initialize()
{
  // Parameters for theta and energy distribution of cosmic rays
  ThetaEnergyDist_Bh = VerticalFlux_B * MWEDepth;
  ThetaEnergyDist_Dh = VerticalFlux_D * MWEDepth;

  // Parameters for theta distribution of cosmic rays
  ThetaDist_B = VerticalFlux_B * MWEDepth;
  ThetaDist_D = VerticalFlux_D * MWEDepth;

  // Parameters for energy distribution of cosmic rays
  G4double EiBh = ExponentialIntegral(-ThetaEnergyDist_Bh);
  EnergyDist_A = -ThetaEnergyDist_A*EiBh;
  EnergyDist_C = ThetaEnergyDist_C
	         *ExponentialIntegral(-ThetaEnergyDist_Dh)/EiBh;

  // Total muon rate per unit horizontal area (flux)
  TotalMuonFlux = Flux(muEnergyMin, muEnergyMax);
}


// Generate random energy, direction, and muon type

G4double SurfaceCosmicMuon::shootKineticEnergy() const {
  return RandomEnergy();
}

G4ThreeVector SurfaceCosmicMuon::shootDirection() const {
  static G4ThreeVector dir;	// Buffer to avoid some memory churn
  dir.setRThetaPhi(1., RandomTheta(), RandomPhi());
  return dir;
}

G4ParticleDefinition* SurfaceCosmicMuon::shootParticle() const {
  static G4ParticleDefinition* muplus  = G4MuonPlus::Definition();
  static G4ParticleDefinition* muminus = G4MuonMinus::Definition();

  return (G4UniformRand() < 0.5) ? muplus : muminus;
}


G4double SurfaceCosmicMuon::RandomTime(G4double area) const
{
  // Returns a random time based upon the distribution given by
  // TimeProb() (for a horizontal surface of the given area at a
  // depth given by MWEDepth).
  
  // Check if positive area
  if (area < 0) {
    G4cerr << "CosmicMuon::TimeProb  ";
    G4cerr << "Area argument must be positive";
    G4cerr << " (" << (area /m2) << " m2)" << G4endl;
    return -1.;
  }

  return 0. - log(1 - G4UniformRand()) / (TotalMuonFlux * area);
}


G4double SurfaceCosmicMuon::RandomPhi() const
{
  // Returns a random phi from -pi to pi.
  return (G4UniformRand()-0.5) * twopi;
}


G4double SurfaceCosmicMuon::RandomTheta() const
{
  // Sample cos(theta) according to eq. 3 of Mei & Hime,
  // Phys. Rev. D 73, 053004 (2006)
  // Code re-written by D.H. Wright (SLAC) 25 June 2011

  // Function increases monotonically with cos(theta) => max at theta = 0 
  G4double fmax = exp(-ThetaDist_B) + ThetaDist_C*exp(-ThetaDist_D);
  G4double cost;
  G4double f;

  do {
    // Sample cos(theta) evenly from 0 < theta < pi/2
    cost = G4UniformRand();
    f = 0.0;
    if (cost > DBL_MIN)
      f = (exp(-ThetaDist_B/cost) + ThetaDist_C*exp(-ThetaDist_D/cost) )/cost;
  } while (G4UniformRand() > f/fmax);

  return acos(cost);
}


G4double SurfaceCosmicMuon::RandomEnergy() const 
{
  // Generate normalized local energy spectrum using eq. 8 and from 
  // Mei and Hime, Phys. Rev. D 73, 053004 (2006)
  // Code re-written by D.H. Wright (SLAC) 25 June 2011

  // Groom parameters
  G4double b = 0.0004/m;
  G4double gamma = 3.77;
  G4double eps = 693*GeV;
  G4double term = eps*(1. - exp(-b*MWEDepth) );

  G4double rand = G4UniformRand();
  if (rand < DBL_MIN) rand = DBL_MIN;

  return (muEnergyMin + term)*pow(rand, 1./(1. - gamma) ) - term;
}


G4double
SurfaceCosmicMuon::ThetaEnergyDist(const G4double theta, const G4double energy) const
{
  // dN/dTheta/dE = A * tan(t) / g(E) * [exp(-Bx - Bh/cos(t))
  //                                     + C*exp(-Dx - Dh/cos(t))]
  //             where t is theta
  //             g(E) is energy loss rate
  //		     and x = x(E) is muon stopping distance
  // This method returns the flux at an angle theta and energy E
  // (integrated over phi) based upon PhysRevD Vol 7 Num 7 pp2022-31
  // (Cassiday et al.)
  // Valid for depths over 1000 mwe
  // Is properly normalized in Geant4 units!!!

  G4double tan_val;
  G4double stop_dist;
  G4double exp_b;
  G4double exp_d;
  G4double dist_val;
  
  // Check if theta is within allowable range
  if ((theta < 0) | (theta > pi)) {
    G4cerr << "CosmicMuon::ThetaEnergyDist  ";
    G4cerr << "Theta argument must be between 0 and pi";
    G4cerr << " (" << theta << ")" << G4endl;
    return -1.;
  }
	
  // Check if positive energy
  if (energy < 0) {
    G4cerr << "CosmicMuon::ThetaEnergyDist  ";
    G4cerr << "Energy argument must be positive";
    G4cerr << " (" << (energy /GeV) << " GeV)" << G4endl;
    return -1.;
  }
	
  // There are no "upward" moving cosmic muons
  // Special case: halfpi would cause divide by zero
  if (theta >= halfpi) return 0.;
  
  tan_val = tan(theta);
  stop_dist = StoppingDistance(energy);
  exp_b = exp(-ThetaEnergyDist_B*stop_dist
              - ThetaEnergyDist_Bh/cos(theta));
  exp_d = exp(-ThetaEnergyDist_D*stop_dist
              - ThetaEnergyDist_Dh/cos(theta));
  
  dist_val = ThetaEnergyDist_A*tan_val/EnergyLossRate(energy)
             * (exp_b + ThetaEnergyDist_C * exp_d);
  
  return dist_val;
}


G4double 
SurfaceCosmicMuon::Flux(const G4double min_energy, const G4double max_energy) const
{
  // This method returns the rate of cosmic muons per unit time
  // PER UNIT AREA between the given energies.
  // The area is taken to be horizontal at the depth given by
  // MWEDepth.
  // This is done by taking the integral of dN/dt/dE (where t is
  // theta) piecewise for the regions defined for g(E).
  // The indefinite integral is given by RateIndefiniteIntegral().

  if (min_energy < 0) {
    G4cerr << "CosmicMuon::Rate  ";
    G4cerr << "Arguments must be positive";
    G4cerr << " (" << (min_energy /GeV) << " GeV)" << G4endl;
    return -1.;
  }
	
  // Check if minimum energy is less than maximum energy
  if (min_energy > max_energy) {
    G4cerr << "CosmicMuon::Rate  ";
    G4cerr << "Minimum energy must be less than maximum energy";
    G4cerr << " (" << (min_energy /GeV) << " GeV > ";
    G4cerr << (max_energy /GeV) << " GeV)" << G4endl;
    return -1.;
  }
	
  G4double region_min;
  G4double region_max;
  G4double total_flux = 0;
	
  if (min_energy < 10 *GeV) {
    region_min = min_energy;
    if (max_energy < 10 *GeV)
      region_max = max_energy;
    else
      region_max = 10 *GeV;
		
    total_flux +=
       FluxIndefiniteIntegral(region_max, EnergyLossRate_1a,
                              EnergyLossRate_1c, StoppingDistance_D1)
     - FluxIndefiniteIntegral(region_min, EnergyLossRate_1a,
                              EnergyLossRate_1c, StoppingDistance_D1);
  }
	
  if ((min_energy < 100 *GeV) & (max_energy > 10 *GeV)) {
    if (min_energy > 10 *GeV)
      region_min = min_energy;
    else
      region_min = 10 *GeV;
    if (max_energy < 100 *GeV)
      region_max = max_energy;
    else
      region_max = 100 *GeV;
		
    total_flux +=
       FluxIndefiniteIntegral(region_max, EnergyLossRate_10a,
                              EnergyLossRate_10c, StoppingDistance_D10)
     - FluxIndefiniteIntegral(region_min, EnergyLossRate_10a,
                              EnergyLossRate_10c, StoppingDistance_D10);
  }
	
  if ((min_energy < 1000 *GeV) & (max_energy > 100 *GeV)) {
    if (min_energy > 100 *GeV)
      region_min = min_energy;
    else
      region_min = 100 *GeV;
    if (max_energy < 1000 *GeV)
      region_max = max_energy;
    else
      region_max = 1000 *GeV;
		
    total_flux += 
       FluxIndefiniteIntegral(region_max, EnergyLossRate_100a,
                              EnergyLossRate_100c, StoppingDistance_D100)
     - FluxIndefiniteIntegral(region_min, EnergyLossRate_100a,
                              EnergyLossRate_100c, StoppingDistance_D100);
  }
	
  if (max_energy > 1000 *GeV) {
    if (min_energy > 1000 *GeV)
      region_min = min_energy;
    else
      region_min = 1000 *GeV;
    region_max = max_energy;
		
    total_flux +=
       FluxIndefiniteIntegral(region_max, EnergyLossRate_1000a,
                              EnergyLossRate_1000c, StoppingDistance_D1000)
     - FluxIndefiniteIntegral(region_min, EnergyLossRate_1000a,
                              EnergyLossRate_1000c, StoppingDistance_D1000);
  }
	
  return total_flux;
}


G4double SurfaceCosmicMuon::StoppingDistance(const G4double energy) const
{
  // g(E) = (1/a)*ln(1 + cE) + C	defined piecewise
  // a and c are same as Muon Energy Loss Rate constants
  // C is determined by continuity at piecewise boundaries
  // This method returns the muon stopping distance at a given energy
  // based upon CERN 85-03 (Lohmann et al.) 
  // Is properly normalized in Geant4 units!!!
  // Returns using mwe

  // Check if positive energy
  if (energy < 0) {
    G4cerr << "CosmicMuon::StoppingDistance  ";
    G4cerr << "Argument must be positive";
    G4cerr << " (" << (energy /GeV) << " GeV)" << G4endl;
    return -1.;
  }
	
  // Energies 0 - 10 GeV
  if (energy < 10 *GeV)
    return log(1 + energy/EnergyLossRate_1c) / EnergyLossRate_1a
		       + StoppingDistance_D1;
	
  // Energies 10 - 100 GeV
  if (energy < 100 *GeV)
    return log(1 + energy/EnergyLossRate_10c) / EnergyLossRate_10a
		       + StoppingDistance_D10;
	
  // Energies 100 - 1000 GeV
  if (energy < 1000 *GeV)
    return log(1 + energy/EnergyLossRate_100c) / EnergyLossRate_100a
		       + StoppingDistance_D100;
	
  // Energies 1000 - 10000+ GeV
  return log(1 + energy/EnergyLossRate_1000c) / EnergyLossRate_1000a
	       + StoppingDistance_D1000;
}


G4double SurfaceCosmicMuon::EnergyLossRate(const G4double energy) const
{
  // g(E) = a(E + c)  (aka dE/dx) defined piecewise
  // This method returns the muon energy loss rate at a given energy
  // based upon CERN 85-03 (Lohmann et al.) 
  // Is properly normalized in Geant4 units!!!
  // Returns using mwe

  // Check if positive energy
  if (energy < 0) {
    G4cerr << "CosmicMuon::EnergyLossRate  ";
    G4cerr << "Argument must be positive";
    G4cerr << " (" << (energy /GeV) << " GeV)" << G4endl;
    return -1.;
  }
	
  // Energies 0 - 10 GeV
  if (energy < 10*GeV)
    return EnergyLossRate_1a*(energy + EnergyLossRate_1c);
	
  // Energies 10 - 100 GeV
  if (energy < 100*GeV)
    return EnergyLossRate_10a*(energy + EnergyLossRate_10c);
	
  // Energies 100 - 1000 GeV
  if (energy < 1000*GeV)
    return EnergyLossRate_100a*(energy + EnergyLossRate_100c);
	
  // Energies 1000 - 10000+ GeV
  return EnergyLossRate_1000a*(energy + EnergyLossRate_1000c);
}


G4double SurfaceCosmicMuon::ExponentialIntegral(const G4double x) const
{
  // E1(-x) = C + ln(x) + Sum[n=1 -> infinity, (-x)**n /n/n!]
  // where C is Euler's constant.
  // This method computes the -x branch of the exponential integral. 
  // Code re-written (DHW, June 2011) because original implementation did 
  // not include enough terms to achieve convergence, resulting in oscillating 
  // positive and  negative function values.  For argument values < -7, 
  // switch to the asymptotic series:
  // E1(-x) = [Sum(n=0 -> inf., (1/x)**n n!]/x/exp(-x)

  G4double ei_val = 1.;
  G4double term;
  G4int N;
  G4double i;
  
  if (x > 0) {
    G4cerr << "CosmicMuon::ExponentialIntegral  ";
    G4cerr << "Can not take positive argument";
    G4cerr << " (" << x << ")" << G4endl;

  } else if (x < -7) {
    // Asymptotic series 
    N = 7;
    term = 1.;
    for (G4int n = 1; n < N; n++) {
      i = G4double(n);
      term = i*term/x;
      ei_val += term;
    }
    ei_val /= x*exp(-x);
 
  } else {
    // Small argument series
    G4double EulersConstant = 0.577215664901532860606;
    ei_val = EulersConstant + log(-x);	
    term = x;
    ei_val += term;
	
    N = G4int(4.*(-x) + 2.);
    for (G4int n = 2; n < N; n++) {
      i = G4double(n);
      term = x*term*(i - 1.)/i/i;
      ei_val += term;
    }
  }
  return ei_val;
}


G4double SurfaceCosmicMuon::FluxIndefiniteIntegral(const G4double energy,
					 const G4double a, const G4double c,
					 const G4double d) const
{
  // Returns the indefinite integral of dN/dt/dE over dt and dE
  // (given by ThetaEnergyDist()), where t is theta, for an energy
  // loss rate g and stopping distance x of the form:
  // g(E) = a(E+c)
  // x(E) = ln(1+E/c)/a + d
  // The definite integral is the difference between the indefinite
  // integral evaluated at the two endpoints.
  // NOTE: since a,c,d are defined piecewise, any overall integral
  // must be done piecewise.
  
  G4double term1;
  G4double term2;
  
  term1 = 1./EnergyDist_B
          * pow(1 + energy/c, -EnergyDist_B/a)
          * exp(-EnergyDist_B * d);
  term2 = EnergyDist_C/EnergyDist_D
          * pow(1 + energy/c, -EnergyDist_D/a)
          * exp(-EnergyDist_D * d);
  
  return 0. - EnergyDist_A * (term1 + term2);
}


G4double SurfaceCosmicMuon::FunctionFMaximum(const G4double a) const
{
  // Returns the maximum of the function: 
  // F(x) = tan(x) * exp(-a/cos(x))
  //   occurs when cos(x) = (sqrt(1+4*a^2) - 1)/(2a)
  // Useful for computing an upper bound on dN/dTheta/dE for a
  // given E.

  G4double g_a;
  
  if (a < 0) {
    G4cerr << "CosmicMuon::FunctionFMaximum  ";
    G4cerr << "Can not take negative argument";
    G4cerr << " (" << a << ")" << G4endl;
    return -1.;
  }

  g_a = sqrt(1 + 4*a*a) + 1;	
  return sqrt(g_a) * exp(0 - g_a /2) /sqrt(2) /a;
}


// Validate internal functions by printing out selected values

void SurfaceCosmicMuon::Test() const
{
  G4cout << G4endl;	
  G4cout << "=========================" << G4endl;
  G4cout << " Distribution Constants: " << G4endl;
  G4cout << "=========================" << G4endl;
  G4cout << " Ei(-1) = " << ExponentialIntegral(-1)
         << ", Ei(-2) = " << ExponentialIntegral(-2) 
         << ", Ei(-3) = " << ExponentialIntegral(-3) << G4endl;
  G4cout << " Ei(-4) = " << ExponentialIntegral(-4)
         << ", Ei(-5) = " << ExponentialIntegral(-5)
         << ", Ei(-6) = " << ExponentialIntegral(-6) << G4endl;
  G4cout << ", Ei(-6.5) = " << ExponentialIntegral(-6.5)
         << ", Ei(-7) = " << ExponentialIntegral(-7) 
         << ", Ei(-7.5) = " << ExponentialIntegral(-7.5) << G4endl;
  G4cout << " Ei(-10) = " << ExponentialIntegral(-10)
         << ", Ei(-15) = " << ExponentialIntegral(-15)
         << ", Ei(-20) = " << ExponentialIntegral(-20) << G4endl;
  G4cout << " MWEDepth = " << MWEDepth/m << " m " << G4endl;

  G4cout << " ThetaDist_A = " << ThetaDist_A*s*cm2 << G4endl;
  G4cout << " ThetaDist_B = " << ThetaDist_B << G4endl;
  G4cout << " ThetaDist_C = " << ThetaDist_C << G4endl;
  G4cout << " ThetaDist_D = " << ThetaDist_D << G4endl;
 
  G4cout << " EnergyDist_A = " << EnergyDist_A*s*cm2*m << G4endl;
  G4cout << " EnergyDist_B = " << EnergyDist_B*m << G4endl;
  G4cout << " EnergyDist_C = " << EnergyDist_C << G4endl;
  G4cout << " EnergyDist_D = " << EnergyDist_D*m << G4endl;

  G4cout << " StoppingDistance_D1 = " << StoppingDistance_D1/m << G4endl;
  G4cout << " StoppingDistance_D10 = " << StoppingDistance_D10/m << G4endl;
  G4cout << " StoppingDistance_D100 = " << StoppingDistance_D100/m << G4endl;
  G4cout << " StoppingDistance_D1000 = " << StoppingDistance_D1000/m << G4endl;

  G4cout << " StoppingDistance_0GeV = " << StoppingDistance_0GeV/m << G4endl;
  G4cout << " StoppingDistance_10GeV = " << StoppingDistance_10GeV/m << G4endl;
  G4cout << " StoppingDistance_100GeV = " << StoppingDistance_100GeV/m << G4endl;
  G4cout << " StoppingDistance_1000GeV = " << StoppingDistance_1000GeV/m << G4endl;

  G4cout << " EnergyLossRate_1a = " << EnergyLossRate_1a*m << G4endl;
  G4cout << " EnergyLossRate_1c = " << EnergyLossRate_1c/GeV << G4endl;
  G4cout << " EnergyLossRate_10a = " << EnergyLossRate_10a*m << G4endl;
  G4cout << " EnergyLossRate_10c = " << EnergyLossRate_10c/GeV << G4endl;
  G4cout << " EnergyLossRate_100a = " << EnergyLossRate_100a*m << G4endl;
  G4cout << " EnergyLossRate_100c = " << EnergyLossRate_100c/GeV << G4endl;
  G4cout << " EnergyLossRate_1000a = " << EnergyLossRate_1000a*m << G4endl;
  G4cout << " EnergyLossRate_1000c = " << EnergyLossRate_1000c/GeV << G4endl;

  G4double energy;
  G4double muTestEnergy[30] =
     {  1.0,    2.0,    5.0,   10.0,   20.0,   40.0,   60.0,   80.0,  100.0,
      150.0,  200.0,  250.0,  300.0,  350.0,  400.0,  450.0,  500.0,  600.0,
      700.0,  800.0,  900.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0,
     4000.0, 4500.0, 5000.0};
  G4double cost = cos(5*deg);
  G4double tant = tan(5*deg);
  G4double num;
  G4double den;

  for (G4int i = 0; i < 30; i++) {
    energy = muTestEnergy[i]*GeV;
    num = ThetaEnergyDist(5*deg, energy)/tant;
    den = ThetaDist_A*( exp(-ThetaDist_B/cost) + ThetaDist_C*exp(-ThetaDist_D/cost) );
    G4cout << " energy = " << energy/GeV << " GeV, localSpec = " << (num/den)*GeV << G4endl;
  }

  // G4cout << " Flux-weighted mean energy = " << EWSum/FluxSum/GeV << G4endl;
 
  G4cout << " TotalMuonFlux = " << TotalMuonFlux*s*m2 << " /s/m2 " << G4endl;

  G4cout << "=========================" << G4endl;
  G4cout << G4endl;

  G4cout << "================" << G4endl;
  G4cout << " Random Thetas: " << G4endl;
  G4cout << "================" << G4endl;

  G4cout << "  First:   " << RandomTheta()/deg << " deg " << G4endl;
  G4cout << "  Second:  " << RandomTheta()/deg << " deg " << G4endl;
  G4cout << "  Third:   " << RandomTheta()/deg << " deg " << G4endl;
  G4cout << "  Fourth:  " << RandomTheta()/deg << " deg " << G4endl;
  G4cout << "  Fifth:   " << RandomTheta()/deg << " deg " << G4endl;

  G4cout << "================" << G4endl;
  G4cout << G4endl;

  G4cout << "==================" << G4endl;
  G4cout << " Random Energies: " << G4endl;
  G4cout << "==================" << G4endl;

  G4cout << "  First:   " << RandomEnergy()/GeV << " GeV" << G4endl;
  G4cout << "  Second:  " << RandomEnergy()/GeV << " GeV" << G4endl;
  G4cout << "  Third:   " << RandomEnergy()/GeV << " GeV" << G4endl;
  G4cout << "  Fourth:  " << RandomEnergy()/GeV << " GeV" << G4endl;
  G4cout << "  Fifth:   " << RandomEnergy()/GeV << " GeV" << G4endl;

  G4cout << "==================" << G4endl;
  G4cout << G4endl;

  G4cout << "====================="  << G4endl;
  G4cout << " Theta Distribution: "  << G4endl;
  G4cout << "====================="  << G4endl;

  G4cout << "  Minimum Theta:  " << muThetaMin/pi << " pi" << G4endl;
  G4cout << "  Maximum Theta:  " << muThetaMax/pi << " pi" << G4endl;
  G4cout << "  Bins:           " << muThetaBins << "" << G4endl;
  G4cout << "  Values:         " << G4endl;

  G4cout << "====================="  << G4endl;
  G4cout << G4endl;

  G4cout << "======================" << G4endl;
  G4cout << " Energy Distribution: " << G4endl;
  G4cout << "======================" << G4endl;

  G4cout << " Minimum Energy: " << muEnergyMin/GeV << " GeV" << G4endl;
  G4cout << " Maximum Energy: " << muEnergyMax/GeV << " GeV" << G4endl;
  G4cout << " Bins:           " << muEnergyBins << "" << G4endl;
  G4cout << " Values:         " << G4endl;

  G4cout << "====================="  << G4endl;
  G4cout << G4endl;
}




// $Id: CosmicMuon.cc,v 1.3 2011/07/21 21:22:18 kelsey Exp $
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:        CosmicMuon.cc                                                //
//  Description: generator of cosmic ray muons according to the distributions //
//               in Cassiday (PRD 7, 2022 (1973) ) which uses the Groom       //
//               parameterization                                             //
//  Author:      unknown, adpated for use in CDMS framework by                //
//               D.H. Wright (SLAC)                                           //
//  Date:        24 May 2011                                                  //
//                                                                            //
//  20110716  M. Kelsey -- Drop unnecessary functions, improve const-ness     //
////////////////////////////////////////////////////////////////////////////////

/*#include "SurfaceCosmicMuon.hh"
#include "G4ThreeVector.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "Randomize.hh"
#include "CLHEP/Random/RandGeneral.h"
#include <math.h>


// To describe vertical flux of cosmic rays
// dN/dx = A * [exp(-B*h) + C*exp(-D*h)]
// where h is slant depth
//
// VerticalFlux_A = I_1 = 12.8e-6 (Cassiday) 
// VerticalFlux_C = I_2/I_1 = 1.13e-6/12.8e-6 (Cassiday)
// VerticalFlux_B = 1/lambda_1 = 1/366
// VerticalFlux_D = 1/lambda_2 = 1/794

SurfaceCosmicMuon::SurfaceCosmicMuon() : CDMSVParticleGenerator("SurfaceCosmicMuon")
{
  G4int nbins = 2000;
  binsize = 0.5;
  minbin = 0.5;
  maxbin = 1000.;
  G4double muonMass = 105.6584;
  G4double thisbin = minbin;  
  G4int binCount = 0;
  G4double test;
  EDistFunc = new G4double[2000];

  G4cout << "SurfaceCosmicMuon::SurfaceCosmicMuon()" << G4cout;
  while(thisbin < maxbin)
  {
    // My random generator computes the momentum of the muon
    // but we need the kinetic energy; do the conversion here
    EDistFunc[binCount] = sqrt(pow(muonMass,2) + pow(MomentumDistribution(thisbin),2)) - muonMass;
    test = MomentumDistribution(thisbin);
    thisbin = thisbin + binsize;
    binCount++;
  }

  randGenEnergy = new CLHEP::RandGeneral(EDistFunc,nbins);

}


// Generate random energy, direction, and muon type
G4double SurfaceCosmicMuon::shootKineticEnergy() const {
  G4cout << "SurfaceCosmicMuon::shootKineticEnergy()" << G4endl;
  return RandomEnergy();
}

G4ThreeVector SurfaceCosmicMuon::shootDirection() const {
  G4cout << "SurfaceCosmicMuon::shootDirection()" << G4endl;
  static G4ThreeVector dir;	// Buffer to avoid some memory churn
  dir.setRThetaPhi(1., RandomTheta(), RandomPhi());
  return dir;
}

G4ParticleDefinition* SurfaceCosmicMuon::shootParticle() const {
  G4cout << "SurfaceCosmicMuon::shootParticle()" << G4endl;
  static G4ParticleDefinition* muplus  = G4MuonPlus::Definition();
  static G4ParticleDefinition* muminus = G4MuonMinus::Definition();

  return (G4UniformRand() < 0.5) ? muplus : muminus;
}


// This might be useful some day, when the measured
// rate is implemented.
G4double SurfaceCosmicMuon::RandomTime(G4double area) const
{
  return 0.;
}


// Assume uniformly distributed polar angle
G4double SurfaceCosmicMuon::RandomPhi() const
{
  // Returns a random phi from -pi to pi.
  return (G4UniformRand()-0.5) * twopi;
}



// Assume a simple cos^2(theta) distribution;
// use rejection sampling
G4double SurfaceCosmicMuon::RandomTheta() const
{
  G4double thetaMin = 0;
  G4double thetaMax = twopi/4.;
  G4double Fmax = 1;
  G4double Fmin = 0;

  G4double theta = G4UniformRand() * (thetaMax - thetaMin);
  G4double F = G4UniformRand() * (Fmax - Fmin);
  G4int count = 0;
  while(F > ThetaDistribution(theta) && count < 1000)
  {
    theta = G4UniformRand() * (thetaMax - thetaMin);
    count++;
  }
  if(count == 1000)
    G4cout << "CosmicMuon::RandomEnergy: "
	   << "Ran rejection sampling for 1000 iterations "
	   << "but did not find suitable point!!" << G4endl;
  G4cout << "thetaDist: " << ThetaDistribution(theta) << G4endl;
  return theta;
}



// Returns a random energy from parameterization in hep-ph/0102042
// Note that this does not take into account the zenith-angle
// dependence of the muon energy spectrum.  This should be added
// at some point. (FIXME!)
G4double SurfaceCosmicMuon::RandomEnergy() const
{
  // Draw random number from energy distribution.
  // Note that this gives a number on [0,1] that we
  // must convert to an energy in MeV.
  G4double randNumber = randGenEnergy->fire();

  return ((randNumber * (maxbin - minbin)) + minbin ) * 1000. * MeV;
}


// This is equation (10) of hep-ph/0102042
G4double SurfaceCosmicMuon::MomentumDistribution(const G4double momentum) const
{
  // Define parameters
  G4double H1 = 0.135;
  G4double H2 = -2.529;
  G4double H3 = -5.76;
  G4double S2 = -2.10;

  G4double y = log10(momentum);

  G4double H = H1 * (pow(y,3.)/2. - 5.*pow(y,2.)/2. + 3.*y) +
               H2 * (-2.*pow(y,3.)/3. + 3.*pow(y,2.) - 10.*y/3. + 1.) +
               H3 * (pow(y,3.)/6. - pow(y,2.)/2 + y/3.) +
               S2 * (pow(y,3.)/3. - 2.*pow(y,2.) + 11.*y/3. - 2.);

  return pow(10.,H);
}

G4double SurfaceCosmicMuon::ThetaDistribution(const G4double theta) const
{
  return pow(cos(theta),2);
}
*/
