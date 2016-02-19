/**
   C++ implementation of the TC builder

   S.Viret, G.Galbit : 09/12/15
**/

#include "../interface/TCBuilder.h"

TCBuilder::TCBuilder():TrackFitter(0){

}

TCBuilder::TCBuilder(int nb):TrackFitter(nb)
{
  m_nMissingHits = 1;                 //Maximum number of missing layers in a TC (from the number of layers in the pattern)
  m_bHardwareSimulation = false;      //Define if the Hardware binning and overflows are emulated
  updateThresholds();
}

TCBuilder::~TCBuilder()
{
}

void TCBuilder::initialize(){

}

void TCBuilder::setHardwareEmulation(bool hardwareEmulation)
{
  m_bHardwareSimulation = hardwareEmulation;
  updateThresholds();
}

void TCBuilder::updateThresholds(){

  if (m_bHardwareSimulation){
    //// Hardware format LUT constants ////
  
    //Barrel
    addThresholds( 0,  1,  2, SEC_BARREL, 0.002117, 0.298828);
    addThresholds( 0,  1,  8, SEC_BARREL, 0.005257, 2.695801);
    addThresholds( 0,  1,  9, SEC_BARREL, 0.010082, 2.861084);
    addThresholds( 0,  1, 10, SEC_BARREL, 0.012878, 3.115967);
    addThresholds( 0,  2,  8, SEC_BARREL, 0.002262, 2.535645);
    addThresholds( 0,  2,  9, SEC_BARREL, 0.005077, 2.578613);
    addThresholds( 0,  2, 10, SEC_BARREL, 0.007915, 2.765625);
    addThresholds( 1,  2,  8, SEC_BARREL, 0.001514, 2.562988);
    addThresholds( 1,  2,  9, SEC_BARREL, 0.004452, 2.700928);
    addThresholds( 1,  2, 10, SEC_BARREL, 0.007317, 2.892578);

    //Hybrid
    addThresholds( 0,  1,  2, SEC_HYBRID, 0.002125, 0.331055);
    addThresholds( 0,  1,  8, SEC_HYBRID, 0.004444, 2.867432);
    addThresholds( 0,  1,  9, SEC_HYBRID, 0.008244, 2.764160);
    addThresholds( 0,  1, 10, SEC_HYBRID, 0.009171, 3.031006);
    addThresholds( 0,  1, 11, SEC_HYBRID, 0.009350, 4.846436);
    addThresholds( 0,  1, 12, SEC_HYBRID, 0.009167, 5.155273);
    addThresholds( 0,  1, 13, SEC_HYBRID, 0.013248, 5.164795);
    addThresholds( 0,  2,  8, SEC_HYBRID, 0.002472, 2.536133);
    addThresholds( 0,  2,  9, SEC_HYBRID, 0.004913, 2.634033);
    addThresholds( 0,  2, 10, SEC_HYBRID, 0.004910, 2.887451);
    addThresholds( 0,  2, 11, SEC_HYBRID, 0.006683, 4.866699);
    addThresholds( 0,  2, 12, SEC_HYBRID, 0.006672, 5.009277);
    addThresholds( 0,  2, 13, SEC_HYBRID, 0.008404, 5.110840);
    addThresholds( 1,  2,  8, SEC_HYBRID, 0.002354, 2.580078);
    addThresholds( 1,  2,  9, SEC_HYBRID, 0.007805, 2.699951);
    addThresholds( 1,  2, 10, SEC_HYBRID, 0.003143, 2.849121);
    addThresholds( 1,  2, 11, SEC_HYBRID, 0.006020, 4.935059);
    addThresholds( 1,  2, 12, SEC_HYBRID, 0.006187, 5.085938);
    addThresholds( 1,  2, 13, SEC_HYBRID, 0.005692, 5.270752);

    //Endcap
    addThresholds( 0,  1,  3, SEC_ENDCAP, 0.003838, 0.457275);
    addThresholds( 0,  1,  4, SEC_ENDCAP, 0.007999, 0.611084);
    addThresholds( 0,  1,  5, SEC_ENDCAP, 0.002743, 0.502197);
    addThresholds( 0,  1, 11, SEC_ENDCAP, 0.006119, 5.451172);
    addThresholds( 0,  1, 12, SEC_ENDCAP, 0.007030, 6.014893);
    addThresholds( 0,  1, 13, SEC_ENDCAP, 0.012161, 7.055176);
    addThresholds( 0,  1, 14, SEC_ENDCAP, 0.017605, 7.991699);
    addThresholds( 0,  1, 15, SEC_ENDCAP, 0.027527, 7.594482);
    addThresholds( 0,  3,  4, SEC_ENDCAP, 0.002586, 0.788574);
    addThresholds( 0,  3,  5, SEC_ENDCAP, 0.003773, 1.233887);
    addThresholds( 0,  3,  6, SEC_ENDCAP, 0.004192, 1.924805);
    addThresholds( 0,  3,  7, SEC_ENDCAP, 0.006290, 2.382568);
    addThresholds( 0,  3, 11, SEC_ENDCAP, 0.004772, 5.415283);
    addThresholds( 0,  3, 12, SEC_ENDCAP, 0.005955, 6.401367);
    addThresholds( 0,  3, 13, SEC_ENDCAP, 0.006142, 6.928223);
    addThresholds( 0,  3, 14, SEC_ENDCAP, 0.009441, 8.484619);
    addThresholds( 0,  3, 15, SEC_ENDCAP, 0.015686, 10.022217);
    addThresholds( 0,  4,  5, SEC_ENDCAP, 0.001667, 0.804443);
    addThresholds( 0,  4,  6, SEC_ENDCAP, 0.003101, 0.898438);
    addThresholds( 0,  4,  7, SEC_ENDCAP, 0.004223, 1.154297);
    addThresholds( 0,  4, 12, SEC_ENDCAP, 0.005821, 6.302246);
    addThresholds( 0,  4, 13, SEC_ENDCAP, 0.004322, 7.031738);
    addThresholds( 0,  4, 14, SEC_ENDCAP, 0.006763, 8.208984);
    addThresholds( 0,  4, 15, SEC_ENDCAP, 0.010582, 10.083984);
    addThresholds( 1,  3,  4, SEC_ENDCAP, 0.001503, 0.739502);
    addThresholds( 1,  3,  5, SEC_ENDCAP, 0.001484, 0.997070);
    addThresholds( 1,  3, 11, SEC_ENDCAP, 0.004833, 5.438965);
    addThresholds( 1,  3, 12, SEC_ENDCAP, 0.005066, 6.509766);
    addThresholds( 1,  3, 13, SEC_ENDCAP, 0.005581, 6.774414);
    addThresholds( 1,  3, 14, SEC_ENDCAP, 0.006779, 8.566895);
    addThresholds( 1,  3, 15, SEC_ENDCAP, 0.008495, 6.973633);
    addThresholds( 1,  4,  5, SEC_ENDCAP, 0.000893, 0.485352);
    addThresholds( 1,  4, 12, SEC_ENDCAP, 0.004780, 6.164062);
    addThresholds( 1,  4, 13, SEC_ENDCAP, 0.004242, 7.127197);
    addThresholds( 1,  4, 14, SEC_ENDCAP, 0.006092, 7.571777);
    addThresholds( 1,  4, 15, SEC_ENDCAP, 0.006893, 6.656494);
    addThresholds( 3,  4,  5, SEC_ENDCAP, 0.001774, 1.175781);
    addThresholds( 3,  4,  6, SEC_ENDCAP, 0.002792, 2.051514);
    addThresholds( 3,  4,  7, SEC_ENDCAP, 0.004890, 3.117188);
    addThresholds( 3,  4, 12, SEC_ENDCAP, 0.006310, 6.310059);
    addThresholds( 3,  4, 13, SEC_ENDCAP, 0.004448, 7.303711);
    addThresholds( 3,  4, 14, SEC_ENDCAP, 0.005508, 8.145264);
    addThresholds( 3,  4, 15, SEC_ENDCAP, 0.006287, 11.077881);


    //// Hardware format LUT constants ////

  }
  else{
    ///// Floating point Thresholds ////

    //Barrel
    addThresholds( 0,  1,  2, SEC_BARREL, 0.002031, 0.295711);
    addThresholds( 0,  1,  8, SEC_BARREL, 0.005259, 2.692444);
    addThresholds( 0,  1,  9, SEC_BARREL, 0.010119, 2.865237);
    addThresholds( 0,  1, 10, SEC_BARREL, 0.012580, 3.116440);
    addThresholds( 0,  2,  8, SEC_BARREL, 0.002157, 2.534317);
    addThresholds( 0,  2,  9, SEC_BARREL, 0.004695, 2.579884);
    addThresholds( 0,  2, 10, SEC_BARREL, 0.007653, 2.771254);
    addThresholds( 1,  2,  8, SEC_BARREL, 0.001197, 2.562876);
    addThresholds( 1,  2,  9, SEC_BARREL, 0.004099, 2.697098);
    addThresholds( 1,  2, 10, SEC_BARREL, 0.007030, 2.893524);

    //Hybrid
    addThresholds( 0,  1,  2, SEC_HYBRID, 0.002266, 0.332948);
    addThresholds( 0,  1,  8, SEC_HYBRID, 0.004518, 2.861103);
    addThresholds( 0,  1,  9, SEC_HYBRID, 0.008576, 2.765429);
    addThresholds( 0,  1, 10, SEC_HYBRID, 0.009470, 3.028944);
    addThresholds( 0,  1, 11, SEC_HYBRID, 0.008760, 4.842227);
    addThresholds( 0,  1, 12, SEC_HYBRID, 0.009025, 5.150895);
    addThresholds( 0,  1, 13, SEC_HYBRID, 0.013822, 5.157996);
    addThresholds( 0,  2,  8, SEC_HYBRID, 0.002430, 2.535672);
    addThresholds( 0,  2,  9, SEC_HYBRID, 0.004603, 2.630091);
    addThresholds( 0,  2, 10, SEC_HYBRID, 0.004610, 2.887653);
    addThresholds( 0,  2, 11, SEC_HYBRID, 0.006391, 4.865986);
    addThresholds( 0,  2, 12, SEC_HYBRID, 0.006385, 5.008283);
    addThresholds( 0,  2, 13, SEC_HYBRID, 0.008599, 5.106528);
    addThresholds( 1,  2,  8, SEC_HYBRID, 0.002456, 2.574467);
    addThresholds( 1,  2,  9, SEC_HYBRID, 0.007859, 2.700417);
    addThresholds( 1,  2, 10, SEC_HYBRID, 0.003278, 2.848935);
    addThresholds( 1,  2, 11, SEC_HYBRID, 0.005641, 4.930833);
    addThresholds( 1,  2, 12, SEC_HYBRID, 0.005814, 5.077612);
    addThresholds( 1,  2, 13, SEC_HYBRID, 0.005388, 5.262541);

    //Endcap
    addThresholds( 0,  1,  3, SEC_ENDCAP, 0.003871, 0.455573);
    addThresholds( 0,  1,  4, SEC_ENDCAP, 0.008008, 0.610674);
    addThresholds( 0,  1,  5, SEC_ENDCAP, 0.001982, 0.502591);
    addThresholds( 0,  1, 11, SEC_ENDCAP, 0.006242, 5.456658);
    addThresholds( 0,  1, 12, SEC_ENDCAP, 0.006971, 6.008507);
    addThresholds( 0,  1, 13, SEC_ENDCAP, 0.012287, 7.060206);
    addThresholds( 0,  1, 14, SEC_ENDCAP, 0.017514, 7.984267);
    addThresholds( 0,  1, 15, SEC_ENDCAP, 0.027412, 7.601592);
    addThresholds( 0,  3,  4, SEC_ENDCAP, 0.002333, 0.790936);
    addThresholds( 0,  3,  5, SEC_ENDCAP, 0.003576, 1.233592);
    addThresholds( 0,  3,  6, SEC_ENDCAP, 0.004491, 1.925231);
    addThresholds( 0,  3,  7, SEC_ENDCAP, 0.006731, 2.384964);
    addThresholds( 0,  3, 11, SEC_ENDCAP, 0.004525, 5.417496);
    addThresholds( 0,  3, 12, SEC_ENDCAP, 0.005669, 6.401115);
    addThresholds( 0,  3, 13, SEC_ENDCAP, 0.005966, 6.933886);
    addThresholds( 0,  3, 14, SEC_ENDCAP, 0.009222, 8.496601);
    addThresholds( 0,  3, 15, SEC_ENDCAP, 0.015258, 10.029259);
    addThresholds( 0,  4,  5, SEC_ENDCAP, 0.001593, 0.806697);
    addThresholds( 0,  4,  6, SEC_ENDCAP, 0.002348, 0.894645);
    addThresholds( 0,  4,  7, SEC_ENDCAP, 0.004056, 1.173814);
    addThresholds( 0,  4, 12, SEC_ENDCAP, 0.005605, 6.305890);
    addThresholds( 0,  4, 13, SEC_ENDCAP, 0.004358, 7.032155);
    addThresholds( 0,  4, 14, SEC_ENDCAP, 0.006293, 8.210576);
    addThresholds( 0,  4, 15, SEC_ENDCAP, 0.010423, 10.089836);
    addThresholds( 1,  3,  4, SEC_ENDCAP, 0.001003, 0.734577);
    addThresholds( 1,  3,  5, SEC_ENDCAP, 0.001096, 1.003146);
    addThresholds( 1,  3, 11, SEC_ENDCAP, 0.004626, 5.441451);
    addThresholds( 1,  3, 12, SEC_ENDCAP, 0.004657, 6.512224);
    addThresholds( 1,  3, 13, SEC_ENDCAP, 0.005205, 6.785536);
    addThresholds( 1,  3, 14, SEC_ENDCAP, 0.006592, 8.608834);
    addThresholds( 1,  3, 15, SEC_ENDCAP, 0.007877, 6.975325);
    addThresholds( 1,  4,  5, SEC_ENDCAP, 0.000597, 0.484591);
    addThresholds( 1,  4, 12, SEC_ENDCAP, 0.004410, 6.162822);
    addThresholds( 1,  4, 13, SEC_ENDCAP, 0.004339, 7.126968);
    addThresholds( 1,  4, 14, SEC_ENDCAP, 0.006339, 7.569645);
    addThresholds( 1,  4, 15, SEC_ENDCAP, 0.007153, 6.650499);
    addThresholds( 3,  4,  5, SEC_ENDCAP, 0.000990, 1.185052);
    addThresholds( 3,  4,  6, SEC_ENDCAP, 0.001868, 2.053982);
    addThresholds( 3,  4,  7, SEC_ENDCAP, 0.003892, 3.117370);
    addThresholds( 3,  4, 12, SEC_ENDCAP, 0.005585, 6.314577);
    addThresholds( 3,  4, 13, SEC_ENDCAP, 0.004198, 7.300415);
    addThresholds( 3,  4, 14, SEC_ENDCAP, 0.005310, 8.142437);
    addThresholds( 3,  4, 15, SEC_ENDCAP, 0.006489, 11.085539);

    //// End of Floating point Thresholds ////
  }
}

void TCBuilder::mergePatterns(){
  //cout<<"Merging of patterns not implemented"<<endl;
}

// Tentative duplicate removal
// Not used for the moment

void TCBuilder::mergeTracks(){
}

/* From the hits of the best TC, process the track parameters, create and fill a Track object and return its pointer */
Track* TCBuilder::createFittedTrack(vector <Hit*> &bestTC)
{
  int lay;
  bool barrelmod;

  int size = bestTC.size();

  float r;
  float rmin,rmax;

  float phi_est,z_est,eta_est,pt_est;

  double parZX[2][2];
  double resZX[2];
  double invZX[2][2];
  double detZX = 0;

  double x,y,z;

  int npt=0;

  float x1,x2,y1,y2;
  float x0,y0;

  for (int i=0;i<2;++i)
    {
      resZX[i] = 0.;
      for (int j=0;j<2;++j) parZX[i][j] = 0.;
      for (int j=0;j<2;++j) invZX[i][j] = 0.;
    }

  rmin = 1000;
  rmax = 0;
  int imax=-1;
  int imin=-1;

  x0=0;
  y0=0;

  x2=0;
  y2=0;

  int n2Sb=0;
  int nPSb=0;
  int n2Se=0;
  int nPSe=0;

  // Loop over stubs in the TC
  // In order to determine the lowest/largest radius of the 
  // TC, and the TC composition
  //
  // The stub with the largest radius is our reference in 
  // the conformal space


  for (int i=0;i<size;++i) 
    {
      lay  = bestTC.at(i)->getLayer();

      if (lay<=7) //If PS
        {
          (lay<=2) ? nPSb++ : nPSe++;
        }
      else  //else 2S
        {
          (lay<11) ? n2Sb++ : n2Se++;
        }

      if (lay>10) continue; // 2S disks stubs not used
      ++npt;

      x = bestTC.at(i)->getX();
      y = bestTC.at(i)->getY();
      r = sqrt(x*x+y*y);

      if (r<rmin)
        {
          rmin = r;
          imin = i;
          x2   = x;
          y2   = y;
        }

      if (r>rmax)
        {
          rmax = r;
          imax = i;
          x0   = x;
          y0   = y;
        }
    }

  float rmax2 = 0;

  x1 = 0;
  y1 = 0;

  int nc=0;

  // Loop 2 over stubs in the TC
  // In order to determine the point with the second largest radius 

  for (int i=0;i<size;++i) // Loop over stubs in the TC
    {
      if (i==imax || i==imin) continue; // Already used

      lay  = bestTC.at(i)->getLayer();

      barrelmod=0;
      if (lay<=2 || (lay>=8 && lay<=10)) barrelmod=1;
      if (!barrelmod && (nPSb+n2Sb)>=3) continue; // Barrel modules have a priority
      if (lay>10 && (nPSb+nPSe)>=3) continue;     // Don't use 2S modules in the disks if possible

      x = bestTC.at(i)->getX();
      y = bestTC.at(i)->getY();
      r = sqrt(x*x+y*y);

      if (r>rmax2)
        {
          rmax2  = r;
          x1     = x;
          y1     = y;
          nc++;
        }
    }

  // Now get the coordinates in the conformal space.

  double sqR1  = binning((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0), 13, 18, UNSIGNED);
  double sqR2  = binning((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0), 13, 18, UNSIGNED);

  x1 = binning((x1-x0)/sqR1, -3, 18, SIGNED);
  y1 = binning((y1-y0)/sqR1, -3, 18, SIGNED);

  x2 = binning((x2-x0)/sqR2, -3, 18, SIGNED);
  y2 = binning((y2-y0)/sqR2, -3, 18, SIGNED);


  double mult_x1_y2 = binning(x1*y2, -8, 18, SIGNED);
  double mult_x2_y1 = binning(x2*y1, -8, 18, SIGNED);

  // Now we got everything for the r/phi plane   

  double sub_of_mult = binning(mult_x1_y2 - mult_x2_y1, -11, 18, SIGNED);

  double divA = binning((y1-y2)/sub_of_mult, 17, 18, SIGNED);
  double divB = binning((x2-x1)/sub_of_mult, 17, 18, SIGNED);

  double a = binning(x0-0.5*divA, 17, 18, SIGNED);
  double b = binning(y0-0.5*divB, 17, 18, SIGNED);

  double a2_b2 = binning(a*a + b*b, 29, 18, UNSIGNED);

  pt_est  = 0.003*3.833*sqrt(a2_b2);

  int charge =-b/fabs(b);

  phi_est = atan2(charge*a,-charge*b);

  //Apply the rotation
  phi_est += sec_phi;

  //Set the value betweend -Pi and Pi
  phi_est = fmod(phi_est + M_PI, 2 * M_PI) - M_PI;

  // Then we do the RZ fit (LS)

  float wght;
  int cnt=0;
  for (int i=0;i<size;++i) // Loop over stubs in the TC
    {
      lay  = bestTC.at(i)->getLayer();

      if (lay>7) continue; // Don't use 2S modules
      if (lay>2 && nPSb>=2) continue; // Don't use PS modules of the disks if number of good point in the barrel is sufficient        

      ++cnt;
      x = bestTC.at(i)->getX();
      y = bestTC.at(i)->getY();
      z = bestTC.at(i)->getZ();        
      r = sqrt(x*x+y*y);

      wght=1;
      if (lay>2) wght=1/7.;

      parZX[0][0] += wght*r*r;
      parZX[1][1] += wght*1;
      parZX[1][0] += wght*r;

      resZX[0] += wght*r*z;
      resZX[1] += wght*z;       
    } // End of stub loop
    

  detZX = parZX[0][0]*parZX[1][1]-parZX[1][0]*parZX[1][0];

  invZX[0][0] =  parZX[1][1]/detZX;
  invZX[1][0] = -parZX[1][0]/detZX;
  invZX[1][1] =  parZX[0][0]/detZX;

  // Finally estimate track parameters in the R/Z plane 

  eta_est = std::asinh((invZX[0][0]*resZX[0] + invZX[1][0]*resZX[1]));
  z_est   = invZX[1][0]*resZX[0] + invZX[1][1]*resZX[1];

  Track* fit_track = new Track();
  fit_track->setCurve(pt_est);
  fit_track->setPhi0(phi_est);
  fit_track->setEta0(eta_est);
  fit_track->setZ0(z_est);
      
  for(unsigned int hitIndex=0;hitIndex < bestTC.size();hitIndex++)
    {
      fit_track->addStubIndex(bestTC[hitIndex]->getID());
    }

  return fit_track;
}

/* Function which simulate the HardWare representation of the values : manage UNSIGNED and SIGNED (2's complement) overflows and accuracy according to the available dynamic of the binary word */
double TCBuilder::binning(double fNumber, int nMSBpowOfTwo, int nBits, HW_SIGN_TYPE signType)
{
  if (!m_bHardwareSimulation)
    //If the Hardware binning simulation is not asked, return directly the original number
    return fNumber;

  if (signType == UNSIGNED && fNumber < 0)
    {
		  //Bad interpretation, a negative number is stored in an UNSIGNED format (sign lost)
		  fNumber = -fNumber;
    }
  
  int nLSBpowOfTwo;
	
	//Process the power of two of the LSB for the binary representation
	if (signType == UNSIGNED)
    {
      //If UNSIGNED
		  nLSBpowOfTwo = nMSBpowOfTwo - (nBits-1);
    }	
  else
    {
      //If SIGNED, 1 bit is used for the sign
      nLSBpowOfTwo = nMSBpowOfTwo - (nBits-2);
		}

  /* Accuracy Simulation */

  //Divide the number by the power of two of the LSB => the integer part of the new number is the value we are looking for
	fNumber = fNumber / pow(2, nLSBpowOfTwo);
	
  //Remove the fractionnal part by rounding down (for both positive and negative values), this simulate the HW truncature
	fNumber = floor(fNumber);
	
	//Multiply the number by the power of two of the LSB to get the correct float value
	fNumber = fNumber * pow(2, nLSBpowOfTwo);


  double fBinnedNumber = fNumber;

	/* Overflow Simulation */

	if (signType == UNSIGNED)
		{
		  //If the number is in UNSIGNED representation
			fNumber = fmod(fNumber, pow(2, nMSBpowOfTwo+1));
    }
	else
    {
		  //If the number is in SIGNED representation (2's complement)
      
      double fTempResult = fNumber - pow(2, nMSBpowOfTwo+1); //substract the possible range to the number

      if (fTempResult >= 0)
        {
          //If there is an overflow, it's a positive one
          fNumber = fmod(fTempResult, pow(2, nMSBpowOfTwo+2)) - pow(2, nMSBpowOfTwo+1);
        }
      else
        {
          //If there is an overflow, it's a negative one (2's complement format has an asymetric range for positive and negative values)
          fNumber = fmod(fTempResult + pow(2, nLSBpowOfTwo), pow(2, nMSBpowOfTwo+2)) - pow(2, nLSBpowOfTwo) + pow(2, nMSBpowOfTwo+1);
        }
		}

  //If the new number is different from the previous one, an HW overflow occured
  if (fNumber != fBinnedNumber)
    {
      cout<<"WARNING HW overflow for the value : "<<fBinnedNumber<<" resulting value : "<<fNumber<<" (diff= "<<fBinnedNumber-fNumber<<")"<<endl;
    }
	
	return fNumber;
}

/* Process the alignment scores (on RPHI and on RZ plan) between the 2 seeds and an other stub */
void TCBuilder::alignScore(Hit& hSeed1, Hit& hSeed2, Hit& hTestStub, double tScores[])
{
  double fRPHI_Score , fRZ_Score;

  double X1, Y1, Z1, R1, PHI1;
  double X2, Y2, Z2, R2, PHI2;
  double X3, Y3, Z3, R3, PHI3;

  double RPHI_S1, RPHI_S2, RZ_S1, RZ_S2;

  //Coordinates X, Y, Z are already binned
  X1 = hSeed1.getX();
  Y1 = hSeed1.getY();
  Z1 = hSeed1.getZ();
	
  X2 = hSeed2.getX();
  Y2 = hSeed2.getY();
  Z2 = hSeed2.getZ();

  X3 = hTestStub.getX();
  Y3 = hTestStub.getY();
  Z3 = hTestStub.getZ();
	
  R1 = binning(sqrt(X1*X1 + Y1*Y1), 6, 18, SIGNED);
  R2 = binning(sqrt(X2*X2 + Y2*Y2), 6, 18, SIGNED);
  R3 = binning(sqrt(X3*X3 + Y3*Y3), 6, 18, SIGNED);

  //RPHI plan
  PHI1 = binning(atan(Y1/X1), 4, 18, SIGNED);
  PHI2 = binning(atan(Y2/X2), 4, 18, SIGNED);
  PHI3 = binning(atan(Y3/X3), 4, 18, SIGNED);

  RPHI_S1 = binning((PHI2 - PHI1) * (R3 - R2), 8, 20, SIGNED);
  RPHI_S2 = binning((PHI2 - PHI3) * (R2 - R1), 8, 20, SIGNED);

  fRPHI_Score = binning(fabs(RPHI_S1 + RPHI_S2), 7, 18, UNSIGNED);

  //RZ plan
  RZ_S1 = binning((Z2 - Z1) * (R3 - R2), 12, 20, SIGNED);
  RZ_S2 = binning((Z2 - Z3) * (R2 - R1), 12, 20, SIGNED);

  fRZ_Score = binning(fabs(RZ_S1 + RZ_S2), 10, 18, UNSIGNED);

  tScores[0] = fRPHI_Score;
  tScores[1] = fRZ_Score;
}


/* Fill thresholds data with new thresholds */
void TCBuilder::addThresholds(int nLaySeed1, int nLaySeed2, int nLayTestStub, SEC_TYPE secType, double fRPHI_Thresh, double fRZ_Thresh)
{

  //Fill the tab corresponding to the sector type with the pair of thresholds
  switch (secType)
    {
    case SEC_BARREL : m_tabBarrelThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0] = fRPHI_Thresh;
      m_tabBarrelThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1] = fRZ_Thresh;
      break;
    case SEC_HYBRID : m_tabHybridThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0] = fRPHI_Thresh;
      m_tabHybridThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1] = fRZ_Thresh;
      break;
    case SEC_ENDCAP : m_tabEndcapThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0] = fRPHI_Thresh;
      m_tabEndcapThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1] = fRZ_Thresh;
    }
}

/* Get the thresholds corresponding to the 3 layerID for a given sector type */
void TCBuilder::getThresholds(int nLaySeed1, int nLaySeed2, int nLayTestStub, SEC_TYPE secType, double tabThresh[])
{
  switch (secType)
    {
    case SEC_BARREL : tabThresh[0] = m_tabBarrelThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0];
      tabThresh[1] = m_tabBarrelThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1];
      break;
    case SEC_HYBRID : tabThresh[0] = m_tabHybridThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0];
      tabThresh[1] = m_tabHybridThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1];
      break;
    case SEC_ENDCAP : tabThresh[0] = m_tabEndcapThresholds[nLaySeed1][nLaySeed2][nLayTestStub][0];
      tabThresh[1] = m_tabEndcapThresholds[nLaySeed1][nLaySeed2][nLayTestStub][1];
    }
}


/* Operate the layerID transcoding (to get the layerID 4bits Hardawre representation) */
char TCBuilder::transcodeLayer(Hit * pHit)
{
  int nOrigLayer = pHit->getLayer();

  //layer transcoding of the disks is based on the radius, it can be optimized by using the ladder_ID
  double X = pHit->getX();
  double Y = pHit->getY();

  int nTransLayer;
		
  if (nOrigLayer <= 10)
    {
      //If the stub is in the barrel

      if (nOrigLayer <= 7)
	{    
	  //If layer 5, 6, 7
	  nTransLayer = nOrigLayer - 5;     //5->0, 6->1, ...
	}    
      else
	{
	  //If layer 8, 9, 10
	  nTransLayer = nOrigLayer;         //no change
	}
    }
  else if (sqrt(X*X + Y*Y) >= 62)
    {
      //If the stub is on an outer ring of a disk !!! (2S modules)

      if (nOrigLayer <= 15)
	{
	  //If layer 11, 12, 13, 14, 15
	  nTransLayer = nOrigLayer;         //no change
	}
      else
	{
	  //If layer 18, 19, 20, 21, 22
	  nTransLayer = nOrigLayer - 7;     //18->11, 19->12, ...
	}
    }
  else
    {
      //If the stub is on an inner ring of a disk !!! (PS modules)

      if (nOrigLayer <= 15)
	{
	  //If layer 11, 12, 13, 14, 15
	  nTransLayer = nOrigLayer - 8;     //11->3, 12->4, ...
	}
      else
	{
	  //If layer 18, 19, 20, 21, 22
	  nTransLayer = nOrigLayer - 15;    //18->3, 19->4, ...
	}
    }

  return nTransLayer;	
}

// TC builder module
/* Take as input the list of stubs contained in a matched road */

void TCBuilder::fit(vector<Hit*> originalHits)
{

  //cout<<"trying to fit "<<originalHits.size()<<" points"<<endl;

  int tow = sector_id; // The tower ID, necessary to get the phi shift

  SEC_TYPE currentSec;

  //Get the sec_type from the sector_id
  if (tow >= 16 && tow <= 31)
    currentSec = SEC_BARREL; 
  else if (tow >=8 && tow <=39)
    currentSec = SEC_HYBRID;
  else
    currentSec = SEC_ENDCAP;

  //Process the starting phi of the tower
  sec_phi = (tow%8) * M_PI / 4.0 - 0.4;

  //cos and sin values for a rotation of an angle -sec_phi
  double ci = cos(-sec_phi);
  double si = sin(-sec_phi);

  double rotatedX, rotatedY;

  //Create a new vector to store the custom hits parameters
  vector <Hit> hits;
  
  Hit* pOrigHit;
  
  //For each hit of the lists
  for (unsigned int origHitIndex = 0; origHitIndex<originalHits.size(); origHitIndex++)
    {
      pOrigHit = originalHits[origHitIndex];

      //Process the rotated coordinnates
      rotatedX = pOrigHit->getX() * ci - pOrigHit->getY() * si;
      rotatedY = pOrigHit->getX() * si + pOrigHit->getY() * ci;

      //Add the modified hit to the hits vector
      hits.push_back( Hit(transcodeLayer(pOrigHit),
			  pOrigHit->getLadder(),
			  pOrigHit->getModule(),
			  pOrigHit->getSegment(),
			  pOrigHit->getStripNumber(),
			  pOrigHit->getID(),
			  pOrigHit->getParticuleID(),
			  pOrigHit->getParticulePT(),
			  pOrigHit->getParticuleIP(),
			  pOrigHit->getParticuleETA(),
			  pOrigHit->getParticulePHI0(),
			  binning(rotatedX, 6, 18, SIGNED),
			  binning(rotatedY, 6, 18, SIGNED),
			  binning(pOrigHit->getZ(), 8, 18, SIGNED),
			  pOrigHit->getX0(),
			  pOrigHit->getY0(),
			  pOrigHit->getZ0(),
			  pOrigHit->getBend())
		      );

    }

  //Sort the hits by ascending order of layerID
  //(using a lambda definition of the sorting criteria which return a boolean)
  sort(hits.begin(), hits.end(), [ ]( const Hit& lhs, const Hit& rhs ) { return lhs.getLayer() < rhs.getLayer(); });


  int nLayersCurrentPattern = 0;
  int lastAddedLayer = -1;
  //Count the number of layers present in the pattern
  for (unsigned int hitIndex=0; hitIndex < hits.size(); hitIndex++)
    {
      if (lastAddedLayer != hits[hitIndex].getLayer())
	{
	  nLayersCurrentPattern++;
	  lastAddedLayer = hits[hitIndex].getLayer();
	}
    }


  vector <Hit*> vecCurrentCandidateHits;
  vector <double> vecCurrentCandidateScore;

  vector <Hit*> vecBestCandidateHits;
  double fBestCandidateScore = 0.0;

  for (unsigned int seed1Index=0; seed1Index<hits.size(); seed1Index++)
    {
      Hit& hSeed1   = hits[seed1Index];
      int nLaySeed1 = hSeed1.getLayer();

      if (nLaySeed1 == 2) continue; //layer 2 can't be the innermost seed stub
      if (nLaySeed1 > 3) break;     //no more possible combinations for this pattern

      //We have a correct Seed1

      //Get the radius of the seed1
      double fRseed1 = binning(sqrt(hSeed1.getX()*hSeed1.getX() + hSeed1.getY()*hSeed1.getY()), 6, 18, SIGNED);

      for (unsigned int seed2Index = seed1Index+1; seed2Index<hits.size(); seed2Index++)
	{
	  Hit& hSeed2   = hits[seed2Index];
	  int nLaySeed2 = hSeed2.getLayer();

	  if (nLaySeed1 == nLaySeed2) continue; //The seed layers have to be differents
	  if (nLaySeed2 > 4) break;             //no more possible combinations for the current seed1


	  //We have a correct Seed1/Seed2 combination !!!

    //Get the radius of the seed2
    double fRseed2 = binning(sqrt(hSeed2.getX()*hSeed2.getX() + hSeed2.getY()*hSeed2.getY()), 6, 18, SIGNED);

	  //Current candidate initialization (the 2 seeds)
	  vecCurrentCandidateHits.clear();
	  vecCurrentCandidateHits.push_back(&hSeed1);
	  vecCurrentCandidateHits.push_back(&hSeed2);

	  vecCurrentCandidateScore.clear();


	  for (unsigned int testStubIndex = seed2Index+1; testStubIndex<hits.size(); testStubIndex++)
	    {

	      Hit& hTestStub    = hits[testStubIndex];
	      int nLayTestStub  = hTestStub.getLayer();

	      if (nLayTestStub == nLaySeed2) continue; //The layers have to be differents

        
	      //Score processing of the Seed1/Seed2/testStub combination
	      double tabScore[2];
	      alignScore(hSeed1, hSeed2, hTestStub, tabScore);
        

	      //Get the thresholds corresponding to the current layer combination
	      double tabNormThresh[2];
	      getThresholds(nLaySeed1, nLaySeed2, nLayTestStub, currentSec, tabNormThresh);


        //Process the real thresholds from the normalized one stored in the LUT
        double fThreshRPHI = binning(fabs(tabNormThresh[0] * (fRseed2 - fRseed1)), 4, 18, SIGNED);
        double fThreshRZ = binning(fabs(tabNormThresh[1] * (fRseed2 - fRseed1)), 8, 18, SIGNED);

	      if (tabScore[0] <= fThreshRPHI && tabScore[1] <= fThreshRZ)
		{
		  //The stub is in the window defined by the seed projection (correct stub candidate !)
          
		  if (nLayTestStub != vecCurrentCandidateHits.back()->getLayer())
		    {
		      //The current testStub layer is not yet in the TC
		      vecCurrentCandidateHits.push_back(&hTestStub);
		      vecCurrentCandidateScore.push_back(tabScore[0]);
		    }
		  else if (tabScore[0] < vecCurrentCandidateScore.back())
		    {
		      //The layer is already in the TC but the Phi score of the current stub is better than the previous one
		      vecCurrentCandidateHits.back()   = &hTestStub;            
		      vecCurrentCandidateScore.back()  = tabScore[0];
		    }
		}
	    }

	  //If the current candidate own more than 6 stubs, the lasts (outtermost) are removed
	  while (vecCurrentCandidateHits.size() > 6)
	    {
	      vecCurrentCandidateHits.pop_back();
	      vecCurrentCandidateScore.pop_back();
	    }

	  //All the stubs have been tested for the current Seeds combination

	  if (int(vecCurrentCandidateHits.size()) >= nLayersCurrentPattern - m_nMissingHits)
	    {
	      //The current candidate has enough stubs to be a candidate
 
	      //Process the score of the track candidate
	      double fCurrentCandidateScore = 0.0;
	      while (vecCurrentCandidateScore.empty() == false)
		{
		  fCurrentCandidateScore += vecCurrentCandidateScore.back();
		  vecCurrentCandidateScore.pop_back();
		}

	      if (vecCurrentCandidateHits.size() > vecBestCandidateHits.size() || (vecCurrentCandidateHits.size() == vecBestCandidateHits.size() && fCurrentCandidateScore < fBestCandidateScore))
		{
		  //The current candidate is better than the previous best one
		  vecBestCandidateHits = vecCurrentCandidateHits;
		  fBestCandidateScore = fCurrentCandidateScore;
		}
	    }
	}
    }

  //All the Seeds combinations have been tested

  if ( (currentSec == SEC_HYBRID && vecBestCandidateHits.size() >= 4) || vecBestCandidateHits.size() >= 5 )
    {
      //If there is a recorded best candidate

      //Fit the parameters and create the corresponding track object
      Track * fit_track;
      fit_track = createFittedTrack(vecBestCandidateHits);
    
      //cout<<"adding one track..."<<endl;
      tracks.push_back(fit_track);
    }
}

void TCBuilder::fit(){
  for(unsigned int i=0;i<patterns.size();i++){
    vector<Hit*> allHits = patterns[i]->getHits();
    fit(allHits);
  }
}

TrackFitter* TCBuilder::clone(){
  TCBuilder* fit = new TCBuilder(nb_layers);
  fit->setPhiRotation(sec_phi);
  fit->setSectorID(sector_id);
  return fit;
}

