/*
C++ implementation of the TC builder

S.Viret, G.Galbit : 09/12/15
*/


#include "../interface/TCBuilder.h"

TCBuilder::TCBuilder():TrackFitter(0){

}

TCBuilder::TCBuilder(int nb):TrackFitter(nb)
{

  //Layer infos
  m_nLayer = 6;

  //
  // Initialisation of the matching windows
  //

  for (int i=0;i<3;++i)
  {
    for (int j=0;j<6;++j)
    {
      fenetre_b_z[i][j]=-1;
      fenetre_b_phi[i][j]=-1;
    }    

    for (int j=0;j<15;++j)
    {
      fenetre_e_r[i][j]=-1;
      fenetre_e_phi[i][j]=-1;
    } 
  }

  // Windows for the barrel towers
  // 1 z window per layer (in cm)
  // 1 phi window per layer (in rad)

  fenetre_b_z[2][0]  = 0.28;  
  fenetre_b_z[2][1]  = 0.14;  
  fenetre_b_z[2][2]  = 0.34;  
  fenetre_b_z[2][3]  = 2.8;  
  fenetre_b_z[2][4]  = 3.0;  
  fenetre_b_z[2][5]  = 3.2;  
  fenetre_b_phi[2][0]= 0.0026;  
  fenetre_b_phi[2][1]= 0.0013;  
  fenetre_b_phi[2][2]= 0.003;  
  fenetre_b_phi[2][3]= 0.0065;  
  fenetre_b_phi[2][4]= 0.0115;  
  fenetre_b_phi[2][5]= 0.016;  
   
  // Windows for the endcap towers
  // For barrel layers
  // 1 z window per layer (in cm)
  // 1 phi window per layer (in rad)

  fenetre_b_z[0][0]  = 1.6;  
  fenetre_b_z[0][1]  = 0.9;  
  fenetre_b_phi[0][0]= 0.0038;  
  fenetre_b_phi[0][1]= 0.0016;  
 
  // For disks
  // 1 r window per ring (in cm)
  // 1 phi window per ring (in rad)
  fenetre_e_r[0][0]  = 0.1;  
  fenetre_e_r[0][1]  = 0.25;  
  fenetre_e_r[0][2]  = 0.33;  
  fenetre_e_r[0][3]  = 0.38;  
  fenetre_e_r[0][4]  = 0.46;  
  fenetre_e_r[0][5]  = 0.5;  
  fenetre_e_r[0][6]  = 0.65;  
  fenetre_e_r[0][7]  = 0.7;  
  fenetre_e_r[0][8]  = 0.75;
  fenetre_e_r[0][9]  = 2.9;
  fenetre_e_r[0][10] = 2.8;
  fenetre_e_r[0][11] = 3.0;
  fenetre_e_r[0][12] = 3.0;
  fenetre_e_r[0][13] = 3.1;
  fenetre_e_r[0][14] = 3.2;
  fenetre_e_phi[0][0]= 0.0006;  
  fenetre_e_phi[0][1]= 0.0018;  
  fenetre_e_phi[0][2]= 0.0025;  
  fenetre_e_phi[0][3]= 0.0038;  
  fenetre_e_phi[0][4]= 0.005;  
  fenetre_e_phi[0][5]= 0.0055;  
  fenetre_e_phi[0][6]= 0.007;  
  fenetre_e_phi[0][7]= 0.008;  
  fenetre_e_phi[0][8]= 0.008; 
  fenetre_e_phi[0][9]= 0.010; 
  fenetre_e_phi[0][10]= 0.010; 
  fenetre_e_phi[0][11]= 0.012; 
  fenetre_e_phi[0][12]= 0.012; 
  fenetre_e_phi[0][13]= 0.015; 
  fenetre_e_phi[0][14]= 0.017; 

 // Windows for the hybrid towers

  fenetre_b_z[1][0]  = 0.28;  
  fenetre_b_z[1][1]  = 0.14;  
  fenetre_b_z[1][2]  = 0.34;  
  fenetre_b_z[1][3]  = 2.8;  
  fenetre_b_z[1][4]  = 3.0;  
  fenetre_b_z[1][5]  = 3.2;  
  fenetre_b_phi[1][0]= 0.0026;  
  fenetre_b_phi[1][1]= 0.0013;  
  fenetre_b_phi[1][2]= 0.003;  
  fenetre_b_phi[1][3]= 0.0065;  
  fenetre_b_phi[1][4]= 0.0115;  
  fenetre_b_phi[1][5]= 0.016; 
  
  fenetre_e_r[1][3]  = 0.08;  
  fenetre_e_r[1][4]  = 0.09;  
  fenetre_e_r[1][5]  = 0.11;  
  fenetre_e_r[1][6]  = 0.12;  
  fenetre_e_r[1][7]  = 0.14;  
  fenetre_e_r[1][8]  = 0.16;
  fenetre_e_r[1][9]  = 2.45;
  fenetre_e_r[1][10] = 2.6;
  fenetre_e_r[1][11] = 2.65;
  fenetre_e_r[1][12] = 2.7;
  fenetre_e_r[1][13] = 2.8;
  fenetre_e_r[1][14] = 3.0;
  fenetre_e_phi[1][3]= 0.00065;  
  fenetre_e_phi[1][4]= 0.00095;  
  fenetre_e_phi[1][5]= 0.0018;  
  fenetre_e_phi[1][6]= 0.0022;  
  fenetre_e_phi[1][7]= 0.0033;  
  fenetre_e_phi[1][8]= 0.0036; 
  fenetre_e_phi[1][9]= 0.007; 
  fenetre_e_phi[1][10]= 0.0095; 
  fenetre_e_phi[1][11]= 0.012; 
  fenetre_e_phi[1][12]= 0.012; 
  fenetre_e_phi[1][13]= 0.014; 
  fenetre_e_phi[1][14]= 0.016; 

  // cout<<"TC builder init done"<<endl;
}

TCBuilder::~TCBuilder()
{
}

void TCBuilder::initialize(){

}

void TCBuilder::mergePatterns(){
  //cout<<"Merging of patterns not implemented"<<endl;
}

// Tentative duplicate removal
// Not used for the moment

void TCBuilder::mergeTracks(){
  unsigned int index = 0;
  vector<Track*>::iterator it = tracks.begin();
  while(it!=tracks.end()){
    Track* newTrack = *it;
    bool found = false;
    for(unsigned int i=0;i<index;i++){
      Track* ref = tracks[i];
      float dpt,dphi,dz,deta;
      dpt = fabs(newTrack->getCurve()-ref->getCurve());
      dphi = fabs(newTrack->getPhi0()-ref->getPhi0());
      dz = fabs(newTrack->getZ0()-ref->getZ0());
      deta = fabs(newTrack->getEta0()-ref->getEta0());
      found = (deta<0.02) &&
	(dphi<0.005) &&
	(dpt<0.1) &&
	(dz<0.3);
      if(found)
	break;
    }
    if(found)
      tracks.erase(it);
    else{
      index++;
      it++;
    }
  }
}

// TC builder module
//
// Take as input the list of stubs contained in a matched road

void TCBuilder::fit(vector<Hit*> hits)
{
  //  cout << "Into the TC builder " << endl;

  int tow = sector_id; // The tower ID, necessary to get the phi shift

  //  cout << "Into the TC builder " << tow << endl;

  float sec_phi = 0;
  if (tow%8==0) sec_phi = 0.4; 
  if (tow%8==1) sec_phi = 1.2; 
  if (tow%8==2) sec_phi = 2.0; 
  if (tow%8==3) sec_phi = 2.7; 
  if (tow%8==4) sec_phi = 3.6; 
  if (tow%8==5) sec_phi = -2.0; 
  if (tow%8==6) sec_phi = -1.2; 
  if (tow%8==7) sec_phi = -0.4; 

  float ci = cos(sec_phi);
  float si = sin(sec_phi);


  int type=1; // Hybrid per default

  if (tow>15 && tow<32) type = 2; // Barrel
  if (tow<8 || tow>39)  type = 0; // Endcap

  idxs_lay.clear();
  list_cand6.clear();
  list_scores6.clear();
  list_cand5.clear();   
  list_scores5.clear();

  std::vector<int> slist, cand;
    
  int s_k,s_l;
  int lad_ref,lad_klm;  
  int idx_j,idx_i,idx_kl;

  float a_rp,b_rp,a_rz,b_rz;  

  float ri,zi,pi,xi,yi;
  float rj,zj,pj,xj,yj;
  float rkl,zkl,pkl,xkl,ykl;

  float dist_p,dist_z,dist_r;
  float dist_pm,dist_rm,dist_zm,idx_klm;

  float score;  
  int nseeds=0;

  for (int k=0;k<16;++k)
  {
    slist.clear();
    slist.push_back(-1);
    idxs_lay.push_back(slist);
  }
   
  // First of all we 
  // order the stubs per layers
  // and count the number of layers touched    

  // Layers are numbered as follows
  // Barrel      : 0,1,2,8,9,10
  // Disk z+/- PS: 3,4,5,6,7
  // Disk z+/- 2S: 11,12,13,14,15

  int lay;
  int patt_size=0;
  bool barrelmod;

  for (unsigned int k=0;k<hits.size();++k) // Loop over stubs in the pattern
  {
    lay  = hits[k]->getLayer();
    
    if (idxs_lay.at(lay).size()==1) ++patt_size;
    
    idxs_lay.at(lay).push_back(k);
  }
    

  if (patt_size>6) patt_size=6; // Max TC_size.

  // Now test the seeds
  // Not all the layers are used of the seeds

  // For barrel and hybrid, layers 0,1,2

  // For endcaps, layers 0/1/3/4    
  // In disk, only PS modules are used

  for (int k=0;k<5;++k) // Seed layer 1
  {
    if (idxs_lay.at(k).size()==1) continue; // No stub in this layer
    if (type>0  && k>2)           continue; // Not a good barrel/hybrid seed
    if (type==0 && k==2)          continue; // Not a good endcap seed
    
    for (int l=k+1;l<5;++l) // Seed layer 2
    {
      if (idxs_lay.at(l).size()==1) continue; // No stub in this layer
      if (type>0  && l>2)           continue; // Not a good barrel/hybrid seed
      if (type==0 && l==2)          continue; // Not a good endcap seed

      //      std::cout << "Seeds for tow type " << type << ": " << k << "/" << l << endl;
            
      s_k = idxs_lay.at(k).size(); // Number of stubs on seed layer 1
      s_l = idxs_lay.at(l).size(); // Number of stubs on seed layer 2

      // Now test all possible seeds            
      for (int i=1;i<s_k;++i)
      {
	for (int j=1;j<s_l;++j)
	{
	  ++nseeds; // We have a seed!

	  idx_i = idxs_lay.at(k).at(i); // Index of stub 1
	  idx_j = idxs_lay.at(l).at(j); // Index of stub 2
                 
	  xi = hits[idx_i]->getX()*ci+ hits[idx_i]->getY()*si;
	  yi = -hits[idx_i]->getX()*si+ hits[idx_i]->getY()*ci;
	  zi = hits[idx_i]->getZ();
	  ri = sqrt(xi*xi+yi*yi);
	  pi = atan2(yi,xi);

	  xj = hits[idx_j]->getX()*ci+ hits[idx_j]->getY()*si;
	  yj = -hits[idx_j]->getX()*si+ hits[idx_j]->getY()*ci;
	  zj = hits[idx_j]->getZ();
	  rj = sqrt(xj*xj+yj*yj);
	  pj = atan2(yj,xj);
	
	  // We now have seed coordinates
						   						   
	  if (pi==pj || zi==zj) continue; // Not done for the moment
                    
	  a_rp = (ri-rj)/(pi-pj);
	  b_rp = (ri*pj-rj*pi)/(pj-pi);
	  a_rz = (ri-rj)/(zi-zj);
	  b_rz = (ri*zj-rj*zi)/(zj-zi);
                    
	  // And tracklet direction in r/phi/z space                    
	  // Try to extend the seed
                
	  cand.clear();
	  score=0;
                    
	  // The 2 seeding stubs are pushed into the pre-TC
	  cand.push_back(idx_i);
	  cand.push_back(idx_j);
                    
	  for (int kk=0;kk<16;++kk) // Loop over all the layers
	  {
	    if (kk==k || kk==l) continue; // Seeding layers are excluded
	    if (idxs_lay.at(kk).size()==1) continue; // No stubs on layer
                        
	    // Select the best candidate in the given layer, and increment
	    
	    dist_pm = 10000;
	    dist_zm = 10000;
	    dist_rm = 10000;
	    idx_klm = -1;
	    lad_klm = -1;
                    
	    barrelmod=0;
	    if (kk<=2 || (kk>=8 && kk<=10)) barrelmod=1;
    
	    // Loop over all stubs contained in the test layer
	    for (unsigned int kl=1;kl<idxs_lay.at(kk).size();++kl)
	    {

	      idx_kl = idxs_lay.at(kk).at(kl);
	      lad_ref = hits[idx_kl]->getLadder();

	      // Get stub coordinates

	      xkl = hits[idx_kl]->getX()*ci+ hits[idx_kl]->getY()*si;
	      ykl = -hits[idx_kl]->getX()*si+ hits[idx_kl]->getY()*ci;
	      zkl = hits[idx_kl]->getZ();
	      rkl = sqrt(xkl*xkl+ykl*ykl);
	      pkl = atan2(ykl,xkl);
    
	      // Get residuals w.r.t. the seed projection

	      dist_p = fabs((rkl - b_rp)/a_rp - pkl);

	      if (dist_p>dist_pm) continue;	      

	      dist_z = fabs((rkl - b_rz)/a_rz - zkl);
	      dist_r = fabs(b_rz+a_rz*zkl - rkl);

	      // Check that the stub is in the extrapolation windows
	      // And that it is closer to the tracklet projection 
	      // if there is another stub already within the windows

	      if (barrelmod) // Barrel
	      {
		lay=kk;
		if(kk>=8) lay=kk-5;

		if (dist_p>fenetre_b_phi[type][lay]) continue;
		if (dist_z>fenetre_b_z[type][lay]) continue;
	      }
	      else
	      {
		if (dist_p>fenetre_e_phi[type][lad_ref]) continue;
		if (dist_r>fenetre_e_r[type][lad_ref]) continue;
	      }

	      // This stub is a good candidate, keep it for the moment
                            
	      dist_pm = dist_p;
	      dist_zm = dist_z;
	      dist_rm = dist_r;
	      
	      idx_klm = idx_kl;
	      lad_klm = lad_ref;
	    }
                        
	    if (idx_klm==-1) continue; // No stub in the window, skip;
                        
	    // Add the best candidate to the pre-TC

	    cand.push_back(idx_klm);
			
	    // Increment the TC score (will be used to select among different TCs)

	    if (barrelmod) // Barrel
	    {
	      lay=kk;
	      if(kk>=8) lay=kk-5;
	      
	      score += dist_pm/fenetre_b_phi[type][lay]+dist_zm/fenetre_b_z[type][lay];
	    }
	    else
	    {
	      score += dist_pm/fenetre_e_phi[type][lad_klm]+dist_rm/fenetre_e_r[type][lad_klm];
	    }
	  }
                    
	  // End of the loop, do we have a respectable TC or not?
                    
	  if (patt_size-cand.size()>1) continue; // For a road of size N, TC must have a size of at least N-1 
	  if(type==0 && cand.size()<5) continue;
	  if(type==1 && cand.size()<4) continue;
	  if(type==2 && cand.size()<5) continue;
          
	  // OK, pass it
          
	  if (patt_size-cand.size()==1) // One missed layer
	  {
	    list_cand5.push_back(cand);
	    list_scores5.push_back(score/cand.size());
	  }
	  else // No missed layer
	  {
	    list_cand6.push_back(cand);
	    list_scores6.push_back(score/cand.size());
	  }
	}
      }
    }
  }

  // End of the loop over all seeds
    
  //cout <<  list_cand5.size() << "///" <<  list_cand6.size() << endl; 

  // OK, end of the loop select the best TC among all the TCs built.

  if (list_cand6.size()==0 && list_cand5.size()==0) return;

  std::vector<Hit*> bestTC;

  bestTC.clear();

  // First look at 6 stubs candidates
    
  float score_min = 10000000000000;
  int idx_cand    = -1;
  int ids;

  if (list_cand6.size()!=0)
  {
    for (unsigned int i=0;i<list_scores6.size();++i)
    {
      if (list_scores6.at(i)<score_min)
      {
	score_min = list_scores6.at(i);
	idx_cand  = i;
      }
    } 
        
    for (unsigned int i=0;i<list_cand6.at(idx_cand).size();++i)
    {
      ids = list_cand6.at(idx_cand).at(i);
      bestTC.push_back(hits[ids]);
    }
  }
  else
  {
    for (unsigned int i=0;i<list_scores5.size();++i)
    {
      if (list_scores5.at(i)<score_min)
      {
	score_min = list_scores5.at(i);
	idx_cand  = i;
      }
    }
        
    for (unsigned int i=0;i<list_cand5.at(idx_cand).size();++i)
    {
      ids = list_cand5.at(idx_cand).at(i);
      bestTC.push_back(hits[ids]);
    }
  }
        
  //
  // Now loop again on the stubs in order to estimate the track parameters
  //
        
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
    for (int jj=0;jj<2;++jj) parZX[i][jj] = 0.;
    for (int jj=0;jj<2;++jj) invZX[i][jj] = 0.;
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


  for (int kk=0;kk<size;++kk) 
  {
    lay  = bestTC.at(kk)->getLayer();
    
    if (lay<=7)
    {
      (lay<=2)
	? nPSb++
	: nPSe++;
    }
    else
    {
      (lay<11)
	? n2Sb++
	: n2Se++;
    }
            
    if (lay>10) continue; // 2S disks stubs not used
    ++npt;
        
    x = bestTC.at(kk)->getX()*ci+bestTC.at(kk)->getY()*si;
    y = -bestTC.at(kk)->getX()*si+bestTC.at(kk)->getY()*ci;
    r = sqrt(x*x+y*y);
    
    if (r<rmin)
    {
      rmin = r;
      imin = kk;
      x2   = x;
      y2   = y;
    }
        
    if (r>rmax)
    {
      rmax = r;
      imax = kk;
      x0   = x;
      y0   = y;
    }
  }

  float rmax2 = 0;
    
  x1 = 0;
  y1 = 0;
  
  int nc=0;
    
  float xtemp1,ytemp1;
  float xtemp2,ytemp2;
  
  xtemp1=0.;
  ytemp1=0.;

  // Loop 2 over stubs in the TC
  // In order to determine the point with the second largest radius 

  for (int kk=0;kk<size;++kk) // Loop over stubs in the TC
  {
    if (kk==imax || kk==imin) continue; // Already used

    lay  = bestTC.at(kk)->getLayer();

    barrelmod=0;
    if (lay<=2 || (lay>=8 && lay<=10)) barrelmod=1;
    if (!barrelmod && (nPSb+n2Sb)>=3) continue; // Barrel modules have a priority
    if (lay>10 && (nPSb+nPSe)>=3) continue;     // Don't use 2S modules in the disks if possible
    
    x = bestTC.at(kk)->getX()*ci+bestTC.at(kk)->getY()*si;
    y = -bestTC.at(kk)->getX()*si+bestTC.at(kk)->getY()*ci;
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

  xtemp1 = (x1-x0)/((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
  ytemp1 = (y1-y0)/((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0));
  
  xtemp2 = (x2-x0)/((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
  ytemp2 = (y2-y0)/((x2-x0)*(x2-x0)+(y2-y0)*(y2-y0));
  
  x1=xtemp1;
  y1=ytemp1;
  x2=xtemp2;
  y2=ytemp2;
  

  // Now we got everything for the r/phi plane   

  double a = x0-0.5*(y1-y2)/(y2*x1-y1*x2);
  double b = y0-0.5*(x1-x2)/(y1*x2-y2*x1);
    
  int charge =-b/fabs(b);
   
  phi_est = atan2(charge*a,-charge*b)+sec_phi;
  pt_est  = 0.003*3.833*sqrt(a*a+b*b);

  // Then we do the RZ fit (LS)

  float wght;
  int cnt=0;
  for (int kk=0;kk<size;++kk) // Loop over stubs in the TC
  {
    lay  = bestTC.at(kk)->getLayer();

    if (lay>7) continue; // Don't use 2S modules
    if (lay>2 && nPSb>=2) continue; // Don't use PS modules of the disks if number of good point in the barrel is sufficient        

    ++cnt;
    x = bestTC.at(kk)->getX()*ci+bestTC.at(kk)->getY()*si;
    y = -bestTC.at(kk)->getX()*si+bestTC.at(kk)->getY()*ci;    
    z = bestTC.at(kk)->getZ();        
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
  
  tracks.push_back(fit_track);
}

void TCBuilder::fit(){

  vector<Hit*> activatedHits;

  //////// Get the list of unique stubs from all the patterns ///////////
  set<int> ids;
  int total=0;
  for(unsigned int i=0;i<patterns.size();i++){
    vector<Hit*> allHits = patterns[i]->getHits();
    total+=allHits.size();
    for(unsigned int j=0;j<allHits.size();j++){
      pair<set<int>::iterator,bool> result = ids.insert(allHits[j]->getID());
      if(result.second==true)
	activatedHits.push_back(allHits[j]);
    }
  }
  fit(activatedHits);
}

TrackFitter* TCBuilder::clone(){
  TCBuilder* fit = new TCBuilder(nb_layers);
  fit->setPhiRotation(sec_phi);
  fit->setSectorID(sector_id);
  return fit;
}

