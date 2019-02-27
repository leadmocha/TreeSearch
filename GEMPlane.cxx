//*-- Author :    Ole Hansen, Jefferson Lab   11-Jan-2010

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TreeSearch::GEMPlane                                                     //
//                                                                          //
// A 1-dimensional readout plane of a GEM chamber.                          //
// Use two of these objects for a 2-d readout plane.                        //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include "GEMPlane.h"
#include "GEMHit.h"
#include "GEMTracker.h"
#include "Projection.h"

#include "THaDetMap.h"
#include "TClonesArray.h"
#include "TCollection.h"
#include "TH1.h"
#include "TError.h"
#include "TString.h"
#include <TMinuit.h>
#include <TF1.h>

#ifndef MCDATA
#include "THaEvData.h"
#else
#include "SimDecoder.h"
//#include "TSolSimDecoder.h"
#endif

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <algorithm>

using namespace std;
using namespace Podd;

//#define PRINTCLUSTER
//#define PRINT_CLUSTER_DETAIL
#define TIME_SPLIT

#define ALL(c) (c).begin(), (c).end()

extern Double_t      fpulsex[20];
extern Double_t      fpulsey[20];

namespace TreeSearch {
  // Got this map from Danning. Maps the read channel number to corresponding
  // APV channel number. (The mapping is due to the multiplexing nature of
  // the readout).
  const int kMPDAPVMAP[128] = {1, 33, 65, 97, 9, 41, 73, 105, 17, 49, 81, 113,
    25, 57, 89, 121, 3, 35, 67, 99, 11, 43, 75, 107, 19, 51, 83, 115, 27, 59,
    91, 123, 5, 37, 69, 101, 13, 45, 77, 109, 21, 53, 85, 117, 29, 61, 93,
    125, 7, 39, 71, 103, 15, 47, 79, 111, 23, 55, 87, 119, 31, 63, 95, 127,
    0, 32, 64, 96, 8, 40, 72, 104, 16, 48, 80, 112, 24, 56, 88, 120, 2, 34,
    66, 98, 10, 42, 74, 106, 18, 50, 82, 114, 26, 58, 90, 122, 4, 36, 68, 100,
    12, 44, 76, 108, 20, 52, 84, 116, 28, 60, 92, 124, 6, 38, 70, 102, 14, 46,
    78, 110, 22, 54, 86, 118, 30, 62, 94, 126};


  // Sanity limit on number of channels (strips)
  static const Int_t kMaxNChan = 20000;

  //_____________________________________________________________________________
  GEMPlane::GEMPlane( const char* name, const char* description,
		      THaDetectorBase* parent )
    : Plane(name,description,parent),
      fMapType(kOneToOne), fMaxClusterSize(0), fMinAmpl(0), fSplitFrac(0),
      fMaxSamp(1), fAmplSigma(0), fADCraw(0), fADC(0), fHitTime(0), fADCcor(0), fMCCharge(0), 
      fGoodHit(0), fDnoise(0), fNrawStrips(0), fNhitStrips(0), fHitOcc(0),
      fOccupancy(0), fADCMap(0)
  {
    fStripsPerChannel = 1;
    fStripOrderedData = 0;
    // Constructor

    static const char* const here = "GEMPlane";

    // c1 = new TCanvas(name,name,800,600); 

    assert( dynamic_cast<GEMTracker*>(fTracker) );
    fTotalPriHit = 0, fCoveredPriHit = 0;

    try {
#ifdef MCDATA
      if( fTracker->TestBit(Tracker::kMCdata) ) // Monte Carlo data mode?
	fHits = new TClonesArray("TreeSearch::MCGEMHit", 200);
      else
#endif
	fHits = new TClonesArray("TreeSearch::GEMHit", 200);
    }
    catch( std::bad_alloc ) {
      Error( Here(here), "Out of memory allocating hit array for readout "
	     "plane %s. Call expert.", name );
      MakeZombie();
      return;
    }
  }

  //_____________________________________________________________________________
  GEMPlane::~GEMPlane()
  {
    // Destructor.

    // Histograms in fHist should be automatically deleted by ROOT when
    // the output file is closed

    if( fIsSetup )
      RemoveVariables();

    // fHits deleted in base class
    delete fGoodHit;
    delete fADCcor;
    delete fMCCharge;
    delete fHitTime;
    delete fADC;
    delete fADCraw;
  }

  //_____________________________________________________________________________
  Int_t GEMPlane::Begin( THaRunBase* run )
  {
    // Book diagnostic histos if not already done

    Plane::Begin( run );

#ifdef TESTCODE
    if( TestBit(kDoHistos) ) {
      //TODO: check if exist
      string hname(fPrefix);
      hname.append("adcmap");
      fADCMap = new TH1F( hname.c_str(), hname.c_str(), fNelem, 0, fNelem );
    }
#endif
    return 0;
  }

  //_____________________________________________________________________________
  void GEMPlane::Clear( Option_t* opt )
  {
    // Clear event-by-event data (hits)

    Plane::Clear(opt);

    if( !IsDummy() ) {
      assert( fADCraw and fADC and fHitTime and fADCcor and fGoodHit and fMCCharge);
      memset( fADCraw, 0, fNelem*sizeof(Float_t) );
      memset( fADC, 0, fNelem*sizeof(Float_t) );
      memset( fHitTime, 0, fNelem*sizeof(Float_t) );
      memset( fADCcor, 0, fNelem*sizeof(Float_t) );
      memset( fMCCharge, 0, fNelem*sizeof(Float_t) );
      memset( fGoodHit, 0, fNelem*sizeof(Byte_t) );
      fSigStrips.clear();
      fStripsSeen.assign( fNelem, false );
      fStrips.clear();
      fStrips.reserve(fNelem);
    }

    fNhitStrips = fNrawStrips = 0;
    fHitOcc = fOccupancy = fDnoise = 0.0;
  }

  //_____________________________________________________________________________
  Int_t GEMPlane::MapChannel( Int_t idx ) const
  {
    // Map hardware channel number to logical strip number based on mapping
    // prescription from database

    assert( idx >= 0 and idx < fNelem );

    Int_t ret = 0;
    switch( fMapType ) {
    case kOneToOne:
      ret = idx;
      break;
    case kReverse:
      ret = fNelem-idx-1;
      break;
    case kGassiplexAdapter1:
      assert( idx < 240 );
      if( idx == 0 )
	ret = 1;
      else if( idx == 239 )
	ret = 238;
      else if( idx % 2 ) // odd
	ret = idx + 2;
      else               // even
	ret = idx - 2;
      break;
    case kGassiplexAdapter2:
      assert( idx < 240 );
      if( idx == 1 )
	ret = 0;
      else if( idx == 238 )
	ret = 239;
      else if( idx % 2 ) // odd
	ret = idx - 2;
      else               // even
	ret = idx + 2;
      break;
    case kTable:
      // Use the mapping lookup table
      assert( fChanMap.size() == static_cast<Vint_t::size_type>(fNelem) );
      ret = fChanMap[idx];
      break;
    default:
      ret = -1;
    }
    assert( ret >= 0 and ret < fNelem );
    return ret;
  }

  //_____________________________________________________________________________
  Int_t GEMPlane::MapChannel( Int_t effChan, Int_t chan) const
  {
    assert(chan>=0 && chan<128);
    // Right now, the only map type that should call this is kMPDAPV/V2
    Int_t nvals = 0;  // Number of map entries per apv
    Int_t impd = 0;   // index of mpd_id
    Int_t iadc = 0;   // index of adc_id
    Int_t ipos = 0;   // index of pos
    Int_t iinv = 0;   // index of invert
    switch(fMapType) {
      case kMPDAPV:
        nvals = 8;
        impd = 2;
        iadc = 4;
        ipos = 6;
        iinv = 7;
        break;
      case kMPDAPV2:
        nvals = 6;
        impd = 0;
        iadc = 2;
        ipos = 4;
        iinv = 5;
        break;
      default:
        break;
    }
    assert(nvals>0); // This ensures only a valid MPD map calls this version
    assert(effChan>=0 && effChan < 65536); // Must fit into 16 bits
    Int_t mpd_id = (effChan>>8)&0xFF;
    Int_t adc_id = effChan&0xFF;
    Int_t n = fChanMap.size()/nvals;
    Int_t ret = -1;
    for(Int_t i = 0; i < n && ret==-1; i++) {
      if( mpd_id == fChanMap[impd + i*nvals] &&
          adc_id == fChanMap[iadc + i*nvals]) {
        ret = kMPDAPVMAP[chan];
        ret += (127-2*ret)*fChanMap[iinv + i*nvals];
        ret += 128*fChanMap[ipos + i*nvals];
      }
    }
    assert(ret>=0 && ret<= fNelem);
    return ret;
  }

  //for fitting pulse shape
Double_t Pdf(Double_t *xPtr, Double_t par[])
{
  Double_t x = *xPtr;
  Double_t A = par[0]; //normalization
  Double_t start = par[2]-par[1]; //Start Time
  Double_t shape = par[1]; //Shapping Time
  Double_t f=0;
  f = A * (x-start)/(shape)*TMath::Exp(-(x-start)/shape);
  return f;
}
Double_t CalcNLL(Int_t nTS, Double_t *fpulsex, Double_t *fpulsey, Double_t Sigma, int& npar, double par[])
  {
  //NLL calculation Function
  //nTS -> # of Time Bins AKA # of data points
  //X -> coordinates of x 
  //Y -> coordinates of y
  //Sigma -> error on y
  //f -> model
  Double_t nll = 0;
  for(int i=0;i<nTS;i++)
    {
      Double_t x = fpulsex[i];
      Double_t y = fpulsey[i];
      Double_t sigma = Sigma;
      Double_t mu = Pdf(&x, par);
      if(mu < 1e-10) mu = 1e-10;  // avoid log(0) problems
      nll -= (-0.5)*TMath::Log(2*TMath::Pi()*sigma) - 0.5*(y-mu)*(y-mu)/sigma/sigma;
    }
  return 2*nll;  //factor of -2 so minuit gets the errors right
}
void fcn(int& npar, double* deriv, double& f, double par[], int flag)
{
  //Function to be passed to TMinuit; This function is in fixed format
  //  for (int i=0; i<npar; i++)
  // {
  //   fparam->SetParameter(i,par[i]);
  // }
  f = CalcNLL(6, fpulsex, fpulsey, 25, npar, par);
}
  
  static Double_t FitPulse(Double_t* adc,Int_t nsample,
			   Double_t& shapingTime ,Double_t& peakTime ,Double_t& amp)
  {
    for(Int_t i=0;i<nsample;i++)
      {
	fpulsex[i]=i*25;
	fpulsey[i]=adc[i];	  
      }
   
    //For the Fitting
    Int_t npar = 3;
    TMinuit minuit(npar);
    minuit.SetPrintLevel(-1);
    minuit.SetFCN(fcn); //fcn cannot be a member function, so define it to be global; fcn calculates Chi^2 or NLL

    //define the model to be used
    Double_t xmin = -500;
    Double_t xmax = (nsample + 3)*25.0;
    TF1 fV("fV", Pdf, xmin, xmax, npar);

    //parameters for the model
  Double_t par[npar];
  Double_t stepSize[npar];
  Double_t minVal[npar];
  Double_t maxVal[npar];
  TString parName[npar];
  //inital value
  par[0] = 4096;//4200;       // guesses for starting the fit
  par[1] = 85;       // this MUST be done by some means to get things started
  par[2] = (3 + 0.5)*25.0 -60.0;//120;       //start time
  //step size
  stepSize[0] = TMath::Abs(par[0]*0.0001);   // usually 10% is OK for an initial step size, YMMV
  stepSize[1] = TMath::Abs(par[1]*0.001);   // step size MUST be positive!
  stepSize[2] = TMath::Abs(par[2]*0.001);
  //parameter range
  minVal[0] = 0;      // if min and max values = 0, parameter is unbounded.
  maxVal[0] = 0;//5000000;
  minVal[1] = 0;//30; 
  maxVal[1] = 100;//300;
  minVal[2] = 0;//10;
  maxVal[2] = 6*25;//300;
  //parameter name
  parName[0] = "A";
  parName[1] = "shape";
  parName[2] = "peak";

  for (int i=0;i<npar;i++)
    {
      minuit.DefineParameter(i, parName[i].Data(), par[i], stepSize[i], minVal[i], maxVal[i]);
    }
  minuit.Migrad(); //Minuit's best minimization do the fitting [4.4s/5000]

  //get fitting results
  Double_t outpar[npar];
  Double_t err[npar];
  for(int i=0;i<npar;i++)
    {
      minuit.GetParameter(i,outpar[i],err[i]);
    }

  fV.SetParameters(outpar);

  //Fill timing variables
  // StartTime = fV.GetParameter(2);
  shapingTime = fV.GetParameter(1);
  peakTime = fV.GetParameter(2);
  //Fitting_ReducedChisquare = fV.GetChisquare()/fV.GetNDF();
  amp = fV.GetMaximum(xmin, xmax);
  //amp = fV.GetMaximumX();

  return 0;
  }

  //_____________________________________________________________________________
  static StripData_t ChargeDep( const vector<Float_t>& amp )
  {
    // Deconvolute signal given by samples in 'amp', return approximate integral.
    // Currently analyzes exactly 3 samples.
    // From Kalyan Allada
    // NIM A326, 112 (1993)

    //FIXME: from database, proper value for Tp
    const Float_t delta_t = 25.0; // time interval between samples (ns)
    const Float_t Tp      = 50.0; // RC filter time constant (ns)

    assert( amp.size() >= 3 );
    Float_t adcraw = delta_t*(amp[0]+amp[1]+amp[2]);
    // Weight factors calculated based on the response of the silicon microstrip
    // detector:
    // v(t) = (delta_t/Tp)*exp(-delta_t/Tp)
    // Need to update this for GEM detector response(?):
    // v(t) = A*(1-exp(-(t-t0)/tau1))*exp(-(t-t0)/tau2)
    // where A is the amplitude, t0 the begin of the rise, tau1 the time
    // parameter for the rising edge and tau2 the for the falling edge.

    Float_t x = delta_t/Tp;

    Float_t w1 = TMath::Exp(x-1)/x;
    Float_t w2 = -2*TMath::Exp(-1)/x;
    Float_t w3 = TMath::Exp(-x-1)/x;

    // Deconvoluted signal samples, assuming measurements of zero before the
    // leading edge
    Float_t sig[3] = { amp[0]*w1,
		       amp[1]*w1+amp[0]*w2,
		       amp[2]*w1+amp[1]*w2+amp[0]*w3 };

    Float_t adc    = delta_t*(sig[0]+sig[1]+sig[2]);
    Float_t time   = 0;     // TODO
    Bool_t pass;
    // Calculate ratios for 3 samples and check for bad signals
    if( amp[2] > 0 ) {
      Float_t r1 = amp[0]/amp[2];
      Float_t r2 = amp[1]/amp[2];
      pass = (r1 < 1.0 and r2 < 1.0 and r1 < r2);
    } else
      pass = false;
    //printf("adcraw = %1.2f, sig[0, 1, 2] = %1.2f, %1.2f, %1.2f, adc = %1.2f \n", adcraw, sig[0], sig[1], sig[2], adc);
    return StripData_t(adcraw,adc,time,time,pass,amp);
  }

  //_____________________________________________________________________________
  bool GEMPlane::AnalyzeStrip( const vector<Float_t>& amp, StripData_t &stripdata )
  {
    //Deconvolution method is not proper for GEM signals since the "zero" time isn't fixed to a certain sample,
    //Especially with APV-25 chip, no syncronization between 40MHz clock and trigger, 
    //making the "starting time" may jitter up to an additional 25 ns.
    //Danning Di --Oct 2017

    //This routine takes the maximum charge in the samples
    Float_t maxAdc=0, maxAdc_t=0, adcSum=0;
    Int_t maxTimeSample=0;
    Bool_t pass = false;
    Double_t shapingtime=0, peaktime=0, adcmax_fit=0;
    Double_t sampleAdc[fMaxSamp];

    Int_t Nsample = amp.size();
    for(Int_t isamp=0;isamp<Nsample;isamp++)
      {
	sampleAdc[isamp] = amp[isamp];
	//if(isamp>0&&isamp<4){
	if(maxAdc<amp[isamp]) 
	  {
	    maxAdc=amp[isamp];
	    maxTimeSample=isamp;
	  }
	if(isamp>0&&isamp<4){
	  if(maxAdc_t<amp[isamp]) 
	    {
	      maxAdc_t=amp[isamp];
	      
	    }
	}
	//}
	adcSum+=amp[isamp];
	//	cout<<adcSum<<" "<<amp[isamp]<<endl;
      }
    Float_t adcAvg = adcSum/Nsample;
    //if(adcSum<=0){cout<<" AA "<<endl;getchar();}
    

    if(adcAvg>fpedestal_sigma*ftmp_pedestal_rms)// pedestal parameters to be put in data base
      {
	++fNhitStrips;
	if(maxTimeSample==0||maxTimeSample==(Nsample-1))
	  pass = true;
	else{
	  pass  = true;
	  FitPulse(sampleAdc,Nsample, shapingtime, peaktime, adcmax_fit);
	}

#ifdef TIME_SPLIT
	
	if(peaktime<=0 || peaktime>=150 || maxAdc_t<100){
	  pass = false;
	  return false;
	}

#endif

	//	stripdata =   StripData_t(maxAdc,adcSum,maxTimeSample,peaktime,pass,amp);
	stripdata.maxAdc = maxAdc_t;
	stripdata.adcSum = adcSum;
	stripdata.maxTimeSample = maxTimeSample;
	stripdata.peaktime = peaktime;
	stripdata.pass = pass;
	stripdata.vADC = amp;
	return true;
      }
    // if(pass){cout<<adcAvg<<" time: "<<peaktime<<endl;getchar();}
    // return StripData_t(maxAdc,adcSum,maxTimeSample,peaktime,pass,amp);
    return false;
  }



  //_____________________________________________________________________________
  void GEMPlane::AddStrip( Int_t istrip, Int_t module )
  {
    // Record a hit on the given strip number in internal arrays.
    // Utility function used by Decode.

    Float_t adc = mStrip[istrip].maxAdc;

    if( fGoodHit[istrip] and mStrip[istrip].pass ) {
      fSigStrips.push_back(istrip);
      fmStripModule[istrip] = module;
    }
#ifdef TESTCODE
    if( TestBit(kDoHistos) ) {
      fHitMap->Fill(istrip);
      fADCMap->Fill(istrip, adc);
    }
#endif
  }

//_____________________________________________________________________________
  void GEMPlane::DoClustering(){ // aim to be able to seperate cluster in high occupancy case by using strip timing information


}


  //_____________________________________________________________________________
  Int_t GEMPlane::GEMDecode( const THaEvData& evData )
  {
    // Extract this plane's hit data from the raw evData.
    //
    // This routine decodes the front-end readout data.
    // Finds clusters of active strips (=above software threshold) and
    // computes weighted average of position. Each such cluster makes one "Hit".
    
    int flagt=0,flagc=0;// to be removed
    
    const char* const here = "GEMPlane::Decode";
#ifdef PRINTCLUSTER
    cout<<"Decode Plane: "<<this->GetName()<<endl;
#endif

#ifdef MCDATA
    bool mc_data = fTracker->TestBit(Tracker::kMCdata);
    assert( !mc_data || dynamic_cast<const SimDecoder*>(&evData) != 0 );
#endif
    assert( fADCraw and fADC and fADCcor and fHitTime and fMCCharge);

#ifdef TESTCODE
    if( TestBit(kDoHistos) )
      assert( fHitMap != 0 and fADCMap != 0 );
#endif
    assert( fPed.empty() or
	    fPed.size() == static_cast<Vflt_t::size_type>(fNelem) );
    assert( fSigStrips.empty() );
    assert( fStripsSeen.size() == static_cast<Vbool_t::size_type>(fNelem) );

    UInt_t nHits = 0;

    bool do_pedestal_subtraction = !fPed.empty();
    bool do_noise_subtraction    = TestBit(kDoNoise);

#ifdef MCDATA
    const TSBSSimDecoder* simdata = 0;
    if( mc_data ) {
      assert((&evData) );
      assert( dynamic_cast<const TSBSSimDecoder*>(&evData) );
      simdata = static_cast<const TSBSSimDecoder*>(&evData);
    }
#endif

    
    // Decode data
    // Each DetMap->Module stands for a MPD
    for( Int_t imod = 0; imod < fDetMap->GetSize(); ++imod ) {
      THaDetMap::Module * d = fDetMap->GetModule(imod);
      // GEM module ID
      Int_t moduleID = d->plane; //0 1 2....
      // Typically, fcModeSize is 128 and each channel in the detmap
      // corresponds to one strip. However, sometimes, one channel corresponds
      // to in APV25, and hence fcModeSize is set to 1, and
      // fStripsPerChannel are set to 128
      assert((d->hi-d->lo+1)%fcModeSize==0);// checking total channels is 128*n
      Int_t Napvs = (d->hi-d->lo+1)/fcModeSize; // Number of groups of 128 strips
      if(Napvs==0) // Just in case less than 128 channels. Should be unusual,
        Napvs = 1; // but better safe than sorry!

      // First we read all the data from evData and organize them like so:
      // [apv][strip][adc_samples], where the actual strip number is recorded
      // in another array.
      std::vector<Vflt_t > samples_data[Napvs];
      Vint_t strips[Napvs];

      // Read the all channels in current DAQ-module
      Int_t nchan = evData.GetNumChan( d->crate, d->slot ); //

      for( Int_t ichan = 0; ichan < nchan; ++ichan ) {
        Int_t chan = evData.GetNextChan( d->crate, d->slot, ichan );
        Int_t chan_idx = d->first + ((d->reverse) ? d->hi - chan : chan - d->lo);
        //	cout<<chan<<"  "<<ichan<<endl;
        if( chan < d->lo or chan > d->hi ) continue; // not part of this detector
        Int_t apv = (int)(chan_idx/fcModeSize);
        assert(apv<Napvs);

        Int_t nsamp = evData.GetNumHits( d->crate, d->slot, chan );
        if( fStripsPerChannel == 1) {
          // In the simple case, each channel in the detmap corresponds to
          // just one strip.
          Int_t istrip = MapChannel(chan_idx);
          // Test for duplicate istrip, if found, warn and skip
          assert( istrip >= 0 and istrip < fNelem );
          if( fStripsSeen[istrip] ) {
            Warning( Here(here), "Duplicate strip number %d in plane %s, event %d. "
                "Ignorning it. Fix your detector crate slot map.",
                istrip, GetName(), evData.GetEvNum() );
            continue;
          }

          fStripsSeen[istrip] = true;
          fStrips.push_back(istrip);
          // For the APV25 analog pipeline, multiple "hits" on a decoder channel
          // correspond to time samples 25 ns apart
          Int_t nsamp = evData.GetNumHits( d->crate, d->slot, chan );
          assert( nsamp > 0 );
          ++fNrawStrips;
          nsamp = TMath::Min( nsamp, static_cast<Int_t>(fMaxSamp) );  
          strips[apv].push_back(istrip);

          // populate vector commonMode[][], all information from data is here
          Vflt_t samps;
          samps.resize(nsamp);
          for( Int_t isamp = 0; isamp < nsamp; ++isamp ) {
            samps[isamp] = static_cast<Float_t> (
                evData.GetData(d->crate, d->slot, chan, isamp) );
          }
          samples_data[apv].push_back(samps);

#ifdef MCDATA
          // If doing MC data, save the truth information for each strip
          if( mc_data ) {
            Double_t fmccharge=0;
            fMCHitList.push_back(istrip);
            fMCHitInfo[istrip] = simdata->GetSBSMCHitInfo(d->crate,d->slot,chan);
            //	fMCCharge[istrip]  = fmccharge;
            //cout<<fmccharge<<endl;
            //fHitTime[istrip] = fMCHitInfo[istrip].fMCTime;
          }
#endif
        } else if (fStripOrderedData) { // Each channel is more than one strip, and is strip ordered
          nsamp /= fStripsPerChannel;
          if(nsamp > fMaxSamp)
            nsamp = fMaxSamp;
          Int_t isamp = -1;
          Int_t strip = -1;
          Int_t ihit = 0;
          Vflt_t samps[fStripsPerChannel];
          while(++isamp < nsamp) {
            for(Int_t istrip = 0; istrip < Int_t(fStripsPerChannel); istrip++,ihit++) {
              strip = MapChannel( chan, istrip );
              assert( strip >= 0 and strip < fNelem);
              if(isamp==0 && fStripsSeen[strip]) {
                Warning( Here(here), "Duplicate strip number %d in plane %s, event %d. "
                    "Ignorning it. Fix your detector crate slot map.",
                    istrip, GetName(), evData.GetEvNum() );
              } else {
                // For the APV25 analog pipeline, with no SSP, multiple "hits"
                // on a decoder channel correspond to time samples 25 ns apart
                // organized as [strip1_ts1, strip2_ts1, .... stripn_ts1][strip1_ts2, ..., [stripn_ts2]....
                // This is a special case where chan is (mpd_id<<8)|adc_id
                if(isamp==0) {
                  samps[istrip].resize(nsamp);
                  fStripsSeen[strip] = true;
                }
                samps[istrip][isamp] = static_cast<Float_t> (
                      evData.GetData(d->crate, d->slot, chan, ihit) );
                if(do_pedestal_subtraction)
                  samps[istrip][isamp]-=fPed[istrip];
              }
            }
          }
          for(Int_t istrip = 0; istrip < fStripsPerChannel; istrip++) {
            Vflt_t lsamps;
            lsamps.reserve(samps[istrip].size());
            for(size_t k = 0; k < samps[istrip].size(); k++) {
              lsamps.push_back(samps[istrip][k]);
            }
            samples_data[apv].push_back(lsamps);
            strip = MapChannel( chan, istrip );
            fStrips.push_back(strip);
            strips[apv].push_back(strip);
            ++fNrawStrips;
          }
        } else { // When each channel corresponds to more than one strip and is sample ordered

          // The raw adc value is the apv_ch_number (from 0-127)
          Int_t isamp = 0;
          Int_t istrip = -1;
          while(isamp < nsamp) {
            istrip = MapChannel( chan, evData.GetRawData(
                  d->crate,d->slot,chan,isamp) );

            // Test for duplicate istrip, if found, warn and skip
            assert( istrip >= 0 and istrip < fNelem );
            if( fStripsSeen[istrip] ) {
              Warning( Here(here), "Duplicate strip number %d in plane %s, event %d. "
                  "Ignorning it. Fix your detector crate slot map.",
                  istrip, GetName(), evData.GetEvNum() );
              while(++isamp<nsamp && istrip ==
                  MapChannel( chan, evData.GetRawData(
                      d->crate,d->slot,chan,isamp) ) ) {
                // Do nothing, just increment isamp
              }
            } else {
              // For the APV25 analog pipeline, multiple "hits" on a decoder
              // channel correspond to time samples 25 ns apart *for one or
              // more strips.* This is a special case where chan is
              // (mpd_id<<8)|adc_id
              Vflt_t samps;
              samps.reserve(fMaxSamp+1);
              Float_t sampv = 0; // sample value
              sampv = evData.GetData(d->crate, d->slot, chan, isamp);//
              if(do_pedestal_subtraction)
                sampv -= fPed[istrip];
              samps.push_back(sampv);
              while(++isamp<nsamp && istrip ==
                  MapChannel( chan, evData.GetRawData(
                      d->crate,d->slot,chan,isamp) ) ) {
                if(samps.size() < fMaxSamp) {
                  sampv = evData.GetData(d->crate, d->slot, chan, isamp);
                  if(do_pedestal_subtraction)
                    sampv-=fPed[istrip];
                  samps.push_back(sampv);
                }
                        
              }
              samples_data[apv].push_back(samps);
              fStripsSeen[istrip] = true;
              fStrips.push_back(istrip);
              strips[apv].push_back(istrip);
              ++fNrawStrips;
            } // end loop while isamp <nsamp
          } // end if fStripsPerChannel > 1
        } // end while loop isamp<nsamp
      } // end loop ichan

      // What follows now is filling the Common mode storage that was here
      // previously, but now takes data from the sample_data and strips
      // vectors instead of evData.

      // (Juan Carlos Cornejo Oct 24, 2018): It doesn't seem like the common
      // mode stuff  is actually used by anything, but I'm leaving it
      // intact in case someone intended to modify and make use of it in the
      // future.
      vector<Float_t> commonMode[Napvs][fMaxSamp]; // storing value for calculating common mode
      Float_t cMode[Napvs][fMaxSamp]; //common mode value of each group of 128 channels and each sample

      // Calculating common mode from vector commonMode[][] and store in Int cMode[][]
      // 1 loop check adc range, could add more loop to get more precise commonMode, this will be O(n), better than sorting O(nlogn)/O(n^2), time matters for online processing
      for(int iapv=0; iapv<Napvs;iapv++)
	{
          // October 24, 2018 (jc2): Fill commonMode here so that the following
          // loop's logic doesn't change. (commonMode is an array with the
          // following structure [apv][sample_number][adc_sample]
          for(size_t istrip = 0; istrip < samples_data[iapv].size(); istrip++) {
            for(UInt_t isamp  = 0; isamp<samples_data[iapv][istrip].size();
                isamp++) {
              assert(isamp<6);
              commonMode[iapv][isamp].push_back(
                  samples_data[iapv][istrip][isamp]);
            }
          }
          // Back to the old logic of using commonMode...
	  // cout<<d->crate<<" "<<d->slot<<" apv: "<<iapv<<endl;
	  for(UInt_t isamp=0; isamp<fMaxSamp; isamp++)
	    {
	      int Commsize = commonMode[iapv][isamp].size();
        Vflt_t common_sort;
        common_sort.insert(common_sort.begin(),&(commonMode[iapv][isamp][0]),&(commonMode[iapv][isamp][Commsize-1]));
        sort(common_sort.begin(),common_sort.end());
     	      int NcommMode = 0;
	      Float_t tempcMode = 0;
        int common_mode_first = 0;
        int common_mode_last = Commsize;
        if(Commsize==128) {
          common_mode_first = 28;
          common_mode_last = 100;
        }
	      for(int is=common_mode_first;is<Commsize&&is<common_mode_last;is++)
		{
		  Float_t tempADC = common_sort[is];
		  if(Commsize!=128&&(tempADC>(ftmp_comm_range)||tempADC<-ftmp_comm_range))continue; // constant 200 to be put in database, this number should come from a clean, no signal run, and be determined by evaluating the common mode distribution.
		  tempcMode+=tempADC;
		  NcommMode++;
		}
	      cMode[iapv][isamp]=tempcMode/NcommMode;
	      //if(NcommMode<120){cout<<cMode[iapv][isamp]<<" # "<<NcommMode<<"        ";getchar();}
	    }
	  //	  cout<<endl;
	}



      // subtracting cMode[][]
      // do zero suppression, fit timing and store hits in fSigStrip
      // (Now uses the samples_data vector already filled above)
      for(int iapv=0; iapv<Napvs;iapv++) {
        for(size_t ist = 0; ist < samples_data[iapv].size(); ist++) {
          StripData_t stripdata;

          Vflt_t samples;
          UInt_t nsamps = (Int_t)(samples_data[iapv][ist].size());
          if( nsamps > 1 )
            samples.reserve(nsamps);

          Int_t istrip = strips[iapv][ist];
	
	for(UInt_t isamp=0; isamp<nsamps; isamp++)
	  {
	    Float_t fsamp = samples_data[iapv][ist][isamp];
            // 20181024 (jc2): Why is this commented out? Doesn't this
            // defeat the whole purpose of computing cMode above?
            if(do_noise_subtraction)
               fsamp-=cMode[iapv][isamp];
	    samples.push_back( fsamp );
	  }

	// Do zero suppression and Analyze the pulse shape
	if(!AnalyzeStrip(samples, stripdata))
	  continue;// Do nothing for strips not passing zero suppression

	// Save strip information for cluster finding, this struct have all information, in principle we don't need the following fADCraw[], fADC[].......
	mStrip[istrip] = stripdata;
      
	// Save results for cluster finding later
	// not necessary anymore, still here because some histograms needs it.
	fADCraw[istrip]  = stripdata.maxAdc;
	fADC[istrip]     = stripdata.adcSum;
	fHitTime[istrip] = stripdata.peaktime;
	fADCcor[istrip] = stripdata.adcSum;
	fGoodHit[istrip] = not TestBit(kCheckPulseShape) or stripdata.pass;

	AddStrip( istrip , imod);
        } // end loop on istrip
      } // end of loop on iapv
      // cout<<GetName()<<" "<<moduleID<<endl;
    }    // end of loop modules


    fHitOcc    = static_cast<Double_t>(fNhitStrips) / fNelem;
    fOccupancy = static_cast<Double_t>(GetNsigStrips()) / fNelem;

    //
    // DoClustering();
    //

    //Draw strips hit from mStrip
    //  cout<<" ### "<<endl;
    TH2F *histo = new TH2F("h2","h2",6000,0,6000,6,0,6);
    histo->Fill(1,1,0.1);
    for(std::map<Int_t,StripData_t>::iterator it=mStrip.begin();it!=mStrip.end();it++){
      int samplesize = (*it).second.vADC.size();
      for(int its=0;its<samplesize;its++){
	histo->Fill((*it).first, its, (*it).second.vADC[its]);
      }
    }
    //   c1->cd();
       //   histo->Draw("LEGO1");
    //histo->Draw("colz");
    //   c1->Draw();
    //  c1->Modified();
    //   c1->Update();
    //  cout<<histo->GetEntries();
    //   getchar();

    // c1->Delete();
     histo->Delete();
    //





    // Find and analyze clusters. Clusters of active strips are considered
    // a "Hit".
    //
    // The cluster analysis is a critical part of the GEM analysis. Various
    // things can and probably need to be done right here already: splitting
    // oversized clusters, detecting noise hits/bogus clusters, detecting and
    // fitting overlapping clusters etc.
    //
    // This analysis may even need to be re-done after preliminary tracking to
    // see if the clustering can be improved using candidate tracks.
    // Additionally, correlated amplitude information from a second readout
    // direction in the same readout plane could be used here. These advanced
    // procedures would require significant redesign of the code:
    // all raw strip info will have to be saved and prcessed at a later point,
    // similar to the finding of hit pairs in like-oriented planes of the MWDC.
    //
    // For the moment, we implement a very simple algorithm: any cluster of
    // strips larger than what a single cluster should be is assumed to be two or
    // more overlapping hits, and the cluster will be split as follows: anything
    // that looks like a local peak followed by a valley will be considered an
    // actual cluster. The parameter frac = fSplitFrac (0.0 ... 1.0) can
    // be used for some crude tuning. frac > 0.0 means that a peak is
    // only a peak if the amplitude drops below (1-frac), so
    // frac = 0.1 means: trigger on a drop below 90% etc. Likewise for the
    // following valley: the bottom is found if the amplitude rises again
    // by (1+frac), so frac = 0.1 means: trigger on a rise above 110% etc.

    // The active strip numbers must be sorted for the clustering algorithm
    sort( ALL(fSigStrips) );

    Double_t frac_down = 1.0 - fSplitFrac, frac_up = 1.0 + fSplitFrac;
#ifndef NDEBUG
    GEMHit* prevHit = 0;
#endif
    typedef Vint_t::iterator viter_t;
    Vint_t splits;  // Strips with ampl split between 2 clusters
    viter_t next = fSigStrips.begin();
    //strip no of last type 3
    Int_t t3_nb = -1;
    //type3Flag;
    //cout<<this->GetName()<<" nstrips fired:  "<<fNhitStrips<<endl;
    while( next != fSigStrips.end() ) {
      
      viter_t start = next, cur = next;
      ++next;
      Int_t moduleID = fmStripModule[*start];
      assert( next == fSigStrips.end() or *next > *cur );
      while( next != fSigStrips.end() and (*next - *cur) == 1  and fmStripModule[*next]==moduleID) {
	++cur;
	++next;
      }
      for( viter_t sstart=start; sstart != next; ++sstart ) {
	Int_t istrip = *sstart;
	Double_t pos = GetStart() + istrip * GetPitch();
	Double_t adc = fADCraw[istrip];
	//cout<<"XXSTRIP: "<<istrip<<" cc: "<<fMCHitInfo[istrip].fMCCharge<<" : "<<adc<<" : time: "<<fMCHitInfo[istrip].fMCTime<<endl;
      }
      // Now the cluster candidate is between start and cur
      assert( *cur >= *start );
      // The "type" parameter indicates the result of the cluster analysis:
      // 0: clean (i.e. smaller than fMaxClusterSize, no further analysis)
      // 1: large, maximum at right edge, not split
      // 2: large, no clear minimum on the right side found, not split
      // 3: split, well-defined peak found (may still be larger than maxsize)
      Int_t  type = 0;
      UInt_t size = *cur - *start + 1;
      if( size > fMaxClusterSize ) {
	Double_t maxadc = 0.0, minadc = kBig;
	viter_t it = start, maxpos = start, minpos = start;
	enum EStep { kFindMax = 1, kFindMin, kDone };
	EStep step = kFindMax;
	while( step != kDone and it != next ) {
	  Double_t adc = fADCcor[*it];
	  switch( step ) {
          case kFindMax:
            // Looking for maximum
            if( adc > maxadc ) {
              maxpos = it;
              maxadc = adc;
            } else if( adc < maxadc * frac_down ) {
              assert( maxadc > 0.0 );
              step = kFindMin;
              continue;
            }
            break;
          case kFindMin:
            // Looking for minimum
            if( adc < minadc ) {
              minpos = it;
              minadc = adc;
            } else if( adc > minadc * frac_up ) {
              assert( minadc < kBig );
              step = kDone;
            }
            break;
          case kDone:
            assert( false );  // should never get here
            break;
	  }
	  ++it;
	}
	if( step == kDone ) {
	  // Found maximum followed by minimum
	  assert( minpos != start );
	  assert( minpos != cur );
	  assert( *minpos > *maxpos );
	  // Split the cluster at the position of the minimum, assuming that
	  // the strip with the minimum amplitude is shared between both clusters
	  cur  = minpos;
	  next = minpos;
	  t3_nb = *minpos;
	  // In order not to double-count amplitude, we split the signal height
	  // of that strip evenly between the two clusters. This is a very
	  // crude way of doing what we really should be doing: "fitting" a peak
	  // shape and using the area and centroid of the curve
	  fADCcor[*minpos] /= 2.0;
	  splits.push_back(*minpos);
	}
	type = step;
	size = *cur - *start + 1;
	assert( *cur >= *start );
      }
      assert( size > 0 );
      // Compute weighted position average. Again, a crude (but fast) substitute
      // for fitting the centroid of the peak.
      Double_t xsum[fMaxSamp], x_sum=0, adcsum = 0.0, adcmax_fit=0.0, shapingtime=0, peaktime=0, prim_ratio=0, p_over_total=0, b_over_total=0, p_over_total_r=0, b_over_total_r=0, p_over_total_o=1, ncross_talk=0, amp_fit=0, pos_sigma=0.0, time_sigma=0.0, adc_sigma=0.0;
      bool prim_registered = false;
      Double_t pos_prim_tmp, time_prim_tmp, adc_prim_tmp;
      map<int,int> bg_registered;
      Int_t stripl = *start, striph = *(cur-1);
      viter_t strip_left = start; 
      viter_t strip_right = cur;
 
      while( *strip_left <= *strip_right ) {
	if(*strip_left == *strip_right || *strip_left == *strip_right-1){
	  shapingtime = 56.0; 
	  peaktime = (fHitTime[ *strip_left] + fHitTime[ *strip_right ])/2;
	  //  cout<<"apeakTime: "<<peaktime<<endl;
	  //	  if(peaktime>90 ||peaktime<30)
	  // getchar();
	  adcmax_fit = fADCraw[ *strip_left ] > fADCraw[ *strip_right ] ? fADCraw[ *strip_left ] : fADCraw[ *strip_right ];
	  break;
	}else{
	  strip_left++;
	  strip_right--;
	}
      }

      //cout<<stripl<<" "<<striph<<"  = "<<size<<endl;getchar();
      Double_t sampleAdc[fMaxSamp];
      for(Int_t i_sampleAdc=0;i_sampleAdc<fMaxSamp;i_sampleAdc++)
	{
	  xsum[i_sampleAdc]=0;
	  sampleAdc[i_sampleAdc]=0;
	}
      Int_t nsample=fMaxSamp;

#ifdef MCDATA
      Double_t mcpos = 0.0, mccharge = 0.0, mctime = kBig;
      Int_t mctrack = 0, num_bg = 0;
#endif
#ifdef PRINTCLUSTER
      // cout<<"new Cluster"<<endl;
#endif

      if(*start == t3_nb)
	type = 3;
      

      for( ; start != next; ++start ) {
	
	Int_t istrip = *start;
	Double_t pos = GetStart() + GetPitch()*(istrip -fDetMap->GetModule(fmStripModule[istrip])->first)+GetModuleOffsets(fmStripModule[istrip]);
	//	cout<<"##@@"<<pos<<endl;
	//	cout<<GetType()<<" "<<fDetMap->GetModule(0)->first<<"= = = "<<(GetPitch()*(fDetMap->GetModule(fmStripModule[istrip])->first)-GetModuleOffsets(fmStripModule[istrip]))<<endl;//getchar();
	//	static_cast<const SimDecoder*>(&evData)


	Double_t adc = fADCraw[istrip];
	//cout<<"STRIP: "<<istrip<<" cc: "<<fMCHitInfo[istrip].fMCCharge<<" : "<<adc<<" : time: "<<fMCHitInfo[istrip].fMCTime<<endl;
	x_sum   += pos * adc;
	adcsum += adc;
	nsample = mStrip[istrip].vADC.size();
#ifdef PRINT_CLUSTER_DETAIL
	cout<<setw(5)<<istrip<<" maxadc: "<<setw(6)<<(Int_t)mStrip[istrip].maxAdc<<" pktime: "<<(Int_t)mStrip[istrip].peaktime<<" sampleADC: ";
#endif
	for(Int_t i_ts=0;i_ts<nsample;i_ts++)
	  {
#ifdef PRINT_CLUSTER_DETAIL
	    cout<<setw(6)<<(Int_t)mStrip[istrip].vADC[i_ts]<<"  ";
#endif
	    sampleAdc[i_ts]+=mStrip[istrip].vADC[i_ts];
	    //cout<<" xsumB "<<xsum[i_ts]<<" ";
	    xsum[i_ts]     +=mStrip[istrip].vADC[i_ts]*pos;
	    //   cout<<" xsum "<<xsum[i_ts]<<" ";
	  }
#ifdef PRINT_CLUSTER_DETAIL
	cout<<endl;
#endif
#ifdef PRINT_CLUSTER_DETAIL
	cout<<setw(6)<<(Int_t)fADC[istrip]<<"  ";
	cout<<" maxTs: "<<mStrip[istrip].maxTimeSample;
#ifdef MCDATA
	cout<<" SigType: "<<fMCHitInfo[istrip].fSigType
	    <<" = ";
	if(fMCHitInfo[istrip].fSigType&0x1){	  cout<<"P "; getchar();	}
	if(fMCHitInfo[istrip].fSigType&0x2)cout<<"S ";
	if(fMCHitInfo[istrip].fSigType&0x4)cout<<"C ";
	cout<<endl;

	for(int i_c=0;i_c<fMCHitInfo[istrip].vClusterID.size();i_c++)
	  {
	    cout<<"                                  "<<setw(6)<<fMCHitInfo[istrip].vClusterID[i_c]<<" #";
	    for(int i_ts=0;i_ts<nsample;i_ts++)
	      cout<<setw(6)<<(Int_t)(100*fMCHitInfo[istrip].vClusterADC[i_ts][i_c]/mStrip[istrip].vADC[i_ts])<<"% ";
	    //cout<<endl;
	    cout<<" peaktime: "<<setw(5)<<(Int_t)fMCHitInfo[istrip].vClusterPeakTime[i_c];
	    cout<<" weight: "<<setw(5)<<std::fixed<<std::setprecision(2)<<fMCHitInfo[istrip].vClusterStripWeight[i_c]<<endl;
	  }
#endif // MCDATA
	cout<<endl;
#endif

#ifdef MCDATA
	if(fMCHitInfo[istrip].fSigType&0x4 ) 
	  ncross_talk++;
	//need to add 1). (pos_prim-pos_bg)^2,
	//            2). (time_prim-time_bg)^2
	//            3). (adc_prim-adc_bg)^2
	// These need to be from big enough bg hit,need threshhold, -----this is one way to do it, we can also include all, which can be used to see for instance
	// b_over_total-vs-pos_sigma  gets better with space p-v-p finding split, not very accurate tho, in case 2bg+1prim and 1 of the bg has same pos as prim. the other off 
	if(fMCHitInfo[istrip].fSigType&0x1 ){// do following only for strip with primary hit
	  for(int ts=0;ts<nsample;ts++){
	    p_over_total_o += mStrip[istrip].vADC[ts];
	  }
	  Double_t pos_prim_tmp, time_prim_tmp, adc_prim_tmp;
	  for(int i_c=0;i_c<fMCHitInfo[istrip].vClusterID.size();i_c++)//i_c<1,i_c=0 means primary hit....the first digitized cluster...bad
	    {
	      //
	      
	      if(fMCHitInfo[istrip].vClusterType[i_c]==0){// if primary(from signal file) hit
		prim_ratio+=fMCHitInfo[istrip].vClusterStripWeight[i_c];  
		for(int i_ts=0;i_ts<nsample;i_ts++){
		  p_over_total_r += fMCHitInfo[istrip].vClusterADC[i_ts][i_c];
		}
		//following line assumed prim hit doesn't overlap()
		if(!prim_registered){//if prim havn't been registered
		  //string st_tmp = this->fName;
		  if(this->fName[2]=='x')
		    {  pos_prim_tmp = fMCHitInfo[istrip].vClusterPos[i_c].X();
		  
		    }
		  else if(this->fName[2]=='y')
		    pos_prim_tmp = fMCHitInfo[istrip].vClusterPos[i_c].Y();
		  else
		    cout<<"ERROR in plane name!"<<endl;
		  time_prim_tmp = fMCHitInfo[istrip].vClusterPeakTime[i_c];
		  adc_prim_tmp = fMCHitInfo[istrip].vClusterCharge[i_c];
		  prim_registered = true;
		}
	      }else{//if secondary(from background file) hit
		for(int i_ts=0;i_ts<nsample;i_ts++){
		  b_over_total_r += fMCHitInfo[istrip].vClusterADC[i_ts][i_c];
		}
	      }
	    }
	
	  for(int i_c=0;i_c<fMCHitInfo[istrip].vClusterID.size();i_c++)//i_c<1,i_c=0 means primary hit....the first digitized cluster...bad
	    {
	      //  cout<<"type: "<<fMCHitInfo[istrip].vClusterType[i_c]<<endl;getchar();
	      if(fMCHitInfo[istrip].vClusterType[i_c]>0){//if bg
		//	cout<<this->fName<<endl;
		//	cout<<fMCHitInfo[istrip].vClusterPos[i_c].X()<<" "
		//	    <<fMCHitInfo[istrip].vClusterPos[i_c].Y()<<" "
		  //	    <<fMCHitInfo[istrip].vClusterPos[i_c].Z()<<" "
		//	cout<<this->fName[2]<<endl;
		//getchar();

		int cluster_id_tmp = fMCHitInfo[istrip].vClusterID[i_c];
		if(bg_registered.find(cluster_id_tmp) == bg_registered.end()){// if this bg not registered
		  Double_t pos_tmp;
		  if(this->fName[2]=='x')
		    { 
		      pos_tmp = fMCHitInfo[istrip].vClusterPos[i_c].X();
		      //cout<<fMCHitInfo[istrip].vClusterPos[i_c].X()<<endl;
		    }
		  else if(this->fName[2]=='y')
		    pos_tmp = fMCHitInfo[istrip].vClusterPos[i_c].Y();
		  else
		    cout<<"ERROR in plane name!"<<endl;
		
		  pos_sigma += pow((pos_tmp-pos_prim_tmp), 2);
		  time_sigma += pow((fMCHitInfo[istrip].vClusterPeakTime[i_c]-time_prim_tmp), 2);
		  adc_sigma += pow((fMCHitInfo[istrip].vClusterCharge[i_c]-adc_prim_tmp), 2);
		  bg_registered[cluster_id_tmp] = 0;
		}
	      }
	    }
	  // cout<< bg_registered.size()<<endl; getchar();
	  //TODO tomorrow,friday,    1) add these "sigma" into mchit, then write it to txt
	  //                         2) see the difference the space split made in these sigma-vs-b_over_total
	  //                         3) this tricky! and important to know how to split in time! check all the strip time distribution in a reconstructed cluster, clean hit, and clean hit with background
	  //                         4) see the difference the time split made! The time split has two level? 1. for some overlaped cluster, side strips will not be overlaped,                                                                                                                2. some totally overlaped in space, have to try treat in time, only 6 time sample....hard...:( de-convolution? no.t likely.. 

	}

	if(fMCHitInfo[istrip].fSigType&0x1){
	  
	  if(!flagt){fTotalPriHit++;flagt=1;}
	  //find out how many "valuable" strips covered by noticable background
	  if(fMCHitInfo[istrip].vClusterStripWeight[0]>0.1){
	    for(int i_ts=0;i_ts<nsample;i_ts++)
	      if(fMCHitInfo[istrip].vClusterADC[i_ts][0]/mStrip[istrip].vADC[i_ts]<0.8){
		if(!flagc){fCoveredPriHit++;flagc=1;}
	      }
	  }
	}


	// If doing MC data, analyze the strip truth information
	if( mc_data ) {
	  TSBSMCHitInfo& mc = fMCHitInfo[istrip];
	  // This may be smaller than the actual total number of background hits
	  // contributing to the entire cluster, but counting them would involve
	  // lists of secondary particle numbers ... overkill for now
	  num_bg = TMath::Max( num_bg, mc.fContam );
	  // All primary particle hits in the cluster are from the same track
	  assert( mctrack == 0 || mc.fMCTrack == 0 || mctrack == mc.fMCTrack );
	  //cout<<"OOSTRIP: "<<istrip<<" cc: "<<mc.fMCCharge<<" : "<<adc<<" : time: "<<mc.fMCTime<<endl;	
	  mccharge += mc.fMCCharge;
	  if( mctrack == 0 ) {
	    if( mc.fMCTrack > 0 ) {
	      // If the cluster contains a signal hit, save its info and be done
	      //printf("istrip = %d, pos = %1.3f, ADC = %1.3f\n", istrip, pos, adc);
	      //printf("mc.fMCTtime = %1.3f\n", mc.fMCTime);
	      mctrack = mc.fMCTrack;
	      mcpos   = mc.fMCPos;
	    
	    
	      mctime  = mc.fMCTime;
	    }
	    else {
	      // If background hits only, compute position average
	      mcpos  += mc.fMCPos;
	      mctime  = TMath::Min( mctime, mc.fMCTime );
	    }
	  }
	}
#endif // MCDATA
      }
      assert( adcsum > 0.0 );
      Double_t sPos[fMaxSamp],pos=0;
      //  cout<<"###########################\nsamplePos: ";
      for(Int_t i_ts=0;i_ts<nsample;i_ts++)
	{	
	  sPos[i_ts] = xsum[i_ts]/sampleAdc[i_ts];
	  //	  	  cout<<100*sPos[i_ts]<<" ";
       
	}
      pos=x_sum/adcsum;
      p_over_total = p_over_total_r/p_over_total_o;
      b_over_total = b_over_total_r/p_over_total_o;
      int number_of_bg = bg_registered.size();
      // cout<<"dd: "<<pos_sigma<<" "<<time_sigma<<" "<<adc_sigma<<endl;
      if(number_of_bg!=0){
	pos_sigma = sqrt(pos_sigma/number_of_bg);
	time_sigma = sqrt(time_sigma/number_of_bg);
	adc_sigma = sqrt(adc_sigma/number_of_bg);
      }
      //	cout<<"pos: "<<pos<<endl;
      //   cout<<"\n##########################"<<endl;getchar();
      
      //FitPulse(sampleAdc,nsample, shapingtime, peaktime, adcmax_fit);
      //cout<<"peakTime: "<<peaktime<<endl;
      // FIXME: What was this line supposed to do? It's missing brackets
      // and breaks several lines below. I commented out this line only,
      // not the getchar() that was already commented out.
      //if(peaktime>90 ||peaktime<30)
	//getchar();
      

      //  cout<<"shapingtime: "<<shapingtime<<" peaktime: "<<peaktime<<" amp: "<<amp_fit<<endl;
      //    cout<<"mctime-peaktime: "<<mctime-peaktime<<endl;
      //  cout<<"adcmax/ampfit: "<<adcsum<<" : "<<amp_fit<<endl;
#ifdef MCDATA
      if( mc_data && mctrack == 0 ) {
	mcpos /= static_cast<Double_t>(size);
      }
#endif
      // The resolution (sigma) of the position measurement depends on the
      // cluster size. In particular, if the cluster consists of only a single
      // hit, the resolution is much reduced
      Double_t resolution = fResolution;
      if( size == 1 ) {
	resolution = TMath::Max( 0.25*GetPitch(), fResolution );
	// The factor of 1/2*pitch is just a guess. Since with real GEMs
	// there _should_ always be more than one strip per cluster, we must
	// assume that the other strip(s) did not fire due to inefficiency.
	// As a result, the error is bigger than it would be if only ever one
	// strip fired per hit.
	//       resolution = TMath::Max( 0.5*GetPitch(), 2.0*fResolution );
	//     } else if( size == 2 ) {
	//       // Again, this is a guess, to be quantified with Monte Carlo
	//       resolution = 1.2*fResolution;
      }

      // Make a new hit

      if(peaktime<30 || peaktime>90){
	continue;
      }
#ifdef PRINT_CLUSTER_DETAIL
      if(p_over_total!=0 && b_over_total>0.3){
	cout<<"Writing a new cluster, prim ratio: "<<prim_ratio<<"  p_over_total: "<<p_over_total<<" bg: "<<b_over_total<<endl;
	getchar();
      }
      cout<<"####################bg: "<<b_over_total<<endl;
#endif

#ifndef NDEBUG
      GEMHit* theHit = 0;
#endif
#ifdef MCDATA
      if( !mc_data ) {
#endif
#ifndef NDEBUG
	theHit =
#endif
	  new( (*fHits)[nHits++] ) GEMHit( moduleID,
					   pos,
					   adcsum,
					   adcmax_fit,
					   peaktime,
					   size,
					   type,
					   resolution,
					   this
					   );

#ifdef MCDATA
      } else {
	// Monte Carlo data
      
#ifndef NDEBUG
	theHit =
#endif
	  new( (*fHits)[nHits++] ) MCGEMHit( moduleID,
					     pos,
					     adcsum,
					     adcmax_fit,
					     prim_ratio,
					     p_over_total,
					     b_over_total,
					     pos_sigma,
					     time_sigma,
					     adc_sigma,
					     bg_registered.size(),
					     peaktime,
					     size,
					     stripl,
					     striph,
					     type,
					     resolution,
					     this,
					     mctrack,
					     mcpos,
					     mccharge,
					     mctime,
					     num_bg
					     );

	//	cout<<pos<<" ### "<<mcpos<<" type: "<<type<<"charge: "<<adcmax_fit<<" "<<GetType()<<" "<<moduleID<<endl;
	//cout<<"cluster charge: "<<mccharge<<" adcSum:"<<adcsum<<endl;
	//printf("hit pos = %1.3f, time = %1.3f, size = %d\n", pos, mctime, size);
#define PRINTCLUSTER_TO_FILE
#ifdef PRINTCLUSTER_TO_FILE
	if(prim_ratio!=0){
	  std::ofstream outfile;
	  outfile.open("Cluster_with_priHit.txt",ios::app);
	  outfile<<moduleID<<" "<<pos<<" "<<adcsum<<" "<<size<<" "<<prim_ratio<<" "<<p_over_total<<" "<<b_over_total<<" "<<type<<" "<<ncross_talk<<" "<<peaktime<<endl;
	  outfile.close();
	}
#endif
#ifdef PRINTCLUSTER
	cout<<"pos: "<<pos<<" adcsum: "<<adcsum<<" peaktime: "<<peaktime<<" size: "<<size<<endl;
#endif
      }
#endif // MCDATA
#ifndef NDEBUG
      // Ensure hits are ordered by position (should be guaranteed by std::map)
      // no sense compare the postion after adding module concept
      //assert( prevHit == 0 or theHit->Compare(prevHit));
      prevHit = theHit;
#endif
    }

    // Undo amplitude splitting, if any, so fADCcor contains correct ADC values
    for( viter_t it = splits.begin(); it != splits.end(); ++it ) {
      fADCcor[*it] *= 2.0;
    }

    // Negative return value indicates potential problem
    if( nHits > fMaxHits )
      nHits = -nHits;
    fHits->Sort(fHits->GetEntriesFast());
    //ValueCmp(fHits);
    //  sort(*fHits,(*fHits).end());
    // cout<<"#########: "<<fCoveredPriHit<<"  "<<fTotalPriHit<<"  "<<(Double_t)(fCoveredPriHit)/fTotalPriHit<<endl;//getchar();
    
    //cout<<"DD: Number of hits in plane "<<" : "<<nHits<<endl;
    return nHits;
   
  }
  /*
  void GEMPlane::ValueCmp(TClonesArray* f)
  {
    Int_t size = f->GetEntriesFast();
    if(size==1) return;
    cout<<size<<"###"<<endl;
    for(int i=0;i<size;i++)
      {
	Hit* h = (Hit*)f->At(i);
	Hit* hi = (Hit*)f->At(i+1);
	//
	if(h->GetPos()>hi->GetPos())
	  {
	    //	cout<<h->GetPos()<<endl;getchar();
	  }
      }
  }
  */
  //_____________________________________________________________________________
  Hit* GEMPlane::AddHitImpl( Double_t pos )
  {
    // Make a dummy hit of the correct type at the given projection coordinate
    // and add it to the hit array

    assert( IsDummy() );

    GEMHit* theHit = 0;

    // Emulate parameters for dummy hits
    const UInt_t size = 1, type = 0;
    const Int_t moduleID = 0;
    const Double_t adcsum = 10.*fMinAmpl, adcmax_fit=adcsum/2, prim_ratio=0, p_over_total=0, b_over_total=0, peaktime=0, resolution = fResolution;
    Int_t stripl = 0, striph = 0;
  
#ifdef MCDATA
    const Int_t mctrack = 1, num_bg = 0;
    const Double_t mcpos = pos, mccharge = 0.0, mctime = 0.0;
    bool mc_data = fTracker->TestBit(Tracker::kMCdata);
    if( mc_data )
      // Monte Carlo data
      theHit = new( (*fHits)[GetNhits()] ) MCGEMHit( moduleID,
						     pos,
						     adcsum,
						     adcmax_fit,
						     prim_ratio,
						     p_over_total,
						     b_over_total,
						     0,
						     0,
						     0,
						     0,
						     peaktime,
						     size,
						     stripl,
						     striph,
						     type,
						     resolution,
						     this,
						     mctrack,
						     mcpos,
						     mccharge,
						     mctime,
						     num_bg
						     );
    else
#endif
      theHit = new( (*fHits)[GetNhits()] ) GEMHit( moduleID,
						   pos,
						   adcsum,
						   adcmax_fit,
						   peaktime,
						   size,
						   type,
						   resolution,
						   this
						   );
    return theHit;
  }

  //_____________________________________________________________________________
  Int_t GEMPlane::Decode( const THaEvData& evData )
  {
    // Convert evData to hits
    if( IsDummy() )
      // Special "decoding" for dummy planes
      return DummyDecode( evData );

    return GEMDecode( evData );
  }

  //_____________________________________________________________________________
  Int_t GEMPlane::DefineVariables( EMode mode )
  {
    // initialize global variables

    if( mode == kDefine && fIsSetup ) return kOK;
    fIsSetup = ( mode == kDefine );

    // Register variables in global list

    RVarDef vars[] = {
      { "strips",         "strips number of seen strips",     "fStrips" },
      { "nrawstrips",     "nstrips with decoder data",        "fNrawStrips" },
      { "nhitstrips",     "nstrips > 0",                      "fNhitStrips" },
      { "nstrips",        "Num strips with hits > adc.min",   "GetNsigStrips()" },
      { "hitocc",         "strips > 0 / n_all_strips",        "fHitOcc" },
      { "occupancy",      "nstrips / n_all_strips",           "fOccupancy" },
      //   { "strip.adcraw_max","raw strip ADC max",                "fADCraw" },
      //   { "strip.adc",      "raw strip ADC sum",                "fADC" },
      //   { "strip.adc_c",    "Pedestal-sub strip ADC sum",       "fADCcor" },
      //{ "strip.time",     "Leading time of strip signal (ns)","fHitTime" },
      // { "strip.good",     "Good pulse shape on strip",        "fGoodHit" },
      { "nhits",          "Num hits (clusters of strips)",    "GetNhits()" },
      { "noise",          "Noise level (avg below adc.min)",  "fDnoise" },
      { "ncoords",        "Num fit coords",                   "GetNcoords()" },
      { "coord.pos",      "Position used in fit (m)",         "fFitCoords.TreeSearch::FitCoord.fPos" },
      { "coord.trkpos",   "Track pos from projection fit (m)","fFitCoords.TreeSearch::FitCoord.fTrackPos" },
      { "coord.trkslope", "Track slope from projection fit",  "fFitCoords.TreeSearch::FitCoord.fTrackSlope" },
      { "coord.resid",    "Residual of trkpos (m)",           "fFitCoords.TreeSearch::FitCoord.GetResidual()" },
      { "coord.3Dpos",    "Crossing position of fitted 3D track (m)", "fFitCoords.TreeSearch::FitCoord.f3DTrkPos" },
      { "coord.3Dresid",  "Residual of 3D trkpos (m)",        "fFitCoords.TreeSearch::FitCoord.Get3DTrkResid()" },
      { "coord.3Dslope",  "Slope of fitted 3D track wrt projection",  "fFitCoords.TreeSearch::FitCoord.f3DTrkSlope" },
      { 0 }
    };
    Int_t ret = DefineVarsFromList( vars, mode );

    if( ret != kOK )
      return ret;

#ifdef MCDATA
    if( !fTracker->TestBit(Tracker::kMCdata) ) {
#endif
      // Non-Monte Carlo hit data
      RVarDef nonmcvars[] = {
	{ "hit.pos",  "Hit centroid (m)",      "fHits.TreeSearch::GEMHit.fPos" },
	{ "hit.adc",  "Hit ADC sum",           "fHits.TreeSearch::GEMHit.fADCsum" },
	{ "hit.max",  "Hit ADC max",           "fHits.TreeSearch::GEMHit.fADCmax" },
	{ "hit.time", "Hit peak time",         "fHits.TreeSearch::GEMHit.fPeaktime" },
	{ "hit.size", "Num strips ",           "fHits.TreeSearch::GEMHit.fSize" },
	{ "hit.type", "Hit analysis result",   "fHits.TreeSearch::GEMHit.fType" },
	{ 0 }
      };
      ret = DefineVarsFromList( nonmcvars, mode );
#ifdef MCDATA
    } else {
      // Monte Carlo hit data includes the truth information
      // For safety, we make sure that all hit variables are referenced with
      // respect to the MCGEMHit class and not just GEMHit - the memory layout
      // of classes under multiple inheritance might be implemetation-dependent
      RVarDef mcvars[] = {
	{ "hit.pos",   "Hit centroid (m)",      "fHits.TreeSearch::MCGEMHit.fPos" },
	{ "hit.adc",   "Hit ADC sum",           "fHits.TreeSearch::MCGEMHit.fADCsum" },
	{ "hit.max",  "Hit ADC max",            "fHits.TreeSearch::GEMHit.fADCmax" },
	{ "hit.time", "Hit peak time",          "fHits.TreeSearch::GEMHit.fPeaktime" },
	{ "hit.size",  "Num strips ",           "fHits.TreeSearch::MCGEMHit.fSize" },
	{ "hit.type",  "Hit analysis result",   "fHits.TreeSearch::MCGEMHit.fType" },
	{ "hit.mctrk", "MC track number",       "fHits.TreeSearch::MCGEMHit.fMCTrack" },
	{ "hit.mcpos", "MC track position (m)", "fHits.TreeSearch::MCGEMHit.fMCPos" },
	{ "hit.mccharge", "MC track charge",    "fHits.TreeSearch::MCGEMHit.fMCCharge" },//
	{ "hit.mctime","MC track time (s)",     "fHits.TreeSearch::MCGEMHit.fMCTime" },
	{ "hit.numbg", "MC num backgr hits",    "fHits.TreeSearch::MCGEMHit.fContam" },
	{ 0 }
      };
      ret = DefineVarsFromList( mcvars, mode );
    }
#endif
    return ret;
  }

  //_____________________________________________________________________________
  Int_t GEMPlane::End( THaRunBase* run )
  {
    // Write diagnostic histograms to file

    Plane::End( run );

#ifdef TESTCODE
    if( TestBit(kDoHistos) ) {
      assert( fADCMap );
      fADCMap->Write();
    }
#endif
    return 0;
  }

  //_____________________________________________________________________________
  Int_t GEMPlane::ReadDatabase( const TDatime& date )
  {
    // Read database

    static const char* const here = "ReadDatabase";

    // Read the database for the base class, quit if error
    Int_t status = ReadDatabaseCommon(date);
    if( status != kOK )
      return status;
    // Now read additional info needed for this derived class
    FILE* file = OpenFile( date );
    if( !file ) return kFileError;

    // Set defaults
    TString mapping;
    Int_t do_noise = 1, check_pulse_shape = 1;
    fMaxClusterSize = kMaxUInt;
    fMinAmpl   = 0.0;
    fSplitFrac = 0.0;
    fMapType   = kOneToOne;
    fMaxSamp   = 1;
    fChanMap.clear();
    fPed.clear();
    fAmplSigma = 0.36; // default, an educated guess
    ftmp_comm_range = 0.; // a very simple guess

    Int_t gbl = GetDBSearchLevel(fPrefix);
    try {
      const DBRequest request[] = {
	{ "nstrips",        &fNelem,          kInt,     0, 0, gbl },
	{ "strip.pos",      &fStart },
	{ "strip.pitch",    &fPitch,          kDouble,  0, 0, gbl },
	{ "maxclustsiz",    &fMaxClusterSize, kUInt,    0, 1, gbl },
	{ "maxsamp",        &fMaxSamp,        kUInt,    0, 1, gbl },
	{ "adc.min",        &fMinAmpl,        kDouble,  0, 1, gbl },
	{ "split.frac",     &fSplitFrac,      kDouble,  0, 1, gbl },
	{ "mapping",        &mapping,         kTString, 0, 1, gbl },
	{ "chanmap",        &fChanMap,        kIntV,    0, 1, gbl },
	{ "pedestal",       &fPed,            kFloatV,  0, 1 },
	{ "do_noise",       &do_noise,        kInt,     0, 1, gbl },
	{ "adc.sigma",      &fAmplSigma,      kDouble,  0, 1, gbl },
	{ "check_pulse_shape",&check_pulse_shape, kInt, 0, 1, gbl },
	{ "pedestal_sigma", &fpedestal_sigma,      kInt,0, 1, gbl },
	{ "tmp_pedestal_rms",&ftmp_pedestal_rms,  kDoubleV, 0, 1, gbl },
	{ "tmp_comm_range", &ftmp_comm_range,  kInt,    0, 1, gbl },
	{ "commonmode_groupsize", &fcModeSize,kInt,     0, 1, gbl },
	{ 0 }
      };
      status = LoadDB( file, date, request, fPrefix );
    }
    // Catch exceptions here so that we can close the file
    catch(...) {
      fclose(file);
      throw;
    }

    // Finished reading the database
    fclose(file);
    if( status != kOK )
      return status;

    // Sanity checks
    if( fNelem <= 0 or fNelem > kMaxNChan ) {
      Error( Here(here), "Invalid number of channels: %d. Must be > 0 and < %d. "
	     "Fix database or recompile code with a larger limit.",
	     fNelem, kMaxNChan );
      return kInitError;
    }

    // Dummy planes ignore all of the parameters that are checked below,
    // so we can return right here.
    if( IsDummy() ) {
      fIsInit = true;
      return kOK;
    }

    SetBit( kDoNoise, do_noise );
    SetBit( kCheckPulseShape, check_pulse_shape );

    SafeDelete(fADCraw);
    SafeDelete(fMCCharge);
    SafeDelete(fADC);
    SafeDelete(fHitTime);
    SafeDelete(fADCcor);
    SafeDelete(fGoodHit);
#ifdef MCDATA
    delete [] fMCHitInfo; fMCHitInfo = 0;
#endif

    // Allocate arrays. The only reason that these are parallel C-arrays is
    // that the global variable system still doesn't support arrays/vectors
    // of structures/objects.
    // Out of memory exceptions from here are caught in Tracker.cxx.
    fADCraw = new Float_t[fNelem];
    fADC = new Float_t[fNelem];
    fHitTime = new Float_t[fNelem];
    fMCCharge = new Float_t[fNelem];
    fADCcor = new Float_t[fNelem];
    fGoodHit = new Byte_t[fNelem];
    fSigStrips.reserve(fNelem);
    fStripsSeen.resize(fNelem);

#ifdef MCDATA
    if( fTracker->TestBit(Tracker::kMCdata) ) {
      fMCHitInfo = new TSBSMCHitInfo[fNelem];
      fMCHitList.reserve(fNelem);
    }
#endif

    TString::ECaseCompare cmp = TString::kIgnoreCase;
    if( !mapping.IsNull() ) {
      if( mapping.Length() >= 3 and
	  TString("one-to-one").BeginsWith(mapping,cmp) )
	fMapType = kOneToOne;
      else if( mapping.Length() >=3 and
	       TString("reverse").BeginsWith(mapping,cmp) )
	fMapType = kReverse;
      else if( mapping.Length() >= 5 and
	       mapping.BeginsWith(TString("gassiplex-adapter"),cmp) ) {
	if( fNelem > 240 ) {
	  Error( Here(here), "Gassiplex adapter mapping allows at most 240 "
		 "strips, but %d configured. Fix database.", fNelem );
	  return kInitError;
	}
	if( fNelem < 240 ) {
	  Warning( Here(here), "Gassiplex adapter mapping expects 240 "
		   "strips, but %d configured. Database may be misconfigured "
		   "(or you know what you are doing).", fNelem );
	}
	if( mapping.BeginsWith(TString("gassiplex-adapter-2"),cmp) ) {
	  fMapType = kGassiplexAdapter2;
	} else {
	  fMapType = kGassiplexAdapter1;
	}
      }
      else if( TString("table").CompareTo(mapping,cmp) == 0 ) {
	if( fChanMap.empty() ) {
	  Error( Here(here), "Channel mapping table requested, but no map "
		 "defined. Specify chanmap in database." );
	  return kInitError;
	}
	if( fChanMap.size() != static_cast<UInt_t>(fNelem) ) {
	  Error( Here(here), "Number of channel map entries (%u) must equal "
		 "number of strips (%d). Fix database.",
		 static_cast<unsigned int>(fChanMap.size()), fNelem );
	  return kInitError;
	}
	// check if entries in channel map are within range
	for( Vint_t::const_iterator it = fChanMap.begin();
	     it != fChanMap.end(); ++it ) {
	  if( (*it) < 0 or (*it) >= fNelem ) {
	    Error( Here(here), "Illegal chanmap entry: %d. Must be >= 0 and "
		   "< %d. Fix database.", (*it), fNelem );
	    return kInitError;
	  }
	}
	fMapType = kTable;
      //} else if( TString("mpd-apv").BeginsWith(mapping,cmp) ) {
      } else if( TString(mapping).BeginsWith("mpd-apv",cmp) ) {
        if( fChanMap.empty() ) {
          Error( Here(here), "Channel mapping %s requested, but no chanmap"
             " was defined. Ensure chanamp is specified in the database "
             " for this plane.",mapping.Data());
          return kInitError;
        }
        TString map_name("mpd-apv"); // Default
        TString map_format("crate slot mpd_id gem_id adc_id i2c pos invert");
        fStripsPerChannel = 128;
        fcModeSize = 1; // Reset fcModeSize back to 1 APV = 128 channels
        Int_t nvals = 8;
        Int_t napvs = fNelem/128;
        // Which version is specified?
        // v1(default): crate slot mpd_id gem_id adc_id i2c pos invert
        // v2: mpd_id gem_id adc_id i2c pos invert
        if( TString("mpd-apv2").CompareTo(mapping,cmp)  == 0 ||
            TString("mpd-apv2-ssp").CompareTo(mapping,cmp) == 0 ) {
          fMapType = kMPDAPV2;
          nvals = 6;
          map_name = "mpd-apv2";
          map_format = "mpd_id gem_id adc_id i2c pos invert";
        } else {
          fMapType = kMPDAPV;
        }
        if(TString(mapping).EndsWith("-nossp",cmp)) {
          fStripOrderedData = 1;
          map_name = mapping;
        }
        // Should be nvals*napv entries
        if(fChanMap.size()%nvals != 0 || ( fChanMap.size() !=
              UInt_t(napvs*nvals) ) ){
          Error( Here(here), "Channel mapping %s requested, which "
              " requires %d entries per channel in the following format:\n%s\n"
              " mpd_id gem_id adc_id i2c pos invert\n "
              " for a total of %d*%d = %d entries Fix database.",
              mapping.Data(),nvals,map_format.Data(),nvals,napvs,nvals*napvs);
          return kInitError;
        }
      } else {
	Error( Here(here), "Unknown channel mapping type %s. Fix database.",
	       mapping.Data() );
	return kInitError;
      }

    } else
      fChanMap.clear();

    Int_t nchan = fStripsPerChannel*fDetMap->GetTotNumChan();
    if( nchan != fNelem ) {
      Error( Here(here), "Number of detector map channels (%d) "
	     "disagrees with number of strips (%d)", nchan, fNelem );
      return kInitError;
    }
    fNapvs = fDetMap->GetTotNumChan()/fcModeSize;


    if( !fPed.empty() and fPed.size() != static_cast<UInt_t>(fNelem) ) {
      Error( Here(here), "Size of pedestal array (%u) must equal "
	     "number of strips (%d). Fix database.",
	     static_cast<unsigned int>(fPed.size()), fNelem );
      return kInitError;
    }

    // Sanity checks on fMaxSamp
    static const UInt_t max_maxsamp = 32; // arbitrary sanity limit
    if( fMaxSamp == 0 )
      fMaxSamp = 1;
    else if( fMaxSamp > max_maxsamp ) {
      Warning( Here(here), "Illegal maximum number of samples: %u. "
	       "Adjusted to maximum allowed = %u.", fMaxSamp, max_maxsamp );
      fMaxSamp = max_maxsamp;
    }
    if( fMaxSamp == 1 )
      ResetBit( kCheckPulseShape );

    if( fAmplSigma < 0.0 ) {
      Warning( Here(here), "Negative adc.sigma = %lf makes no sense. Adjusted "
	       "to positive.", fAmplSigma );
      fAmplSigma = TMath::Abs( fAmplSigma );
    }
    if( fAmplSigma < 0.001 ) {
      Warning( Here(here), "adc.sigma = %lf is extremely small. "
	       "Double-check database.", fAmplSigma );
    }

    fIsInit = true;
    return kOK;
  }
  //_____________________________________________________________________________
  void GEMPlane::Print( Option_t* opt ) const
  {
    // Print plane info

    Plane::Print( opt );
    //TODO: add detailed info
  }

  //_____________________________________________________________________________
  Double_t GEMPlane::GetAmplSigma( Double_t ampl ) const
  {
    // Return expected sigma of the hit amplitude distribution at the given
    // amplitude 'ampl'.
    // For starters, calculate a constant relative uncertainty
    // from a database parameter.

    return fAmplSigma * ampl;
  }

  //_____________________________________________________________________________

}

ClassImp(TreeSearch::GEMPlane)

///////////////////////////////////////////////////////////////////////////////

