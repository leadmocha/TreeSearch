//*-- Author :    Ole Hansen, Jefferson Lab   06-Feb-2012

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// SBS::GEMTracker                                                         //
//                                                                           //
// This is class GEMTracker in namespace SBS. It inherits from class       //
// GEMTracker in namespace TreeSearch. Although confusing, this is fine.     //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "SBSGEMTracker.h"
#include "SBSGEMPlane.h"
#include "SBSSpec.h"
#include "Helper.h"
#include "TClonesArray.h"

using namespace std;

namespace SBS {

typedef vector<Plane*>::size_type vrsiz_t;

//_____________________________________________________________________________
GEMTracker::GEMTracker( const char* name, const char* desc, THaApparatus* app )
  : TreeSearch::GEMTracker(name,desc,app), fSector(-1), fXOffset(0)
{
  // Constructor
}

//_____________________________________________________________________________
GEMTracker::~GEMTracker()
{
  // Destructor
}

//_____________________________________________________________________________
const char* GEMTracker::GetDBFileName() const
{
  // Return database file name prefix. For SBS trackers, this is the detector
  // name with any trailing sector number removed. In this way, all trackers
  // use the same database file.

  // fDBPrefix is set in MakePrefix()
  return fDBPrefix.Data();
}

//_____________________________________________________________________________
void GEMTracker::MakePrefix()
{
  // Set up name prefixes for global variables and database.
  // Global variables and database keys get the standard prefix,
  // e.g. "solid.tracker.3."
  // The database file name gets a different, shorter prefix, allowing
  // the trackers to share a single database file. For the example above,
  // "solid.tracker."

  TreeSearch::GEMTracker::MakePrefix();

  TString prefix( GetPrefix() );
  assert( prefix.EndsWith(".") );
  Int_t ndot = prefix.CountChar('.');
  if( ndot > 1 ) {
    prefix.Chop();
    if( !prefix.EndsWith(".") )
      prefix.Remove( prefix.Last('.')+1 );
    else
      Warning( Here("MakePrefix"), "Double dot in detector prefix = "
	       "\"%s\"?", GetPrefix() );
  }
  fDBPrefix = prefix;
}

//_____________________________________________________________________________
Plane* GEMTracker::MakePlane( const char* name, const char* description,
			      THaDetectorBase* parent ) const
{
  // Create an object of the plane class used by this implementation

  return new SBS::GEMPlane( name, description, parent );
}

//_____________________________________________________________________________
THaAnalysisObject::EStatus GEMTracker::PartnerPlanes()
{
  // Before doing PartnerPlanes() for GEM trackers, do additional geometry
  // calculations. This is done here because the planes need to have read
  // their database and be sorted by z.

  // TODO: protect against repeated calls? (Don't keep subtracting vectors)

  // Take the origin of the first plane as the origin of this Tracker and
  // make all the plane coordinates relative to the Tracker origin
  if( !fPlanes.empty() ) {
    // Bypass the const& - presumably we know what we're doing...
    TVector3& org = const_cast<TVector3&>( fPlanes.front()->GetOrigin() );
    fOrigin = org;
    org.SetXYZ( 0.0, 0.0, 0.0 );  // update first plane
    for( vrsiz_t iplane = 1; iplane < fPlanes.size(); ++iplane ) {
      Plane* thePlane = fPlanes[iplane];
      assert( thePlane->GetZ() >= fOrigin.Z() ); // else not sorted
      TVector3& other_org = const_cast<TVector3&>( thePlane->GetOrigin() );
      other_org -= fOrigin;
    }
  }
  // fOrigin needs to be given in the global (lab) reference frame
  fOrigin *= fRotation;

  // Now continue with standard PartnerPlanes() of the base class
  return TreeSearch::GEMTracker::PartnerPlanes();
}

//_____________________________________________________________________________
Int_t GEMTracker::ReadGeometry( FILE* file, const TDatime& date,
				Bool_t required )
{
  // Read basic geometry for a SBS GEM tracker sector
  //
  // The only geometry parameter needed for the a SBS Tracker (which
  // represents the collection of GEM trackers in a sector) is 'phi', the
  // central phi angle of the sector. The origin of the tracker coordinate
  // system will be determined automatically in PartnerPlanes() using the
  // positions of the GEM planes that are defined for this tracker.
  //
  // 'phi' is always required. The 'required' argument is ignored.

  //  static const char* const here = "ReadGeometry";

  DBRequest request[] = {
    { "xoff",    &fXOffset, kDouble, 0, 0, 0, "" },
    { 0 }
  };
  Int_t err = LoadDB( file, date, request, fPrefix );
  if( err )
    return err;

  // fOrigin will be set later PartnerPlanes (after Plane init)

  // Define this detector's rotation wrt the global coordinate system.
  // The rotation rotates the axes and not the vectors, hence it is a rotation
  // about z by _positive_ phi.
  // Any vector in this Tracker's reference frame can be transformed to the
  // global system simply by multiplication with the rotation matrix:
  // v_global = fRotation * v

  // fRotation.SetToIdentity();
  // fRotation.RotateZ( fPhi );
  // fInvRot.RotateZ( -fPhi );
  // assert( fRotation == fInvRot.Inverse() ); // might fail due to rounding
  
  return kOK;
}

//_____________________________________________________________________________
Int_t GEMTracker::NewTrackCalc( Int_t idx, THaTrack*, const TVector3& pos,
				const TVector3& dir, const FitRes_t& )
{
  // For every new track, convert track coordinates to the cylindrical system
  // appropriate for SBS. This is a temporary solution; the right way to do
  // this is to have the SBS spectrometer use its own track class.

  assert( dynamic_cast<SBSSpec*>(GetApparatus()) );
  SBSSpec* solid = static_cast<SBSSpec*>(GetApparatus());

  TClonesArray* trackInfo = solid->GetTrackInfo();

#ifdef MCDATA
  SBSTrackInfo* theInfo =
#endif
  new( (*trackInfo)[idx] ) SBSTrackInfo( fSector, pos, dir);

#ifdef MCDATA
  using TreeSearch::NumberOfSetBits;
  assert( static_cast<vrsiz_t>(idx) < fMCHitBits.size() );
  theInfo->SetMCHitBits( fMCHitBits[idx] );
  theInfo->SetNMCHits( NumberOfSetBits(theInfo->GetMCHitBits()) );
#endif

  return 0;
}

#ifdef MCDATA
//_____________________________________________________________________________
Int_t GEMTracker::FitMCPoints( Podd::MCTrack* mctrk ) const
{
  // SBS version of FitMCPoints to a MC track. In addition to fitting
  // the track, also reconstruct the track to the target.
  // Currently assumes stright tracks.

  Int_t npt = TreeSearch::Tracker::FitMCPoints( mctrk );
  if( npt < 0 )
    return npt;

  // See SBSSpec::FindVertices. Find distance of closest approach of
  // the fitted MC track to the beam (assumed to be exactly along z)
  // At this point, the track parameters have been converted to the lab frame
  Double_t* coef = mctrk->fMCFitPar;
  TVector3 g0( coef[0], coef[2], fOrigin.Z() );
  TVector3 g1( coef[1], coef[3], 1.0  );
  g1 = g1.Unit();
  TVector3 h0( 0, 0, 0  );
  TVector3 h1( 0, 0, 1. );
  Double_t gh = g1*h1;
  Double_t denom = 1.-gh*gh;
  if( TMath::Abs(denom) > 1e-6 ) {
    TVector3 D0 = g0-h0;
    Double_t tc = -D0*(g1-h1*gh)/denom;
    //Double_t sc =  D0*(h1-g1*gh)/denom;
    TVector3 vertex = g0 + g1*tc;
    // Save the results as additional fit results with mctrk
    coef[6] = vertex.X();
    coef[7] = vertex.Y();
    coef[8] = vertex.Z();
  }
  return npt;
}
#endif // MCDATA

///////////////////////////////////////////////////////////////////////////////

} // end namespace SBS

ClassImp(SBS::GEMTracker)
