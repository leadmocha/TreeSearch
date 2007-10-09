#ifndef ROOT_TreeSearch_Projection
#define ROOT_TreeSearch_Projection

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// TreeSearch::Projection                                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "THaAnalysisObject.h"
#include "TreeWalk.h"
#include "TMath.h"
#include <vector>
#include <cassert>

using std::vector;
using std::string;

class THaDetectorBase;

namespace TreeSearch {

  class Hitpattern;
  class PatternTree;
  class WirePlane;

  class Projection : public THaAnalysisObject {
  public:

    Projection( Int_t type, const char* name, Double_t angle,
		THaDetectorBase* parent = 0 );
    Projection( const Projection& orig );
    const Projection& operator=( const Projection& rhs );
    virtual ~Projection();

    virtual void    Clear( Option_t* opt="" );
    virtual Int_t   Decode( const THaEvData& );
    EStatus         InitLevel2( const TDatime& date );
    virtual void    Print( Option_t* opt="" ) const;
    void            Reset();

    Int_t           FillHitpattern();
    Int_t           TreeSearch();

    void            AddPlane( WirePlane* wp );
    Int_t           GetType() const { return fType; }
    UInt_t          GetNlevels()  const { return fNlevels; }
    Double_t        GetMaxSlope() const { return fMaxSlope; }
    Double_t        GetZsize() const;
    Double_t        GetWidth() const { return fWidth; }
    Double_t        GetAngle() const;
    Double_t        GetSinAngle() const { return fSinAngle; }
    Double_t        GetCosAngle() const { return fCosAngle; }
    Hitpattern*     GetHitpattern() const { return fHitpattern; }
    UInt_t          GetNplanes() const
    { return static_cast<UInt_t>( fPlanes.size() ); }

    void            SetMaxSlope( Double_t m ) { fMaxSlope = m; }
    void            SetPatternTree( PatternTree* pt ) { fPatternTree = pt; }
    void            SetWidth( Double_t width ) { fWidth = width; }

    //FIXME: for testing
    vector<TreeSearch::WirePlane*>& GetListOfPlanes() { return fPlanes; }

  protected:
    Int_t               fType;        // Type of plane (u,v,x,y...)
    vector<WirePlane*>  fPlanes;      // Wire planes of this type
    UInt_t              fNlevels;     // Number of levels of search tree
    Double_t            fMaxSlope;    // Maximum physical track slope (0=perp)
    Double_t            fWidth;       // Width of tracking region (m)
    Double_t            fSinAngle;    // Sine of wire angle
    Double_t            fCosAngle;    // Cosine of wire angle

    Hitpattern*         fHitpattern;  // Hitpattern of current event
    PatternTree*        fPatternTree; // Precomputed template database

    THaDetectorBase*    fDetector;    //! Parent detector

    //FIXME: TEST
    UInt_t  n_test, n_found;
    Double_t search_time;
    void  SetAngle( Double_t a );

    // Podd interface
    virtual Int_t ReadDatabase( const TDatime& date );
    virtual Int_t DefineVariables( EMode mode = kDefine );
    virtual const char* GetDBFileName() const;
    virtual void MakePrefix();

    // NodeVisitor class for comparing patterns in the tree with the
    // hitpattern. Matches represent candidates for track roads and are
    // added to the list of roads for further analysis
    class ComparePattern {
    public:
      ComparePattern( Projection* proj ) 
	: fProj(proj) { assert(fProj); fHitpattern = fProj->GetHitpattern(); }
      TreeWalk::ETreeOp operator() ( const NodeDescriptor& nd );
    private:
      const Hitpattern* fHitpattern;    // Hitpattern to compare to
      Projection* fProj;                // Projection we belong to
//       std::list<Road>&  fRoads;         // List of results
    };
    friend class ComparePattern;

    ClassDef(Projection,0)  // A track projection plane
  };

  //___________________________________________________________________________
  inline
  Double_t Projection::GetAngle() const {
    // Return wire angle in rad, normalized to [-pi,pi]
    Double_t a = TMath::ASin(fSinAngle);
    if( fCosAngle < 0.0 )
      return (fSinAngle > 0.0) ? TMath::TwoPi() - a : -TMath::TwoPi() - a;
  
    return a;
  }

///////////////////////////////////////////////////////////////////////////////

} // end namespace TreeSearch

#endif
