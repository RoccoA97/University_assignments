#include "AnalysisSteering.h"

using namespace std;

AnalysisSteering::AnalysisSteering( const AnalysisInfo* info ):
 aInfo( info ) {
}


AnalysisSteering::~AnalysisSteering() {
}


// function to be called at execution start / end
void AnalysisSteering::update( const AnalysisInfo::AnalysisStatus& status ) {
  switch ( status ) {
  case AnalysisInfo::begin:
    beginJob();
    break;
  case AnalysisInfo::end:
    endJob();
    break;
  }
  return;
}

