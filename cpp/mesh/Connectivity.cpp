#include "Connectivity.h"
#include <map>
#include <algorithm>
#include <iostream>
#include <Eigen/Core>
#include <vector>
#include <algorithm>
#include <chrono>
#include "ElementTypes.h"

using namespace std;

namespace mesh
{
int get_dIndex( vector<vector<unsigned int>> &connec_Dd,
    vector<vector<unsigned int>> &connec_d0,
    int j, vector<unsigned int> &vsorted){

  for (int idcell : connec_Dd[j]) {
    vector<unsigned int> vtest( connec_d0[idcell].size() );
    partial_sort_copy( connec_d0[idcell].begin(), connec_d0[idcell].end(), vtest.begin(), vtest.end() );
    if (vtest == vsorted) {
      return idcell;
    }
  }
  return -1;
}

void addEntry( vector<vector<unsigned int>> &auxCon, int i, int j ) {
  auto it = find( auxCon[i].begin(), auxCon[i].end(), j);
  if (it == auxCon[i].end() ) {
    auxCon[i].push_back( j );
  }
}

Connectivity transpose(Connectivity inCon) {
  vector<vector<unsigned int>> tCon;
  tCon.resize( inCon.nels_tDim );
  for (int ient_odim = 0; ient_odim < inCon.nels_oDim; ient_odim++ ) {
    const vector<unsigned int> *localCon = inCon.getLocalCon( ient_odim );
    for ( auto p_ient_tdim = localCon->begin(); p_ient_tdim != localCon->end(); ++p_ient_tdim ) {
      addEntry( tCon, *p_ient_tdim, ient_odim );
    }
  }
  
  return Connectivity( tCon, inCon.tDim, inCon.oDim, inCon.nels_tDim, inCon.nels_oDim, inCon.telType, inCon.oelType );
}

Connectivity intersect(Connectivity inCon1, Connectivity inCon2) {
  vector<vector<unsigned int>> intersecCon;
  intersecCon.resize( inCon1.nels_oDim );

  if (inCon1.oDim != inCon2.tDim ) {
    cout << "not implemented yet" << endl;
    exit(-1);
  }
  for (int enti = 0; enti < inCon1.nels_oDim; enti++) {
    int idx_entj = 0;
    const vector<unsigned int>* locCon1 = inCon1.getLocalCon(enti);
    for ( auto pentk = locCon1->begin(); (pentk != locCon1->end())&&(*pentk != -1); ++pentk ) {
      const vector<unsigned int>* locCon2 = inCon2.getLocalCon(*pentk);
      for ( auto pentj = locCon2->begin(); (pentj != locCon2->end())&&(*pentj != -1); ++pentj ) {
        if (enti != *pentj){
          addEntry( intersecCon, enti, *pentj );
        }
      }
    }
  }
  return Connectivity( intersecCon, inCon1.oDim, inCon2.tDim, inCon1.nels_oDim, inCon2.nels_tDim, inCon1.oelType, inCon2.telType );
}

std::tuple<Connectivity, Connectivity> buildBoundaryConnectivities( Connectivity DO_connec, Connectivity DD_connec) {
  int D = DO_connec.oDim;
  int d = D - 1;
  vector<vector<unsigned int>> connec_Dd;
  vector<vector<unsigned int>> connec_d0;

  connec_Dd.resize( DD_connec.nels_oDim );

  int k = 0;
  vector<vector<unsigned int>> Vi, Vj;
  vector<bool> newEntity;
  for (int icell = 0; icell < DD_connec.nels_oDim; icell++){
    Vi = getFacetVertexSets(*DO_connec.getLocalCon(icell), DO_connec.oelType);
    newEntity.resize(Vi.size());
    fill(newEntity.begin(), newEntity.end(), true);

    const vector<unsigned int>* locCon = DD_connec.getLocalCon(icell);
    for ( auto pjcell = locCon->begin(); (pjcell != locCon->end())&&(*pjcell != -1); ++pjcell ) {
      if (*pjcell > icell) continue;

      Vj = getFacetVertexSets(*DO_connec.getLocalCon(*pjcell), DO_connec.oelType);

      for (int ient = 0; ient < Vi.size(); ient++) {
        vector<unsigned int> vi = Vi[ient];
        vector<unsigned int> sorted_vi(vi.size());
        partial_sort_copy(begin(vi), end(vi),
            begin(sorted_vi), end(sorted_vi));
        for (vector<unsigned int> vj : Vj ) {
          vector<unsigned int> sorted_vj(vj.size());
          partial_sort_copy(begin(vj), end(vj),
              begin(sorted_vj), end(sorted_vj));
          if (sorted_vi==sorted_vj) {
            int l = get_dIndex(connec_Dd, connec_d0, *pjcell, sorted_vi);
            addEntry( connec_Dd, icell, l );
            newEntity[ient] = false;
          }
        }
      }

    }

    for (int ient = 0; ient < Vi.size(); ient++) {
      if (newEntity[ient]) {
        addEntry( connec_Dd, icell, k );
        connec_d0.push_back( Vi[ient] );
        k++;
      }
    }
  }

  return std::make_tuple(
      Connectivity( connec_d0, d, 0, connec_d0.size(),
      DO_connec.nels_tDim,
      getFacetElType( DO_connec.oelType ), point1 ),
      Connectivity( connec_Dd, DD_connec.oDim, d,
      DD_connec.nels_oDim, connec_d0.size(),
      DO_connec.oelType, getFacetElType( DO_connec.oelType) )
      );
}
}
