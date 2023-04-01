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
vector<vector<int>> getVertexSets_Dd( Eigen::VectorXi localCon_D0,
    int d, ElementType entD_elType) {

  vector<vector<int>> vertexSets;

  ElementType entd_elType = getIncidentElType( entD_elType, d );
  int nnodes = localCon_D0.size();
  int windowSize = getNnodesElType(entd_elType);

  vector<int> vSet;
  vSet.resize( windowSize );

  switch (entD_elType){
    case line2:
      for (int inode = 0; inode < nnodes; inode++) {
        vSet[0] = localCon_D0[inode];
        vertexSets.push_back( vSet );
      }
      break;
    case triangle3: case quad4:
      for (int ipoin = 0; ipoin < nnodes; ipoin++) {
        for (int jpoin = 0; jpoin < windowSize; jpoin++) {
          vSet[jpoin] = localCon_D0( (ipoin + jpoin)%localCon_D0.size() );
        }
        std::sort( vSet.begin(), vSet.end() );
        vertexSets.push_back( vSet );
      }
      break;
    default:
      cout << "ERROR: getBoundary is not implemented yet for " << entD_elType << "element." << endl;
      exit(-1);
  }
  return vertexSets;
}

int get_dIndex( vector<vector<int>> &auxConnec_Dd,
    vector<vector<int>> &auxConnec_d0,
    int j, vector<int> &v){

  for (int idcell : auxConnec_Dd[j]) {
    vector<int> vtest = auxConnec_d0[idcell];
    if (vtest == v) {
      return idcell;
    }
  }
  return -1;
}

Connectivity auxCon2Con(vector<vector<int>> auxCon, int oDim, int tDim, int nels_oDim, int nels_tDim,
    ElementType oelType, ElementType telType) {

  auto begin = std::chrono::steady_clock::now();

  //Compute max tDim per oDim
  int max_tDim = -1;
  for (int ioDim = 0; ioDim < auxCon.size(); ioDim++) {
    max_tDim = std::max( max_tDim, int(auxCon[ioDim].size()) );
  }
  Eigen::MatrixXi con;
  con.resize( auxCon.size(), max_tDim );
  con.setOnes();
  con *= -1;

  for (int i_odim = 0; i_odim < auxCon.size(); i_odim++) {
    for (int i_tdim = 0; i_tdim < auxCon[i_odim].size(); i_tdim++) {
      con( i_odim, i_tdim ) = auxCon[i_odim][i_tdim];
    }
  }

  auto end = std::chrono::steady_clock::now();
  std::cout << "Converting aux con -> con took " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
  return Connectivity(con, oDim, tDim, nels_oDim, nels_tDim, oelType, telType);
}

void addEntry( vector<vector<int>> &auxCon, int i, int j ) {
  auto it = find( auxCon[i].begin(), auxCon[i].end(), j);
  if (it == auxCon[i].end() ) {
    auxCon[i].push_back( j );
  }
}

Connectivity transpose(Connectivity inCon) {
  vector<vector<int>> tCon;
  tCon.resize( inCon.nels_tDim );
  for (int ient_odim = 0; ient_odim < inCon.nels_oDim; ient_odim++ ) {
    for ( auto p_ient_tdim = inCon.con.row(ient_odim).begin();
        (p_ient_tdim != inCon.con.row(ient_odim).end())&&(*p_ient_tdim != -1);
          ++p_ient_tdim ) {
      addEntry( tCon, *p_ient_tdim, ient_odim );
    }
  }
  Connectivity outCon = auxCon2Con( tCon, inCon.tDim, inCon.oDim, inCon.nels_tDim, inCon.nels_oDim, inCon.telType, inCon.oelType );
  return outCon;
}

Connectivity intersect(Connectivity inCon1, Connectivity inCon2) {
  vector<vector<int>> intersecCon;
  intersecCon.resize( inCon1.nels_oDim );

  if (inCon1.oDim != inCon2.tDim ) {
    cout << "not implemented yet" << endl;
    exit(-1);
  }
  for (int enti = 0; enti < inCon1.nels_oDim; enti++) {
    int idx_entj = 0;
    for ( auto pentk = inCon1.con.row(enti).begin();
        (pentk != inCon1.con.row(enti).end())&&(*pentk != -1);
          ++pentk ) {
      for ( auto pentj = inCon2.con.row(*pentk).begin();
          (pentj != inCon2.con.row(*pentk).end())&&(*pentj != -1);
            ++pentj ) {
        if (enti != *pentj){
          addEntry( intersecCon, enti, *pentj );
        }
      }
    }
  }
  Connectivity outCon = auxCon2Con( intersecCon, inCon1.oDim, inCon2.tDim, inCon1.nels_oDim, inCon2.nels_tDim, inCon1.oelType, inCon2.telType );
  return outCon;
}

std::tuple<Connectivity, Connectivity> build(int d, Connectivity DO_connec, Connectivity DD_connec) {
  int D = DO_connec.oDim;
  vector<vector<int>> auxConnec_Dd;
  vector<vector<int>> auxConnec_d0;

  auxConnec_Dd.resize( DD_connec.nels_oDim );

  int k = 0;
  vector<vector<int>> Vi, Vj;
  vector<bool> newEntity;
  for (int icell = 0; icell < DD_connec.nels_oDim; icell++){
    Vi = getVertexSets_Dd(DO_connec.con.row(icell), d, DO_connec.oelType);
    newEntity.resize(Vi.size());
    fill(newEntity.begin(), newEntity.end(), true);

    for ( auto pjcell = DD_connec.con.row(icell).begin();
        (pjcell != DD_connec.con.row(icell).end())&&(*pjcell != -1);
          ++pjcell ) {
      if (*pjcell > icell) continue;

      Vj = getVertexSets_Dd(DO_connec.con.row(*pjcell), d, DO_connec.oelType);

      for (int ient = 0; ient < Vi.size(); ient++) {
        vector<int> vi = Vi[ient];
        for (vector<int> vj : Vj ) {
          if (vi==vj) {
            int l = get_dIndex(auxConnec_Dd, auxConnec_d0, *pjcell, vi);
            addEntry( auxConnec_Dd, icell, l );
            newEntity[ient] = false;
          }
        }
      }

    }

    for (int ient = 0; ient < Vi.size(); ient++) {
      if (newEntity[ient]) {
        addEntry( auxConnec_Dd, icell, k );
        auxConnec_d0.push_back( Vi[ient] );
        k++;
      }
    }
  }

  Connectivity connec_d0 = auxCon2Con( auxConnec_d0, d, 0, auxConnec_d0.size(),
      DO_connec.nels_tDim,
      getIncidentElType( DO_connec.oelType, d ), point1 );
  Connectivity connec_Dd = auxCon2Con( auxConnec_Dd, DD_connec.oDim, d,
      DD_connec.nels_oDim, connec_d0.nels_oDim, DO_connec.oelType,
      getIncidentElType( DO_connec.oelType , d));
  return std::make_tuple( connec_d0, connec_Dd );
}
}
