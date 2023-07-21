#ifndef LASERPATH
#define LASERPATH
#include <iostream>
#include <vector>
#include <Eigen/Core>

namespace heat {
class Track {
  public:
    Eigen::Vector3d *p1;//origin
    Eigen::Vector3d *p2;//destination
    bool hasDeposition = false;
    double speed = 10;// mm/s
    double power = 100;// W
    double length = -1;

    Track(Eigen::Vector3d *p1, Eigen::Vector3d *p2,
        double speed, double power, bool hasDeposition) {
      this->p1 = p1;
      this->p2 = p2;
      this->speed = speed;
      this->power = power;
      this->hasDeposition = hasDeposition;
      length = (*p2-*p1).norm();
    }

    Eigen::Vector3d getSpeed() {
      return (*p2 - *p1).normalized()*speed;
    }
};

class Path {
  public:
    std::vector<Eigen::Vector3d> coordinates;
    std::vector<double> times;
    std::vector<Track> tracks;
    Track *currentTrack = NULL;

    Path( std::vector<Eigen::Vector3d> &coordinates,
          std::vector<double> &speeds,
          std::vector<double> &powers,
          std::vector<int> &arePrinting ) {
      this->coordinates = std::move( coordinates );
      tracks.reserve( this->coordinates.size()-1 );
      for (int itrack = 0; itrack < this->coordinates.size() -1; ++itrack) {
        tracks.push_back( Track( &this->coordinates[ itrack ], &this->coordinates[ itrack+1 ],
              speeds[ itrack + 1], powers[ itrack + 1], bool(arePrinting[ itrack+1 ]) ) );
      }
      setTimes();
      currentTrack = &tracks[0];
    }

    void setTimes() {
      times.resize(coordinates.size());
      times[0] = 0.0;
      for (int i = 0; i < tracks.size(); ++i) {
        Track *t = &tracks[i];
        times[i+1] = times[i] + t->length / t->speed;
      }
    }

    void updateCurrentTrack(double t) {
      auto isBigger = [t](double time){ return (t < time); };
      auto it = std::find_if( times.begin(), times.end(), isBigger );
      if (it == times.end()) {
        currentTrack = NULL;
      } else {
        currentTrack = &tracks[*it];
      }
      return;
    }
};
}
#endif
