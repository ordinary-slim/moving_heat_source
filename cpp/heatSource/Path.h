#ifndef LASERPATH
#define LASERPATH
#include <iostream>
#include <vector>
#include <Eigen/Core>

namespace heat {
class Track {
  public:
    Eigen::Vector3d *p0;//origin
    Eigen::Vector3d *p1;//destination
                        //these pointers are not ours
    bool hasDeposition = false;
    double speed = 10;// mm/s
    double power = 100;// W
    double length = -1;
    double startTime, endTime;

    Track(Eigen::Vector3d *p0, Eigen::Vector3d *p1,
        double startTime, double speed, double power, bool hasDeposition) {
      this->p0 = p0;
      this->p1 = p1;
      this->startTime = startTime;
      this->speed = speed;
      this->power = power;
      this->hasDeposition = hasDeposition;
      this->length = (*p1-*p0).norm();
      this->endTime = this->startTime + this->length / this->speed;
    }

    Eigen::Vector3d getSpeed() const { return (*p1 - *p0).normalized()*speed; }
    bool isOver(double t) { return (t + 1e-7 >= endTime); }
};

class Path {
  public:
    std::vector<Track> tracks;
    std::vector<Eigen::Vector3d> coordinates;
    std::vector<double> times;
    double endTime;

    Path( std::vector<Eigen::Vector3d> &coordinates,
          std::vector<double> &speeds,
          std::vector<double> &powers,
          std::vector<int> &arePrinting ) {
      this->coordinates = std::move( coordinates );
      this->times.resize( this->coordinates.size() );
      this->times[0] = 0.0;
      tracks.reserve( this->coordinates.size()-1 );
      for (int itrack = 0; itrack < this->coordinates.size() -1; ++itrack) {
        tracks.push_back( Track( &this->coordinates[ itrack ], &this->coordinates[ itrack+1 ], times[itrack],
              speeds[ itrack + 1], powers[ itrack + 1], bool(arePrinting[ itrack+1 ]) ) );
        times[itrack+1] = tracks[itrack].endTime;
      }
      endTime = times.back();
    }

    const Track* interpolateTrack(double t) {
      double tol = 1e-7;
      const Track* track = NULL;
      auto it = std::find_if( times.begin(), times.end(),
          [t, tol](double time){ return (t <= time+tol); } );
      if (it != times.end()) {
        int idxTrack = it - times.begin() - 1;
        if (idxTrack == - 1) { idxTrack = 0; }
        track = &tracks[idxTrack];
      }
      return track;
    }

    Eigen::Vector3d interpolatePosition(double t) {
      const Track* track = interpolateTrack( t );
      if (track == NULL) {
        throw std::invalid_argument("Invalid time does not belong to path.");
      } else {
        double trackFraction = (t - track->startTime)/(track->endTime - track->startTime);
        return (*track->p0 + trackFraction*(*track->p1 - *track->p0));
      }
    }

    bool isOver(double t) { return (t + 1e-7 >= endTime); }
};
}
#endif
