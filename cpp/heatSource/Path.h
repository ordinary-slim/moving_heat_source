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
    bool hasDeposition = false;
    double speed = 10;// mm/s
    double power = 100;// W
    double length = -1;
    double startTime;
    double endTime;

    Track(Eigen::Vector3d *p0, Eigen::Vector3d *p1,
        double speed, double power, bool hasDeposition) {
      this->p0 = p0;
      this->p1 = p1;
      this->speed = speed;
      this->power = power;
      this->hasDeposition = hasDeposition;
      length = (*p1-*p0).norm();
    }

    Eigen::Vector3d getSpeed() {
      return (*p1 - *p0).normalized()*speed;
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
        t->startTime = times[i];
        t->endTime = times[i+1];
      }
    }

    Track* interpolateTrack(double t) {
      double tol = 1e-7;
      Track* track = NULL;
      auto it = std::find_if( times.begin(), times.end(),
          [t, tol](double time){ return (t <= time+tol); } );
      if (it != times.end()) {
        int idxTrack = it - times.begin() - 1;
        track = &tracks[idxTrack];
      }
      return track;
    }

    Eigen::Vector3d interpolatePosition(double t) {
      Track* track = interpolateTrack( t );
      if (track == NULL) {
        throw std::invalid_argument("Invalid time does not belong to path.");
      } else {
        double trackFraction = (t - track->startTime)/(track->endTime - track->startTime);
        return (*track->p0 + trackFraction*(*track->p1 - *track->p0));
      }
    }
};
}
#endif
