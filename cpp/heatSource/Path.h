#ifndef LASERPATH
#define LASERPATH
#include <iostream>
#include <vector>
#include <Eigen/Core>

namespace heat {

enum TrackType {
  printing,
  cooling,
  dwelling,
  recoating,
};

class Track {
  public:
    Eigen::Vector3d *p0;//origin
    Eigen::Vector3d *p1;//destination
                        //these pointers are not ours
    TrackType type = cooling;
    bool isNewX = false;
    bool isNewY = false;
    bool isNewZ = false;
    double speed = 10;// mm/s
    double power = 100;// W
    double length = -1;
    double startTime, endTime;
    int index = -1;

    Track(Eigen::Vector3d *p0, Eigen::Vector3d *p1,
        double startTime, double endTime,
        double speed, double power, TrackType trackType, bool isNewX, bool isNewY, bool isNewZ, int index = -1 ) {
      this->p0 = p0;
      this->p1 = p1;
      this->startTime = startTime;
      this->endTime = endTime;
      //this->endTime = this->startTime + this->length / this->speed;
      this->speed = speed;
      this->power = power;
      this->type  = trackType;
      this->isNewX = isNewX;
      this->isNewY = isNewY;
      this->isNewZ = isNewZ;
      this->length = (*p1-*p0).norm();
      this->index = index;
    }

    Eigen::Vector3d getSpeed() const { return (*p1 - *p0).normalized()*speed; }
    bool isOver(double t) { return (t + 1e-7 >= endTime); }
};

class Path {
  public:
    std::vector<Track> tracks;
    std::vector<Eigen::Vector3d> coordinates;
    std::vector<double> times;

    Path( std::vector<Eigen::Vector3d> &coordinates,
          std::vector<double> &times,
          std::vector<double> &speeds,
          std::vector<double> &powers,
          std::vector<TrackType> &trackTypes ) {
      /*
       * Receiving coordinates, times and speeds makes system
       * overdetermined. Responsability of consistency is on gcode reader
       */
      this->coordinates = std::move( coordinates );
      this->times = std::move( times );
      tracks.reserve( this->coordinates.size()-1 );

      double currentX = -1e9;
      double currentY = -1e9;
      double currentZ = -1e9;
      for (int itrack = 0; itrack < this->coordinates.size() -1; ++itrack) {
        bool isNewX = false;
        bool isNewY = false;
        bool isNewZ = false;
        if (currentX /= this->coordinates[ itrack ][0]) {
          isNewX = true;
        }
        if (currentY /= this->coordinates[ itrack ][1]) {
          isNewY = true;
        }
        if (currentZ /= this->coordinates[ itrack ][2]) {
          isNewZ = true;
        }
        currentX = this->coordinates[ itrack ][0];
        currentY = this->coordinates[ itrack ][1];
        currentZ = this->coordinates[ itrack ][2];

        tracks.push_back( Track( &this->coordinates[ itrack ], &this->coordinates[ itrack+1 ],
              this->times[itrack], this->times[itrack+1],
              speeds[ itrack + 1], powers[ itrack + 1],
              trackTypes[ itrack+1 ], isNewX, isNewY, isNewZ, itrack ) );
      }
    }

    const Track* interpolateTrack(double t) const {
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

    bool isOver(double t) { return (t + 1e-7 >= this->times.back() ); }
};
}
#endif
