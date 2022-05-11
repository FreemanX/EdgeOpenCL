#ifndef EDCL_INCLUDE_TRACKABLE_H_
#define EDCL_INCLUDE_TRACKABLE_H_
#include <unordered_map>

static int TRACK_CNT = 0;

class Trackable {
 public:
  Trackable() {
	trackNum = TRACK_CNT++;
  }

  Trackable(const Trackable &t) {
	*this = t;
	this->trackNum = TRACK_CNT++;
  }

  int getTrackNum() const { return trackNum; }
 protected:
  int trackNum;
};

#define ADD_TO_TRACK_MAP(map, trackable_ptr) map[trackable_ptr->getTrackNum()] = trackable_ptr

#endif //EDCL_INCLUDE_TRACKABLE_H_
