#include "TimerStatistic.h"

using namespace ProtoMol::Report;

namespace ProtoMol {

  //_________________________________________________________________ TimerStatistic
  Timer TimerStatistic::timer[static_cast<int>(LAST)-static_cast<int>(FIRST)];
  bool TimerStatistic::myIsParallel = false;

  Report::MyStreamer& operator<<(Report::MyStreamer& os, 
				 const TimerStatistic& ){
    os.setf(std::ios::showpoint|std::ios::fixed);
    os.precision(5);
    os << "wall: "<<TimerStatistic::timer[TimerStatistic::WALL].getTime()
       << ", run: "<<TimerStatistic::timer[TimerStatistic::RUN].getTime().getRealTime()<<"[s]"
       << ", integration: "<<TimerStatistic::timer[TimerStatistic::INTEGRATOR].getTime().getRealTime()<<"[s]"
       << ", forces: "<<TimerStatistic::timer[TimerStatistic::FORCES].getTime().getRealTime()<<"[s]";
    if(TimerStatistic::isParallel()){
      os << ", com "<<TimerStatistic::timer[TimerStatistic::COMMUNICATION].getTime().getRealTime() 
	 << "[s], idle "<<TimerStatistic::timer[TimerStatistic::IDLE].getTime().getRealTime()<<"[s]";
    }
    os.reset();
    return os;
  }
}
