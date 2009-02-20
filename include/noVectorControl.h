#include "transmissionModel.h"

//! Base transmission model, as used in Phase A
class NoVectorControl :
public TransmissionModel { 
 public:

  NoVectorControl() {};
  ~NoVectorControl() {};
  //! initialise the main simulation 
  void initMainSimulation (int populationSize); 	

 private:

  int i;

};
