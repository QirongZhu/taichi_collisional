/*
  definitions of structs
*/
struct Domain {

  double Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;

  void setDomain(double xmin, double xmax, 
		 double ymin, double ymax, 
		 double zmin, double zmax) {
    this->Xmin = xmin;
    this->Xmax = xmax;
    this->Ymin = ymin;
    this->Ymax = ymax;
    this->Zmin = zmin;
    this->Zmax = zmax;
  }

};


extern Domain domain;
