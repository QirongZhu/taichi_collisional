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


struct Particle {
    double X[3];
    uint64_t id;
    uint64_t hp_key;
    uint64_t top_node;
    bool done=false;
};

typedef std::vector<Particle> Particles;

extern Domain domain;
