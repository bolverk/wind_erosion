#include <iostream>
#include "newtonian/two_dimensional/hdsim2d.hpp"
#include "tessellation/VoronoiMesh.hpp"
#include "newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "newtonian/two_dimensional/physical_geometry.hpp"

using namespace std;

namespace {
  vector<Vector2D> create_grid(void)
  {
    vector<Vector2D> res;
    for(double x=-1;x<1;x+=0.1){
      for(double y=-1;y<1;y+=0.1)
	res.push_back(Vector2D(x,y));
    }
  return res;
  }

  vector<ComputationalCell> create_init_cond
  (const Tessellation& tess)
  {
    vector<ComputationalCell> res(tess.GetPointNo());
    for(size_t i=0;i<res.size();++i){
      const Vector2D& r = tess.GetMeshPoint(static_cast<int>(i));
      res.at(i).density = 1;
      if(r.x*r.x+r.y*r.y>0.2*0.2)
	res.at(i).pressure = 1;
      else
	res.at(i).pressure = 10;
      res.at(i).velocity = Vector2D(0,0);
    }
    return res;
  }

  class SimData
  {
  public:

    SimData(void):
      bottom_left_(-1,-1),
      upper_right_(1,1),
      box_(bottom_left_, upper_right_),
      tess_(create_grid(),
	    box_),
      pg_() {}

  private:
    const Vector2D bottom_left_;
    const Vector2D upper_right_;
    const SquareBox box_;
    VoronoiMesh tess_;
    const SlabSymmetry pg_;
  };
}

int main(void)
{
  return 0;
}
