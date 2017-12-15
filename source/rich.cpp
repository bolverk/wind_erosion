#ifdef RICH_MPI
#include "source/mpi/MeshPointsMPI.hpp"
#endif
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "source/tessellation/geometry.hpp"
#include "source/newtonian/two_dimensional/hdsim2d.hpp"
#include "source/tessellation/tessellation.hpp"
#include "source/newtonian/common/hllc.hpp"
#include "source/newtonian/common/ideal_gas.hpp"
#include "source/tessellation/VoronoiMesh.hpp"
#include "source/newtonian/two_dimensional/spatial_distributions/uniform2d.hpp"
#include "source/newtonian/two_dimensional/point_motions/eulerian.hpp"
#include "source/newtonian/two_dimensional/point_motions/lagrangian.hpp"
#include "source/newtonian/two_dimensional/point_motions/round_cells.hpp"
#include "source/newtonian/two_dimensional/source_terms/zero_force.hpp"
#include "source/newtonian/two_dimensional/geometric_outer_boundaries/SquareBox.hpp"
#include "source/newtonian/test_2d/random_pert.hpp"
#include "source/newtonian/two_dimensional/diagnostics.hpp"
#include "source/misc/simple_io.hpp"
#include "source/misc/mesh_generator.hpp"
#include "source/newtonian/test_2d/main_loop_2d.hpp"
#include "source/newtonian/two_dimensional/hdf5_diagnostics.hpp"
#include "source/tessellation/shape_2d.hpp"
#include "source/newtonian/test_2d/piecewise.hpp"
#include "source/newtonian/two_dimensional/simple_flux_calculator.hpp"
#include "source/newtonian/two_dimensional/simple_cell_updater.hpp"
#include "source/newtonian/two_dimensional/simple_extensive_updater.hpp"
#include "source/newtonian/two_dimensional/stationary_box.hpp"
#include "source/tessellation/right_rectangle.hpp"
#include "source/newtonian/test_2d/clip_grid.hpp"
#include "source/newtonian/test_2d/multiple_diagnostics.hpp"
#include "source/misc/vector_initialiser.hpp"
#include "source/newtonian/test_2d/consecutive_snapshots.hpp"
#include "source/newtonian/two_dimensional/condition_action_sequence.hpp"
#include "source/newtonian/two_dimensional/source_terms/cylindrical_complementary.hpp"
#include "source/newtonian/two_dimensional/source_terms/SeveralSources.hpp"

using namespace std;
using namespace simulation2d;

namespace {

  vector<Vector2D> centered_hexagonal_grid(double r_min,
					   double r_max)
  {
    const vector<double> r_list = arange(0,r_max,r_min);
    vector<Vector2D> res;
    for(size_t i=0;i<r_list.size();++i){
      const size_t angle_num = max<size_t>(6*i,1);
      vector<double> angle_list(angle_num,0);
      for(size_t j=0;j<angle_num;++j)
	angle_list.at(j) = 2*M_PI*static_cast<double>(j)/static_cast<double>(angle_num);
      for(size_t j=0;j<angle_num;++j)
	res.push_back(r_list.at(i)*Vector2D(cos(angle_list.at(j)),
					    sin(angle_list.at(j))));
    }
    return res;
  }

  vector<Vector2D> centered_logarithmic_spiral(double r_min,
					       double r_max,
					       double alpha,
					       const Vector2D& center)
  {
    const double theta_max = log(r_max/r_min)/alpha;
    const vector<double> theta_list = 
      arange(0,theta_max,2*M_PI*alpha/(1-0.5*alpha));
  
    vector<double> r_list(theta_list.size(),0);
    for(size_t i=0;i<r_list.size();++i)
      r_list.at(i) = r_min*exp(alpha*theta_list.at(i));
  
    vector<Vector2D> res(r_list.size());
    for(size_t i=0;i<res.size();++i)
      res[i] = center+r_list[i]*Vector2D(cos(theta_list.at(i)),
					 sin(theta_list.at(i)));
    return res;
  }

  vector<Vector2D> complete_grid(double r_inner,
				 double r_outer,
				 double alpha)
  {
    const vector<Vector2D> inner = 
      centered_hexagonal_grid(r_inner*alpha*2*M_PI,
			      r_inner);
    const vector<Vector2D> outer =
      centered_logarithmic_spiral(r_inner,
				  r_outer,
				  alpha,
				  Vector2D(0,0));
    return join(inner, outer);
  }

#ifdef RICH_MPI

  vector<Vector2D> process_positions(const SquareBox& boundary)
  {
    const Vector2D lower_left = boundary.getBoundary().first;
    const Vector2D upper_right = boundary.getBoundary().second;
	int ws=0;
	MPI_Comm_size(MPI_COMM_WORLD,&ws);
    return RandSquare(ws,lower_left.x,upper_right.x,lower_left.y,upper_right.y);
  }

#endif

  vector<ComputationalCell> calc_init_cond(const Tessellation& tess)
  {
    vector<ComputationalCell> res(static_cast<size_t>(tess.GetPointNo()));
    for(size_t i=0;i<res.size();++i){
      const Vector2D& r = tess.GetMeshPoint(static_cast<int>(i));
      res[i].density = abs(r) > 0.1 ? 1 : 1e2;
      res[i].pressure = 1e-9;
      res[i].velocity = abs(r) > 0.1 ? Vector2D(1,0) : Vector2D(0,0);
    }
    return res;
  }

  class LeftBoundaryEdge : public ConditionActionSequence::Condition
  {
  public:

    LeftBoundaryEdge(void) {}

    pair<bool, bool> operator()
    (const Edge& edge,
     const Tessellation& tess,
     const vector<ComputationalCell>& /*cells*/,TracerStickerNames const& /*ts*/) const
    {
      if (tess.GetOriginalIndex(edge.neighbors.first) == tess.GetOriginalIndex(edge.neighbors.second))
	{
	  if (abs(edge.vertices.first.x-edge.vertices.second.x)<abs(edge.vertices.first.y-edge.vertices.second.y))
	    {
	      if(edge.neighbors.second<tess.GetPointNo())
		{
		  Vector2D const& r = tess.GetMeshPoint(edge.neighbors.second);
		  if(r.x>edge.vertices.first.x)
		    return pair<bool, bool>(true, false);
		}
	      else
		{
		  Vector2D const& r = tess.GetMeshPoint(edge.neighbors.first);
		  if(r.x>edge.vertices.first.x)
		    return pair<bool, bool>(true, true);
		}
	    }
	  return pair<bool, bool>(false, false);
	}
      return pair<bool, bool>(false, false);
    }
  };

  void conserved_to_extensive
  (const Conserved& c, const ComputationalCell& /*cell*/, Extensive &res)
  {
    res.mass = c.Mass;
    res.momentum = c.Momentum;
    res.energy = c.Energy;
  }

  ComputationalCell calc_ghost(void)
  {
    ComputationalCell res;
    res.density = 1;
    res.pressure = 1e-9;
    res.velocity = Vector2D(1, 0);
    res.stickers.push_back(false);
    return res;
  }


  class ConstantGhost : public ConditionActionSequence::Action
  {
  public:

    ConstantGhost(const ComputationalCell& ghost,
		  const RiemannSolver& rs) :
      ghost_(ghost), rs_(rs) {}

    void operator()
    (const Edge& edge,
     const Tessellation& tess,
     const Vector2D& /*edge_velocity*/,
     const vector<ComputationalCell>& cells,
     const EquationOfState& eos,
     const bool aux,
     Extensive &res,double /*time*/,TracerStickerNames const& ts) const
    {
      if (aux)
	assert(edge.neighbors.first < tess.GetPointNo());
      else
	assert(edge.neighbors.second < tess.GetPointNo());
      const Vector2D p = normalize
	(edge.vertices.second -
	 edge.vertices.first);
      const Vector2D n = normalize
	(remove_parallel_component
	 (aux ?
	  edge.vertices.first - tess.GetMeshPoint(edge.neighbors.first) :
	  tess.GetMeshPoint(edge.neighbors.second) - edge.vertices.first,
	  p));
      const double v = 0;
      const pair<ComputationalCell, ComputationalCell> cc_left_righ =
	aux ?
	pair<ComputationalCell, ComputationalCell>
	(cells.at(static_cast<size_t>(edge.neighbors.first)), ghost_) :
	pair<ComputationalCell, ComputationalCell>
	(ghost_, cells.at(static_cast<size_t>(edge.neighbors.second)));
      const pair<Primitive, Primitive> left_right =
	pair<Primitive, Primitive>
	(convert_to_primitive(cc_left_righ.first, eos, ts),
	 convert_to_primitive(cc_left_righ.second, eos, ts));
      const Conserved c = rotate_solve_rotate_back
	(rs_,
	 left_right.first,
	 left_right.second,
	 v, n, p);
      conserved_to_extensive(c, ghost_, res);
    }

  private:
    const ComputationalCell ghost_;
    const RiemannSolver& rs_;
  };

  class SelectiveCenterGravity: public Acceleration
  {
  public:
    SelectiveCenterGravity
    (double M,
     double Rmin,
     const Vector2D& centre):
      M_(M),
      Rmin_(Rmin),
      centre_(centre) {}

    Vector2D operator()
    (const Tessellation& tess,
     const vector<ComputationalCell>& cells,
     const vector<Extensive>& /*fluxes*/,
     const double /*time*/,
     const int point,
     TracerStickerNames const& /*tracerstickernames*/) const
    {
      if(cells.at(static_cast<size_t>(point)).density<1e-8)
	return Vector2D(0,0);
      const Vector2D pos(tess.GetCellCM(point)-centre_);
      const double r = abs(pos);
      if(abs(pos-Vector2D(0,-0.04))<0.04)
	return Vector2D(0,0);
      return (-1)*pos*M_/(r*r*r+Rmin_*Rmin_*Rmin_);
    }

  private:
    const double M_;
    const double Rmin_;
    const Vector2D centre_;
  };

  class PressureFloor: public CellUpdater
  {
  public:

    PressureFloor(void) {}

    vector<ComputationalCell> operator()
      (const Tessellation& tess,
       const PhysicalGeometry& /*pg*/,
       const EquationOfState& eos,
       vector<Extensive>& extensives,
       const vector<ComputationalCell>& old,
       const CacheData& cd,
       TracerStickerNames const& tracerstickernames) const
    {
      size_t N = static_cast<size_t>(tess.GetPointNo());
      vector<ComputationalCell> res(N, old[0]);
      for(size_t i=0;i<N;++i){
	Extensive& extensive = extensives[i];
	const double volume = cd.volumes[i];
	res[i].density = extensive.mass / volume;
	if(res[i].density<0)
	  throw UniversalError("Negative density");
	res[i].velocity = extensive.momentum / extensive.mass;
	const double energy = extensive.energy / extensive.mass - 
	  0.5*ScalarProd(res[i].velocity, res[i].velocity);
	try{
	  if(energy>0)
	    res[i].pressure = eos.de2p(res[i].density,
				       energy,
				       res[i].tracers,
				       tracerstickernames.tracer_names);
	  else
	    res[i].pressure = 1e-9;
	}
	catch(UniversalError& eo){
	  eo.AddEntry("cell density", res[i].density);
	  eo.AddEntry("cell energy", energy);
	  throw;
	}
      }	
      return res;
    }
  };

  class SimData
  {
  public:

    SimData(void):
      pg_(Vector2D(0,0), Vector2D(1,0)),
      width_(1e0),
      outer_(-width_,width_,width_,1e-6),
#ifdef RICH_MPI
	  vproc_(process_positions(outer_),outer_),
		init_points_(SquareMeshM(50,50,vproc_,outer_.getBoundary().first,outer_.getBoundary().second)),
		tess_(vproc_,init_points_,outer_),
#else
      init_points_(clip_grid
		   (RightRectangle(Vector2D(-width_,1e-6), Vector2D(width_, width_)),
		    complete_grid(0.12,
				  2*width_,
				  0.005))),
		tess_(init_points_, outer_),
#endif
      eos_(5./3.),
      point_motion_(),
      //      bpm_(),
      //      point_motion_(bpm_, eos_),
      sb_(),
      rs_(),
      gravity_acc_(1e-2,
		   0.1,
		   Vector2D(0,0)),
      gravity_force_(gravity_acc_),
      geom_force_(pg_.getAxis()),
      force_(VectorInitialiser<SourceTerm*>
	     (&gravity_force_)
	     (&geom_force_)()),
      tsf_(0.3),
      fc_(VectorInitialiser<pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*> >
	  (pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*>
	   (new LeftBoundaryEdge, new ConstantGhost(calc_ghost(), rs_)))
	  (pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*>
	   (new IsBoundaryEdge, new FreeFlowFlux(rs_)))
	  (pair<const ConditionActionSequence::Condition*, const ConditionActionSequence::Action*>
	   (new IsBulkEdge, new RegularFlux(rs_)))()),
      eu_(),
      cu_(),
      sim_(
#ifdef RICH_MPI
		  vproc_,
#endif
		  tess_,
	   outer_,
	   pg_,
	   calc_init_cond(tess_),
	   eos_,
	   point_motion_,
	   sb_,
	   force_,
	   tsf_,
	   fc_,
	   eu_,
	   cu_) {}

    hdsim& getSim(void)
    {
      return sim_;
    }

  private:
    const CylindricalSymmetry pg_;
    const double width_;
    const SquareBox outer_;
#ifdef RICH_MPI
	VoronoiMesh vproc_;
#endif
    const vector<Vector2D> init_points_;
    VoronoiMesh tess_;
    const IdealGas eos_;
#ifdef RICH_MPI
    //Eulerian point_motion_;
    //	Lagrangian point_motion_;
#else
    Eulerian point_motion_;
    //    Lagrangian bpm_;
    //    RoundCells point_motion_;
#endif
    const StationaryBox sb_;
    const Hllc rs_;
    SelectiveCenterGravity gravity_acc_;
    ConservativeForce gravity_force_;
    CylindricalComplementary geom_force_;
    SeveralSources force_;
    const SimpleCFL tsf_;
    const ConditionActionSequence fc_;
    const SimpleExtensiveUpdater eu_;
    //    const SimpleCellUpdater cu_;
    const PressureFloor cu_;
    hdsim sim_;
  };

  class WriteCycle: public DiagnosticFunction
  {
  public:

    WriteCycle(const string& fname):
      fname_(fname) {}

    void operator()(const hdsim& sim)
    {
      write_number(sim.getCycle(),fname_);
    }

  private:
    const string fname_;
  };

  class AtmosphereCooling: public Manipulate
  {
  public:

    AtmosphereCooling(void) {}

    void operator()(hdsim& sim)
    {
      vector<ComputationalCell>& cells = sim.getAllCells();
      for(size_t i=0;i<cells.size();++i){
	const Vector2D& r = sim.getTessellation().GetMeshPoint(static_cast<int>(i));
	if(r.y<0)
	  continue;
	ComputationalCell& cell = cells.at(i);
	const double q = 1 - 1e-6*r.y;
	//	cell.density = max(1e-6 + q*(cell.density-1e-6),
	//			   1e-6);
	cell.pressure = max(1e-4 + q*(cell.pressure-1e-4),
			    1e-4);
	cell.velocity *= q;
      }
      sim.recalculateExtensives();
    }
  };
}

int main(void)
{
#ifdef RICH_MPI
	MPI_Init(NULL,NULL);
	int rank=0;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
  SimData sim_data;
  hdsim& sim = sim_data.getSim();
  write_snapshot_to_hdf5(sim, "output/initial.h5");

  const double tf = 1e1;
  SafeTimeTermination term_cond(tf,1e6);
  MultipleDiagnostics diag
  (VectorInitialiser<DiagnosticFunction*>()
   (new ConsecutiveSnapshots(new ConstantTimeInterval(tf/100),
			     new Rubric("output/snapshot_",".h5")))
   (new WriteTime("time.txt"))
   (new WriteCycle("cycle.txt"))());
  AtmosphereCooling ac;
  main_loop(sim,
	    term_cond,
	    &hdsim::TimeAdvance,
	    &diag,
	    &ac);
	    

#ifdef RICH_MPI
  write_snapshot_to_hdf5(sim, "process_"+int2str(rank)+"_final"+".h5");
  MPI_Finalize();
#else
  write_snapshot_to_hdf5(sim, "output/final.h5");
#endif


  return 0;
}

