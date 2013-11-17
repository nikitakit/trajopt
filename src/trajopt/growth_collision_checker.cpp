#include "trajopt/collision_checker.hpp"
#include "boost/multi_array.hpp"
#include <cassert>

#include <vector>
#include <boost/foreach.hpp>
#include "utils/logging.hpp"

#include "utils/eigen_conversions.hpp"
#include "sco/expr_op_overloads.hpp"

using namespace util;
using namespace std;
using namespace trajopt;
using namespace OpenRAVE;

// typedef boost::multi_array<float, 3> array_type;
// typedef array_type::index index;

class ModelCache {
public:
  bool have_model;
  ModelPtr m;
  std::vector<Var> cachedVars;
  ModelCache(): have_model(false) {}
  virtual ~ModelCache() {}

  ModelPtr GetModel() {
    if (!have_model) {
      m = createModel();
      have_model = true;
    }
    return m;
  }

  Var addVar() {
    if (! cachedVars.empty()) {
      Var res = cachedVars.back();
      cachedVars.pop_back();
      return res;
    }
    return m->addVar("TODO:name", 0.0, 1e26);
  }

  void ReturnModel(std::vector<Cnt> cnts) {
    m->removeCnts(cnts);
    cachedVars = m->getVars();
  }
};

static ModelCache modelcache;

class GrowthCollisionChecker : public CollisionChecker {
  // btCollisionWorld* m_world;
  // btBroadphaseInterface* m_broadphase;
  // btCollisionDispatcher* m_dispatcher;
  // btCollisionConfiguration* m_coll_config;
  // typedef map<const OR::KinBody::Link*, CollisionObjectWrapper*> Link2Cow;
  // Link2Cow m_link2cow;
  double m_contactDistance;
  vector<KinBodyPtr> m_prevbodies;
  typedef std::pair<const KinBody::Link*, const KinBody::Link*> LinkPair;
  set< LinkPair > m_excludedPairs;
  // Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> m_allowedCollisionMatrix;

  bool m_initialized;

public:
  GrowthCollisionChecker(OR::EnvironmentBaseConstPtr env) : CollisionChecker(env) {
    m_initialized = false;
    SetContactDistance(.05);
  }
  ~GrowthCollisionChecker() {}

  ///////// public interface /////////
  virtual void SetContactDistance(float distance) {
    m_contactDistance = distance;
  }
  virtual double GetContactDistance() {return m_contactDistance;}
  virtual void PlotCollisionGeometry(vector<OpenRAVE::GraphHandlePtr>& handles) {
    LOG_DEBUG("Plotting collision geometry is not implemented");
    // NOT IMPLEMENTED
  }
  virtual void ExcludeCollisionPair(const KinBody::Link& link0, const KinBody::Link& link1) {
    m_excludedPairs.insert(LinkPair(&link0, &link1));
    // COW *cow0 = GetCow(&link0), *cow1 = GetCow(&link1);
    // if (cow0 && cow1) m_allowedCollisionMatrix(cow0->m_index, cow1->m_index) = 0;
  }
  // collision checking
  virtual void AllVsAll(vector<Collision>& collisions) {
    LOG_WARN("AllVsAll is not implemented!!!");
  }
  virtual void LinksVsAll(const vector<KinBody::LinkPtr>& links, vector<Collision>& collisions, short filterMask) {
    LOG_WARN("Link collision checking not implemented!!!");
    for (int i=0; i < links.size(); ++i) {
      LinkVsAll(*links[i], collisions, filterMask);
    }

    if (! m_initialized){
      m_initialized = true;
      IgnoreZeroStateSelfCollisions();
    }
  }

// #if __OPENRAVE_VERSION_MINOR__ <= 8
//     #define GT_Box KinBody::Link::GEOMPROPERTIES::GeomBox 
//     #define GT_Sphere KinBody::Link::GEOMPROPERTIES::GeomSphere 
//     #define GT_Cylinder KinBody::Link::GEOMPROPERTIES::GeomCylinder 
//     #define GT_TriMesh KinBody::Link::GEOMPROPERTIES::GeomTrimesh 
//     #define TriMesh KinBody::Link::TRIMESH
// #endif

  virtual void LinkVsAll(const KinBody::Link& link, vector<Collision>& collisions, short filterMask) {
    CheckLinkVsAll(link, false, Transform(), Transform(), collisions, filterMask);
  }
  virtual void CastLinkVsAll(const KinBody::Link& link, const Transform& tbefore, const Transform& tafter, vector<Collision>& collisions, short filterMask) {
    CheckLinkVsAll(link, true, tbefore, tafter, collisions, filterMask);
  }

  virtual void ContinuousCheckTrajectory(const TrajArray& traj, Configuration& rad, vector<Collision>&) {
    LOG_ERROR("ContinuousCheckTrajectory is only used by unit tests, and therefore unimplemented");
    throw std::runtime_error("not implemented");
  }
  virtual void CastVsAll(Configuration& rad, const vector<KinBody::LinkPtr>& links, const DblVec& startjoints, const DblVec& endjoints, vector<Collision>& collisions) {
    LOG_WARN("Cast collision checking not implemented!!!");

    Configuration::SaverPtr saver = rad.Save();
    rad.SetDOFValues(startjoints);
    int nlinks = links.size();
    vector<Transform> tbefore(nlinks), tafter(nlinks);
    for (int i=0; i < nlinks; ++i) {
      tbefore[i] = links[i]->GetTransform();
    }
    rad.SetDOFValues(endjoints);

    for (int i=0; i < nlinks; ++i) {
      tafter[i] = links[i]->GetTransform();
    }
    rad.SetDOFValues(startjoints);

    for (int i=0; i < nlinks; ++i) {
      CastLinkVsAll(*links[i].get(), tbefore[i], tafter[i], collisions, /*filtermask=*/-1); // TODO: figure out filtermask
    }
    LOG_DEBUG("CastVsAll checked %li links and found %li collisions", links.size(), collisions.size());
  }

  // ------------------------------------------------------------------

  class GrowthObject {
  public:
    Eigen::Array<float,1,3> min_pt;
    Eigen::Array<float,1,3> max_pt;
    Eigen::Array<float,1,3> center;
    Eigen::Array<float,1,3> sum;
    Eigen::Matrix<float,Eigen::Dynamic,3> points;

    int actual_size;

    GrowthObject(): center(0,0,0), sum(0,0,0),
                    min_pt(-INFINITY, -INFINITY, -INFINITY),
                    max_pt(INFINITY, INFINITY, INFINITY),
                    points(), actual_size(0){}
    virtual ~GrowthObject() {}

    void AddPoint(OR::Vector pt) {
      Eigen::Array<float,1,3> ept(pt.x, pt.y, pt.z);

      if (actual_size > 6) {
        Eigen::Array<float, Eigen::Dynamic, 3> distances = (points.array().rowwise() - ept);
        float dist = distances.matrix().rowwise().norm().minCoeff();
        if (dist < 0.1) {
          return;
        }
      }

      if (actual_size == 0) {
        points.conservativeResize(64, Eigen::NoChange);
      } else if (actual_size == points.rows()) {
        points.conservativeResize(points.rows() * 2, Eigen::NoChange);
      }

      points.row(actual_size) = ept.matrix();
      actual_size += 1;
      // points.bottomRows(1) = 
      sum += ept;
      min_pt = min_pt.cwiseMin(ept);
      max_pt = max_pt.cwiseMax(ept);
    }

    void Recenter() {
      points.conservativeResize(actual_size, Eigen::NoChange);
      center += sum / points.rows();
      sum = Eigen::Array<float,1,3>(0,0,0);
      points.rowwise() -= center.matrix();
      min_pt -= center;
      max_pt -= center;
    }

    OR::AABB ComputeAABB() {
      Eigen::Array<float,1,3> extents = (max_pt - min_pt) / 2.0;
      extents -= Eigen::Array<float,1,3>(0.1, 0.1, 0.1);
      extents = extents.cwiseMax(Eigen::Array<float,1,3>(0,0,0));
      return OR::AABB(OR::Vector(center(0,0), center(0,1), center(0,2)),
                      OR::Vector(extents(0,0), extents(0,1), extents(0,2)));
    }

    Eigen::Matrix<Var, 1, Eigen::Dynamic> CreateVars3(ModelPtr m) {
      Eigen::Matrix<Var, 1, Eigen::Dynamic> vars(1, points.rows());
      for (int i=0; i<points.rows(); i++) {
        vars(0,i) = modelcache.addVar();
      }
      return vars;
    }

    Eigen::Matrix<Var, 1, Eigen::Dynamic> CreateVars(ModelPtr m) {
      Eigen::Matrix<Var, 1, Eigen::Dynamic> vars(1, points.rows());
      for (int i=0; i<points.rows(); i++) {
        vars(0,i) = m->addVar("TODO:variable_name", 0.0, 1e26); // TODO:infinity
      }
      return vars;
    }

    double Radius() {
      return points.rowwise().norm().maxCoeff(); // TODO: is this right????
      // return points.colwise().norm().maxCoeff(); // TODO: is this right????
    }

    Vector Center() {
      return OR::Vector(center(0,0), center(0,1), center(0,2));
    }

    static double CalculateDistance3(GrowthObject& a, GrowthObject& b, OR::Vector& direction_result) {
      ModelPtr m = modelcache.GetModel();
      Var sigma = modelcache.addVar();
      m->update();
      m->setObjective(AffExpr(sigma));

      Eigen::Matrix<Var, 1, Eigen::Dynamic> vars_a = a.CreateVars(m);
      Eigen::Matrix<Var, 1, Eigen::Dynamic> vars_b = b.CreateVars(m);
      m->update();

      std::vector<Cnt> cnts;

      Eigen::Matrix<AffExpr, 1, 3> mtx;
      for (int i=0; i<3; i++) {
        for (int j=0; j<vars_a.cols(); j++) {
          mtx(0,i) += vars_a(0,j) * a.points(j,i);
        }
        mtx(0, i) += a.center(0,i);
        for (int j=0; j<vars_b.cols(); j++) {
          mtx(0,i) -= vars_b(0,j) * b.points(j,i);
        }
        mtx(0, i) -= b.center(0,i);

        cnts.push_back(m->addEqCnt(mtx(0,i), "TODO:eq-cnt"));
      }

      AffExpr a_sum = vars_a.cast<AffExpr>().sum();
      cnts.push_back(m->addEqCnt(a_sum - sigma, "shape_a"));

      AffExpr b_sum = vars_b.cast<AffExpr>().sum();
      cnts.push_back(m->addEqCnt(b_sum - sigma, "shape_b"));

      m->update();


      m->optimize();

      vector<double> values = m->getVarValues(m->getVars());

      Eigen::Matrix<double, 1, 3> direction(0,0,0);
      for (int i=0; i<vars_a.cols(); i++) {
        for (int j=0; j<3; j++) {
          direction(0, j) = vars_a(0, i).value(values) * a.points(i, j);
        }
      }

      // std::cout << m->getVarValue(vars_a(0,0)) << endl
      //           << m->getVarValue(vars_a(0,1)) << endl
      //           << m->getVarValue(vars_a(0,2)) << endl
      //           << m->getVarValue(vars_a(0,3)) << endl
      //           << m->getVarValue(vars_a(0,4)) << endl
      //           << m->getVarValue(vars_a(0,5)) << endl
      //           << m->getVarValue(vars_a(0,6)) << endl
      //           << m->getVarValue(vars_a(0,7)) << endl;


      // std::cout << "Sigma is " << m->getVarValue(sigma) << endl;
      // std::cout << "Direction is " << direction(0,0) << "," << direction(0, 1) << "," << direction(0,2) << endl;

      double distance = m->getVarValue(sigma) - 1.0;
      direction_result = OR::Vector(direction(0,0), direction(0,1), direction(0,2));

      distance *= a.Radius() + b.Radius();
      direction_result *= 0.15;

      modelcache.ReturnModel(cnts);

      return distance;
      // return 100;
    }


    static double CalculateDistance2(GrowthObject& a, GrowthObject& b, OR::Vector& direction_result) {
      ModelPtr m = createModel();
      Var sigma = m->addVar("sigma", 0.0, 1e26); // TODO:infinity
      m->update();
      m->setObjective(AffExpr(sigma));

      Eigen::Matrix<Var, 1, Eigen::Dynamic> vars_a = a.CreateVars(m);
      Eigen::Matrix<Var, 1, Eigen::Dynamic> vars_b = b.CreateVars(m);
      m->update();

      Eigen::Matrix<AffExpr, 1, 3> mtx;
      for (int i=0; i<3; i++) {
        for (int j=0; j<vars_a.cols(); j++) {
          mtx(0,i) += vars_a(0,j) * a.points(j,i);
        }
        mtx(0, i) += a.center(0,i);
        for (int j=0; j<vars_b.cols(); j++) {
          mtx(0,i) -= vars_b(0,j) * b.points(j,i);
        }
        mtx(0, i) -= b.center(0,i);

        m->addEqCnt(mtx(0,i), "TODO:eq-cnt");
      }

      AffExpr a_sum = vars_a.cast<AffExpr>().sum();
      m->addEqCnt(a_sum - sigma, "shape_a");

      AffExpr b_sum = vars_b.cast<AffExpr>().sum();
      m->addEqCnt(b_sum - sigma, "shape_b");

      m->update();


      m->optimize();

      vector<double> values = m->getVarValues(m->getVars());

      Eigen::Matrix<double, 1, 3> direction(0,0,0);
      for (int i=0; i<vars_a.cols(); i++) {
        for (int j=0; j<3; j++) {
          direction(0, j) = vars_a(0, i).value(values) * a.points(i, j);
        }
      }

      // std::cout << m->getVarValue(vars_a(0,0)) << endl
      //           << m->getVarValue(vars_a(0,1)) << endl
      //           << m->getVarValue(vars_a(0,2)) << endl
      //           << m->getVarValue(vars_a(0,3)) << endl
      //           << m->getVarValue(vars_a(0,4)) << endl
      //           << m->getVarValue(vars_a(0,5)) << endl
      //           << m->getVarValue(vars_a(0,6)) << endl
      //           << m->getVarValue(vars_a(0,7)) << endl;


      // std::cout << "Sigma is " << m->getVarValue(sigma) << endl;
      // std::cout << "Direction is " << direction(0,0) << "," << direction(0, 1) << "," << direction(0,2) << endl;

      double distance = m->getVarValue(sigma) - 1.0;
      direction_result = OR::Vector(direction(0,0), direction(0,1), direction(0,2));

      distance *= a.Radius() + b.Radius();
      direction_result *= 0.15;

      return distance;
      // return 100;
    }
  };

  virtual void CheckLinkVsAll(const KinBody::Link& link, bool cast, const Transform& tbefore, const Transform& tafter, vector<Collision>& collisions, short filterMask) {
    if (link.GetGeometries().empty()) return;

    const std::vector<boost::shared_ptr<OpenRAVE::KinBody::Link::Geometry> > & link_geometries=link.GetGeometries();
    std::vector<GrowthObject> growth_objects;
    BOOST_FOREACH(const boost::shared_ptr<OpenRAVE::KinBody::Link::Geometry>& geom, link_geometries) {
          switch (geom->GetType()) {
          case OpenRAVE::GT_Box:
          case OpenRAVE::GT_Sphere:
          case OpenRAVE::GT_Cylinder:
          case OpenRAVE::GT_TriMesh:
          {
            OpenRAVE::TriMesh mesh;
            GrowthObject obj;
            if (cast) {
              mesh = geom->GetCollisionMesh();
              mesh.ApplyTransform(tbefore * geom->GetTransform());
              for (size_t i = 0; i < mesh.indices.size(); i++) {
                obj.AddPoint(mesh.vertices[mesh.indices[i]]);
              }

              mesh = geom->GetCollisionMesh();
              mesh.ApplyTransform(tafter * geom->GetTransform());
              for (size_t i = 0; i < mesh.indices.size(); i++) {
                obj.AddPoint(mesh.vertices[mesh.indices[i]]);
              }
            } else {
              mesh = geom->GetCollisionMesh();
              mesh.ApplyTransform(link.GetTransform() * geom->GetTransform());
              for (size_t i = 0; i < mesh.indices.size(); i++) {
                obj.AddPoint(mesh.vertices[mesh.indices[i]]);
              }
            }

            obj.Recenter();
            growth_objects.push_back(obj);
            break;
          }
          default:
          {
            assert(0 && "unrecognized collision shape type");
            break;
          }
          }
    }

    vector<OR::KinBodyPtr> bodies;
    m_env->GetBodies(bodies);
    BOOST_FOREACH(const KinBodyPtr& body, bodies) {
      const vector<OR::KinBody::LinkPtr> body_links = body->GetLinks();
      BOOST_FOREACH(const OR::KinBody::LinkPtr& body_link, body_links) {
        // if (body_link->GetGeometries().size() == 0) {
        //   LOG_WARN("ignoring link %s", body_link->GetName().c_str());
        //   continue;
        // }

        if (m_excludedPairs.count(LinkPair(&link, body_link.get())) > 0) {
          // LOG_WARN("Excluded pair found %s %s", link.GetName().c_str(), body_link->GetName().c_str());
          continue;
        }

        const std::vector<boost::shared_ptr<OpenRAVE::KinBody::Link::Geometry> > & geometries=body_link->GetGeometries();
        BOOST_FOREACH(const boost::shared_ptr<OpenRAVE::KinBody::Link::Geometry>& geom, geometries) {
          switch (geom->GetType()) {
          case OpenRAVE::GT_Box:
          case OpenRAVE::GT_Sphere:
          case OpenRAVE::GT_Cylinder:
          case OpenRAVE::GT_TriMesh:
          {
            GrowthObject obj;
            OpenRAVE::TriMesh mesh = geom->GetCollisionMesh();
            mesh.ApplyTransform(body_link->GetTransform() * geom->GetTransform());
            for (size_t i = 0; i < mesh.indices.size(); i++) {
              obj.AddPoint(mesh.vertices[mesh.indices[i]]);
            }
            obj.Recenter();

            OR::Vector g_direction;
            double g_distance = 100000000.0;
            OR::Vector g_center;
            BOOST_FOREACH(GrowthObject& other_obj, growth_objects) {
              if (! OR::geometry::AABBCollision(obj.ComputeAABB(), other_obj.ComputeAABB())) {
                  continue;
                }
              OR::Vector direction;
              // TODO: which order?
              double distance = GrowthObject::CalculateDistance(obj, other_obj, direction);
              if (distance < g_distance) {
                g_distance = distance;
                g_direction = direction;
                g_center = other_obj.Center();
              }
            }

            if (g_distance < 0) {
              LOG_WARN("Collision detected %s %s", link.GetName().c_str(), body_link->GetName().c_str());


              // link a, link b, ptA, ptB, gradient(B to A), distance
              collisions.push_back(Collision(&link, body_link.get(),
                                             g_center, g_center+g_direction, g_direction, g_distance)
                );
            }
            break;
          }
          default:
          {
            assert(0 && "unrecognized collision shape type");
            break;
          }
          }
        }
      }
    }
  }

  AABB AABBUnion(const AABB& ab1, const AABB& ab2)
    {
      Vector min1 = ab1.pos - ab1.extents;
      Vector min2 = ab2.pos - ab2.extents;
      Vector max1 = ab1.pos + ab1.extents;
      Vector max2 = ab2.pos + ab2.extents;

      Vector minimum = Vector(fmin(min1.x, min2.x), fmin(min1.y, min2.y), fmin(min1.z, min2.z));
      Vector maximum = Vector(fmax(max1.x, max2.x), fmax(max1.y, max2.y), fmax(max1.z, max2.z));

      Vector pos = 0.5 * (minimum + maximum);
      Vector extents = 0.5 * (maximum - minimum);

      return AABB(pos, extents);

      // return AABB( 0.5f * (ab1.pos + ab2.pos), 0.5f * (ab1.extents + ab2.extents +
      //                                                    Vector(fabsf(ab2.pos.x-ab1.pos.x), fabsf(ab2.pos.y-ab1.pos.y), fabsf(ab2.pos.z-ab1.pos.z))) );
    }



  ////
  ///////

  // CollisionObjectWrapper* GetCow(const KinBody::Link* link) {
  //   Link2Cow::iterator it = m_link2cow.find(link);
  //   return (it == m_link2cow.end()) ? 0 : it->second;
  // }
  // void SetCow(const KinBody::Link* link, COW* cow) {m_link2cow[link] = cow;}
  // void LinkVsAll_NoUpdate(const KinBody::Link& link, vector<Collision>& collisions, short filterMask);
  // void UpdateBulletFromRave();
  // void AddKinBody(const OR::KinBodyPtr& body);
  // void RemoveKinBody(const OR::KinBodyPtr& body);
  // void AddAndRemoveBodies(const vector<OR::KinBodyPtr>& curVec, const vector<OR::KinBodyPtr>& prevVec, vector<KinBodyPtr>& addedBodies);
  // bool CanCollide(const CollisionObjectWrapper* cow0, const CollisionObjectWrapper* cow1) {
  //   return m_allowedCollisionMatrix(cow0->m_index, cow1->m_index);
  // }
  // void SetLinkIndices();
  // void UpdateAllowedCollisionMatrix();
  // void CheckShapeCast(btCollisionShape* shape, const btTransform& tf0, const btTransform& tf1,
  //     CollisionObjectWrapper* cow, btCollisionWorld* world, vector<Collision>& collisions);


};


namespace trajopt {



CollisionCheckerPtr CreateCollisionChecker(OR::EnvironmentBaseConstPtr env) {
  CollisionCheckerPtr checker(new GrowthCollisionChecker(env));
  return checker;
}
}
