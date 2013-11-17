#include "trajopt/collision_terms.hpp"
#include "trajopt/collision_checker.hpp"
#include "trajopt/rave_utils.hpp"
#include "trajopt/utils.hpp"
#include "sco/expr_vec_ops.hpp"
#include "sco/expr_ops.hpp"
#include "sco/sco_common.hpp"
#include <boost/foreach.hpp>
#include "utils/eigen_conversions.hpp"
#include "sco/modeling_utils.hpp"
#include "utils/stl_to_string.hpp"
#include "utils/logging.hpp"
#include <boost/functional/hash.hpp>

#include <btBulletCollisionCommon.h>
#include <BulletCollision/CollisionShapes/btShapeHull.h>
#include <BulletCollision/CollisionDispatch/btConvexConvexAlgorithm.h>
btVector3 toBt(const OpenRAVE::Vector& v){
  return btVector3(v[0], v[1], v[2]);
}


using namespace OpenRAVE;
using namespace sco;
using namespace util;
using namespace std;

#define BOOST_FOREACH_PAIR(KEY, VALUE, COL) BOOST_FOREACH_PREAMBLE() \
if (boost::foreach_detail_::auto_any_t BOOST_FOREACH_ID(_foreach_col) \
= BOOST_FOREACH_CONTAIN(COL)) {} else if \
(boost::foreach_detail_::auto_any_t BOOST_FOREACH_ID(_foreach_cur) = \
BOOST_FOREACH_BEGIN(COL)) {} else if \
(boost::foreach_detail_::auto_any_t BOOST_FOREACH_ID(_foreach_end) = \
BOOST_FOREACH_END(COL)) {} else for (bool \
BOOST_FOREACH_ID(_foreach_continue) = true, \
BOOST_FOREACH_ID(_foreach_key_loop) = true; \
BOOST_FOREACH_ID(_foreach_continue) && !BOOST_FOREACH_DONE(COL); \
BOOST_FOREACH_ID(_foreach_continue) ? BOOST_FOREACH_NEXT(COL) : \
(void)0) if (boost::foreach_detail_::set_false(BOOST_FOREACH_ID(_foreach_continue))) \
{} else if (boost::foreach_detail_::set_false(BOOST_FOREACH_ID(_foreach_key_loop))) \
{} else for (KEY = BOOST_FOREACH_DEREF(COL).first; \
!BOOST_FOREACH_ID(_foreach_key_loop); \
BOOST_FOREACH_ID(_foreach_key_loop) = true) for (VALUE = \
BOOST_FOREACH_DEREF(COL).second; !BOOST_FOREACH_ID(_foreach_continue); \
BOOST_FOREACH_ID(_foreach_continue) = true) 

namespace trajopt {


void CollisionsToDistances(const vector<Collision>& collisions, const Link2Int& m_link2ind,
    DblVec& dists) {
  // Note: this checking (that the links are in the list we care about) is probably unnecessary
  // since we're using LinksVsAll
  dists.clear();
  dists.reserve(collisions.size());
  BOOST_FOREACH(const Collision& col, collisions) {
    Link2Int::const_iterator itA = m_link2ind.find(col.linkA);
    Link2Int::const_iterator itB = m_link2ind.find(col.linkB);
    if (itA != m_link2ind.end() || itB != m_link2ind.end()) {
      dists.push_back(col.distance);
    }
  }
}

void CollisionsToDistanceExpressions(const vector<Collision>& collisions, Configuration& rad,
    const Link2Int& link2ind, const VarVector& vars, const DblVec& dofvals, vector<AffExpr>& exprs, bool isTimestep1) {

  exprs.clear();
  exprs.reserve(collisions.size());
  rad.SetDOFValues(dofvals); // since we'll be calculating jacobians
  BOOST_FOREACH(const Collision& col, collisions) {
    AffExpr dist(col.distance);
    Link2Int::const_iterator itA = link2ind.find(col.linkA);
    if (itA != link2ind.end()) {
      VectorXd dist_grad = toVector3d(col.normalB2A).transpose()*rad.PositionJacobian(itA->second, col.ptA);
      exprInc(dist, varDot(dist_grad, vars));
      exprInc(dist, -dist_grad.dot(toVectorXd(dofvals)));
    }
    Link2Int::const_iterator itB = link2ind.find(col.linkB);
    if (itB != link2ind.end()) {
      VectorXd dist_grad = -toVector3d(col.normalB2A).transpose()*rad.PositionJacobian(itB->second, (isTimestep1 && (col.cctype == CCType_Between)) ? col.ptB1 : col.ptB);
      exprInc(dist, varDot(dist_grad, vars));
      exprInc(dist, -dist_grad.dot(toVectorXd(dofvals)));
    }
    if (itA != link2ind.end() || itB != link2ind.end()) {
      exprs.push_back(dist);
    }
  }
  LOG_DEBUG("%ld distance expressions\n", exprs.size());
}

void CollisionsToDistanceExpressions(const vector<Collision>& collisions, Configuration& rad, const Link2Int& link2ind,
    const VarVector& vars0, const VarVector& vars1, const DblVec& vals0, const DblVec& vals1,
    vector<AffExpr>& exprs) {
  vector<AffExpr> exprs0, exprs1;
  CollisionsToDistanceExpressions(collisions, rad, link2ind, vars0, vals0, exprs0, false);
  CollisionsToDistanceExpressions(collisions, rad, link2ind, vars1, vals1, exprs1,true);

  exprs.resize(exprs0.size());
  for (int i=0; i < exprs0.size(); ++i) {
    exprScale(exprs0[i], (1-collisions[i].time));
    exprScale(exprs1[i], collisions[i].time);
    exprs[i] = AffExpr(0);
    exprInc(exprs[i], exprs0[i]);
    exprInc(exprs[i], exprs1[i]);
    cleanupAff(exprs[i]);
  }
}

inline size_t hash(const DblVec& x) {
  return boost::hash_range(x.begin(), x.end());
}

void CollisionEvaluator::GetCollisionsCached(const DblVec& x, vector<Collision>& collisions) {
  double key = hash(getDblVec(x, GetVars()));
  vector<Collision>* it = m_cache.get(key);
  if (it != NULL) {
    LOG_DEBUG("using cached collision check\n");
    collisions = *it;
  }
  else {
    LOG_DEBUG("not using cached collision check\n");
    CalcCollisions(x, collisions);
    m_cache.put(key, collisions);
  }
}

SingleTimestepCollisionEvaluator::SingleTimestepCollisionEvaluator(ConfigurationPtr rad, const VarVector& vars) :
  m_env(rad->GetEnv()),
  m_cc(CollisionChecker::GetOrCreate(*m_env)),
  m_rad(rad),
  m_vars(vars),
  m_link2ind(),
  m_links(),
  m_filterMask(-1) {
  vector<KinBody::LinkPtr> links;
  vector<int> inds;
  rad->GetAffectedLinks(m_links, true, inds);
  for (int i=0; i < m_links.size(); ++i) {
    m_link2ind[m_links[i].get()] = inds[i];
  }
  // TODO add argument
}


void SingleTimestepCollisionEvaluator::CalcCollisions(const DblVec& x, vector<Collision>& collisions) {
  DblVec dofvals = getDblVec(x, m_vars);
  m_rad->SetDOFValues(dofvals);
  m_cc->LinksVsAll(m_links, collisions, m_filterMask);
}

void SingleTimestepCollisionEvaluator::CalcDists(const DblVec& x, DblVec& dists) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  CollisionsToDistances(collisions, m_link2ind, dists);
}


void SingleTimestepCollisionEvaluator::CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  DblVec dofvals = getDblVec(x, m_vars);
  CollisionsToDistanceExpressions(collisions, *m_rad, m_link2ind, m_vars, dofvals, exprs, false);
}


/////////////////////////////////////////////////

  DPCollisionEvaluator::DPCollisionEvaluator(ConfigurationPtr rad, const VarArray& vars, bool continuous) :
  m_env(rad->GetEnv()),
  m_cc(CollisionChecker::GetOrCreate(*m_env)),
  m_rad(rad),
  m_vars(vars),
  m_link2ind(),
  m_links(),
  m_filterMask(-1) {
  // assume continuous is false
    vector<KinBody::LinkPtr> links;
    vector<int> inds;
    rad->GetAffectedLinks(m_links, true, inds);
    for (int i=0; i < m_links.size(); ++i) {
      m_link2ind[m_links[i].get()] = inds[i];
    }
  // TODO add argument
  }


void DPCollisionEvaluator::CalcCollisions(const DblVec& x, vector<Collision>& collisions) {
  // ------------------
  //
  // TODO: get collisions for each timestep
  // then keep only a consistent collision set
  // put those in a vector, and output that to the collisions argument
  //
  // ------------------
  typedef std::map<int /*t*/, vector<Collision> > subsubmap;
  // typedef std::map<const OR::KinBody::Link* /*b*/, subsubmap > submap;
  typedef std::map<const btCollisionShape* /*b*/, subsubmap > sub2map;
  typedef std::map<const btCollisionShape* /*a*/, sub2map > submap;
  typedef std::map<const OR::KinBody::Link* /*a*/, submap >   mainmap;
  mainmap main_collision_map;


  bool continuous = true;

  if (continuous) { //cont-time checking
    for (int t=1; t<m_vars.rows(); t++) {
      VarVector row = m_vars.row(t);
      VarVector row_past = m_vars.row(t-1);
      DblVec dofvals = getDblVec(x, row);
      DblVec dofvals_past = getDblVec(x, row_past);
      m_rad->SetDOFValues(dofvals_past);
      vector<Collision> t_collisions;// = new vector<Collision>();
      m_cc->CastVsAll(*m_rad, m_links, dofvals_past, dofvals, t_collisions);

      BOOST_FOREACH(Collision& collision, t_collisions) {
        main_collision_map[collision.linkB][collision.shapeB][collision.shapeA][t].push_back(collision);
      }
    }
  }


  for (int i=0; i<m_links.size(); i++) {
    submap collision_map;
    if (!continuous) { // discrete checking
      for (int t=0; t<m_vars.rows(); t++) {
        VarVector row = m_vars.row(t);
        DblVec dofvals = getDblVec(x, row);
        m_rad->SetDOFValues(dofvals);
        vector<Collision> t_collisions;// = new vector<Collision>();
        m_cc->LinkVsAll(*m_links[i], t_collisions, m_filterMask);

        // cout << "found " <<t_collisions->size() << " collisions" << endl;
        BOOST_FOREACH(Collision& collision, t_collisions) {
          collision_map[collision.shapeB][collision.shapeA][t].push_back(collision);
        }
      }
    } else {
      collision_map = main_collision_map[m_links[i].get()];
      // cout << "size is" << collision_map.size() << endl;
    }

    BOOST_FOREACH_PAIR(const btCollisionShape* sa,  sub2map val1, collision_map) {
      BOOST_FOREACH_PAIR(const btCollisionShape* b,  subsubmap val2, val1) {
        int time = m_vars.rows();
        vector<double> costs[time];
        vector<Collision> cols[time];
        vector<int> back_pointers[time];

        costs[0].push_back(0);
        cols[0].push_back(Collision(-1.0));
        back_pointers[0].push_back(0);

        BOOST_FOREACH_PAIR(int t, vector<Collision> t_collisions, val2) {
          if (t == 0) continue;
          BOOST_FOREACH(Collision& collision, t_collisions) {
            cols[t].push_back(collision);
            double cost = collision.distance;
            cost = max(0.0, 0.0 - cost);
            double min_cost = INFINITY;
            int back_ptr = 0;
            for(int j=0; j < cols[t-1].size(); j++) {
              if (toBt(cols[t-1][j].normalB2A * (1.0 / cols[t-1][j].distance)).angle(toBt(collision.normalB2A * (1.0 / collision.distance))) > 2.0) {
                continue;
              }
              // OR::Vector diff = cols[t-1][j].normalB2A * (1.0 / fabs(cols[t-1][j].distance)) - collision.normalB2A * (1.0 / fabs(collision.distance));
              // if (fabs(diff.x) < 0.01 && fabs(diff.y) < 0.01 && fabs(diff.z) < 0.01) {
              //   continue;
              // }
              double path_cost = costs[t-1][j] + cost;
              if (path_cost < min_cost) {
                min_cost = path_cost;
                back_ptr = j;
              }
            }
            costs[t].push_back(min_cost);
            back_pointers[t].push_back(back_ptr);
          }
        }

        double min_cost = INFINITY;
        int ptr = -1;
        int current_time = time-1;
        for(int j=0; j < costs[current_time].size(); j++) {
          if (costs[current_time][j] < min_cost) {
            min_cost = costs[current_time][j];
            ptr = j;
          }
        }
        // cout << "Min cost is" << min_cost << endl;
        assert (ptr != -1);
        while (current_time != 0) {
          // cout << "Adding collision" << cols[current_time][ptr].distance << endl;
          // if (cols[current_time][ptr].distance < m_cc->GetContactDistance()) {
          if (cols[current_time][ptr].distance < 0) {
            Collision c = cols[current_time][ptr];
            c.weight = current_time;
            collisions.push_back(c);
            // cout << "collision between " << c.linkA->GetName() << " and " << c.linkB->GetName() << endl;
          }
          ptr = back_pointers[current_time][ptr];
          current_time--;
        }
      }
    }
  }



  // for (int i=0; i < m_vars.rows(); i++) {
  //   VarVector row = m_vars.row(i);
  //   DblVec dofvals = getDblVec(x, row);
  //   m_rad->SetDOFValues(dofvals);
  //   vector<Collision>* t_collisions = new vector<Collision>();
  //   m_cc->LinksVsAll(m_links, *t_collisions, m_filterMask);
  //   cout << "Number of collisions is" << t_collisions->size() << endl;

  //   // collisions.reserve(collisions.size() + distance(t_collisions->begin(),t_collisions->end()));
  //   // collisions.insert(collisions.end(),t_collisions->begin(),t_collisions->end());
  // }


  // Run DP


}

void DPCollisionEvaluator::CalcDists(const DblVec& x, DblVec& dists) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  CollisionsToDistances(collisions, m_link2ind, dists);
}


void DPCollisionEvaluator::CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  map<int, vector<Collision> > timed_collisions;
  BOOST_FOREACH(Collision c, collisions) {
    timed_collisions[c.weight].push_back(c);
  }
  for (int i=0; i < m_vars.rows(); i++) {
    VarVector row = m_vars.row(i);
    DblVec dofvals = getDblVec(x, row);
    vector<AffExpr> exprs_;
    CollisionsToDistanceExpressions(timed_collisions[i], *m_rad, m_link2ind, row, dofvals, exprs_, false);
    exprs.insert(exprs.end(), exprs_.begin(), exprs_.end());
  }
}


////////////////////////////////////////

CastCollisionEvaluator::CastCollisionEvaluator(ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1) :
  m_env(rad->GetEnv()),
  m_cc(CollisionChecker::GetOrCreate(*m_env)),
  m_rad(rad),
  m_vars0(vars0),
  m_vars1(vars1),
  m_link2ind(),
  m_links() {
  vector<KinBody::LinkPtr> links;
  vector<int> inds;
  rad->GetAffectedLinks(m_links, true, inds);
  for (int i=0; i < m_links.size(); ++i) {
    m_link2ind[m_links[i].get()] = inds[i];
  }
}

void CastCollisionEvaluator::CalcCollisions(const DblVec& x, vector<Collision>& collisions) {
  DblVec dofvals0 = getDblVec(x, m_vars0);
  DblVec dofvals1 = getDblVec(x, m_vars1);
  m_rad->SetDOFValues(dofvals0);
  m_cc->CastVsAll(*m_rad, m_links, dofvals0, dofvals1, collisions);
}
void CastCollisionEvaluator::CalcDistExpressions(const DblVec& x, vector<AffExpr>& exprs) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  DblVec dofvals0 = getDblVec(x, m_vars0);
  DblVec dofvals1 = getDblVec(x, m_vars1);
  CollisionsToDistanceExpressions(collisions, *m_rad, m_link2ind, m_vars0, m_vars1, dofvals0, dofvals1, exprs);
}
void CastCollisionEvaluator::CalcDists(const DblVec& x, DblVec& dists) {
  vector<Collision> collisions;
  GetCollisionsCached(x, collisions);
  CollisionsToDistances(collisions, m_link2ind, dists);
}


//////////////////////////////////////////



typedef OpenRAVE::RaveVector<float> RaveVectorf;

void PlotCollisions(const std::vector<Collision>& collisions, OR::EnvironmentBase& env, vector<OR::GraphHandlePtr>& handles, double safe_dist) {
  BOOST_FOREACH(const Collision& col, collisions) {
    RaveVectorf color;
    if (col.distance < 0) color = RaveVectorf(1,0,0,1);
    else if (col.distance < safe_dist) color = RaveVectorf(1,1,0,1);
    else color = RaveVectorf(0,1,0,1);
    if (col.cctype == CCType_Between) {
      handles.push_back(env.drawarrow(col.ptB, col.ptB1, .002, RaveVectorf(0,0,0,1)));
    }
    OR::Vector ptB = (col.cctype == CCType_Between)  ? ((1-col.time)* col.ptB +col.time*col.ptB1) : col.ptB;
    handles.push_back(env.drawarrow(col.ptA, ptB, .0025, color));
  }
}

CollisionCost::CollisionCost(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars) :
    Cost("collision"),
    m_calc(new SingleTimestepCollisionEvaluator(rad, vars)), m_dist_pen(dist_pen), m_coeff(coeff)
{}

CollisionCost::CollisionCost(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1) :
    Cost("cast_collision"),
    m_calc(new CastCollisionEvaluator(rad, vars0, vars1)), m_dist_pen(dist_pen), m_coeff(coeff)
{}
ConvexObjectivePtr CollisionCost::convex(const vector<double>& x, Model* model) {
  ConvexObjectivePtr out(new ConvexObjective(model));
  vector<AffExpr> exprs;
  m_calc->CalcDistExpressions(x, exprs);
  for (int i=0; i < exprs.size(); ++i) {
    AffExpr viol = exprSub(AffExpr(m_dist_pen), exprs[i]);
    out->addHinge(viol, m_coeff);
  }
  return out;
}
double CollisionCost::value(const vector<double>& x) {
  DblVec dists;
  m_calc->CalcDists(x, dists);
  double out = 0;
  for (int i=0; i < dists.size(); ++i) {
    out += pospart(m_dist_pen - dists[i]) * m_coeff;
  }
  return out;
}

void CollisionCost::Plot(const DblVec& x, OR::EnvironmentBase& env, std::vector<OR::GraphHandlePtr>& handles) {
  vector<Collision> collisions;
  m_calc->GetCollisionsCached(x, collisions);
  PlotCollisions(collisions, env, handles, m_dist_pen);
}

// DP version


  DPCollisionCost::DPCollisionCost(double dist_pen, DblVec coeffs, ConfigurationPtr rad, const VarArray& vars, bool continuous) :
  Cost("DP collision"),
  m_calc(new DPCollisionEvaluator(rad, vars, continuous)), m_dist_pen(dist_pen), m_coeff(coeffs[0])
{}

ConvexObjectivePtr DPCollisionCost::convex(const vector<double>& x, Model* model) {
  ConvexObjectivePtr out(new ConvexObjective(model));
  vector<AffExpr> exprs;
  m_calc->CalcDistExpressions(x, exprs);
  for (int i=0; i < exprs.size(); ++i) {
    AffExpr viol = exprSub(AffExpr(m_dist_pen), exprs[i]);
    out->addHinge(viol, m_coeff);
  }
  return out;
}
double DPCollisionCost::value(const vector<double>& x) {
  DblVec dists;
  m_calc->CalcDists(x, dists);
  double out = 0;
  for (int i=0; i < dists.size(); ++i) {
    out += pospart(m_dist_pen - dists[i]) * m_coeff;
  }
  return out;
}

void DPCollisionCost::Plot(const DblVec& x, OR::EnvironmentBase& env, std::vector<OR::GraphHandlePtr>& handles) {
  vector<Collision> collisions;
  m_calc->GetCollisionsCached(x, collisions);
  PlotCollisions(collisions, env, handles, m_dist_pen);
}

// ALMOST EXACTLY COPIED FROM CollisionCost

CollisionConstraint::CollisionConstraint(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars) :
    m_calc(new SingleTimestepCollisionEvaluator(rad, vars)), m_dist_pen(dist_pen), m_coeff(coeff)
{
  name_="collision";
}

CollisionConstraint::CollisionConstraint(double dist_pen, double coeff, ConfigurationPtr rad, const VarVector& vars0, const VarVector& vars1) :
    m_calc(new CastCollisionEvaluator(rad, vars0, vars1)), m_dist_pen(dist_pen), m_coeff(coeff)
{
  name_="collision";
}
ConvexConstraintsPtr CollisionConstraint::convex(const vector<double>& x, Model* model) {
  ConvexConstraintsPtr out(new ConvexConstraints(model));
  vector<AffExpr> exprs;
  m_calc->CalcDistExpressions(x, exprs);
  for (int i=0; i < exprs.size(); ++i) {
    AffExpr viol = exprSub(AffExpr(m_dist_pen), exprs[i]);
    out->addIneqCnt(exprMult(viol,m_coeff));
  }
  return out;
}
DblVec CollisionConstraint::value(const vector<double>& x) {
  DblVec dists;
  m_calc->CalcDists(x, dists);
  DblVec out(dists.size());
  for (int i=0; i < dists.size(); ++i) {
    out[i] = pospart(m_dist_pen - dists[i]) * m_coeff;
  }
  return out;
}



}
