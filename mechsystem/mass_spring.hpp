#ifndef MASS_SPRING_HPP
#define MASS_SPRING_HPP

#include <nonlinfunc.hpp>
#include <timestepper.hpp>

using namespace ASC_ode;

#include <vector.hpp>
using namespace nanoblas;


template <int D>
class Mass
{
public:
  double mass;
  Vec<D> pos;
  Vec<D> vel = 0.0;
  Vec<D> acc = 0.0;
};


template <int D>
class Fix
{
public:
  Vec<D> pos;
};


class Connector
{
public:
  enum CONTYPE { FIX=1, MASS=2 };
  CONTYPE type;
  size_t nr;
};

std::ostream & operator<< (std::ostream & ost, const Connector & con)
{
  ost << "type = " << int(con.type) << ", nr = " << con.nr;
  return ost;
}

class Spring
{
public:
  double length;  
  double stiffness;
  std::array<Connector,2> connectors;
};

class Chain
{
public:
  double length;
  double stiffness; // we approximate the chain by a very stiff spring when the constraints are enforced
  std::array<Connector,2> connectors;
};

class DistanceConstraint
{
public:
    double length;
    std::array<Connector,2> connectors;
};

template <int D>
class MassSpringSystem
{
  std::vector<Fix<D>> m_fixes;
  std::vector<Mass<D>> m_masses;
  std::vector<Spring> m_springs;
  std::vector<Chain> m_chains;
  std::vector<DistanceConstraint> m_constraints;

  Vec<D> m_gravity=0.0;
public:
  void setGravity (Vec<D> gravity) { m_gravity = gravity; }
  Vec<D> getGravity() const { return m_gravity; }

  Connector addFix (Fix<D> p)
  {
    m_fixes.push_back(p);
    return { Connector::FIX, m_fixes.size()-1 };
  }

  Connector addMass (Mass<D> m)
  {
    m_masses.push_back (m);
    return { Connector::MASS, m_masses.size()-1 };
  }
  
  size_t addSpring (Spring s) 
  {
    m_springs.push_back (s); 
    return m_springs.size()-1;
  }

  size_t addChain (Chain c)
  {
    m_chains.push_back (c);
    return m_chains.size()-1;
  }

  

  size_t addConstraint(DistanceConstraint c) {
      m_constraints.push_back(c);
      return m_constraints.size()-1;
  }

  


  auto & fixes() { return m_fixes; } 
  auto & masses() { return m_masses; } 
  auto & springs() { return m_springs; }
  auto & chains() { return m_chains; }
  auto& constraints() { return m_constraints; }

  void getState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {

    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        valmat.row(i) = m_masses[i].pos;
        dvalmat.row(i) = m_masses[i].vel;
        ddvalmat.row(i) = m_masses[i].acc;
      }
  }

  void setState (VectorView<> values, VectorView<> dvalues, VectorView<> ddvalues)
  {
    auto valmat = values.asMatrix(m_masses.size(), D);
    auto dvalmat = dvalues.asMatrix(m_masses.size(), D);
    auto ddvalmat = ddvalues.asMatrix(m_masses.size(), D);

    for (size_t i = 0; i < m_masses.size(); i++)
      {
        m_masses[i].pos = valmat.row(i);
        m_masses[i].vel = dvalmat.row(i);
        m_masses[i].acc = ddvalmat.row(i);
      }
  }
};

template <int D>
std::ostream & operator<< (std::ostream & ost, MassSpringSystem<D> & mss)
{
  ost << "fixes:" << std::endl;
  for (auto f : mss.fixes())
    ost << f.pos << std::endl;

  ost << "masses: " << std::endl;
  for (auto m : mss.masses())
    ost << "m = " << m.mass << ", pos = " << m.pos << std::endl;

  ost << "springs: " << std::endl;
  for (auto sp : mss.springs())
    ost << "length = " << sp.length << ", stiffness = " << sp.stiffness
        << ", C1 = " << sp.connectors[0] << ", C2 = " << sp.connectors[1] << std::endl;

  ost << "chains: " << std::endl;
  for (auto ch : mss.chains())
    ost << "length = " << ch.length << ", stiffness = " << ch.stiffness
        << ", C1 = " << ch.connectors[0] << ", C2 = " << ch.connectors[1] << std::endl;

  ost << "Constraints:\n";
  for (auto c : mss.constraints())
      ost << "L=" << c.length
          << " C1=" << c.connectors[0]
          << " C2=" << c.connectors[1] << "\n";

  return ost;
}

template <int D>
class MSS_Function : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  MSS_Function (MassSpringSystem<D> & _mss)
    : mss(_mss) { }

    // unkowns: positions + constraints
  virtual size_t dimX() const override { return D*mss.masses().size()+ mss.constraints().size(); }
  virtual size_t dimF() const override{ return D*mss.masses().size()+ mss.constraints().size(); }

  virtual void evaluate (VectorView<double> x, VectorView<double> f) const override
  {
    f = 0.0;

    size_t nm = mss.masses().size();
    size_t nc = mss.constraints().size();

    // split x into coordinates and multipliers
    auto xmat = x.range(0, nm*D).asMatrix(nm, D);
    auto fmat = f.range(0, nm*D).asMatrix(nm, D);
    auto lambda = x.range(nm*D, nm*D + nc);

    // gravity
    for (size_t i = 0; i < nm; i++)
      fmat.row(i) = mss.masses()[i].mass*mss.getGravity();

    for (auto spring : mss.springs())
      {
        auto [c1,c2] = spring.connectors;
        Vec<D> p1, p2; // positions of the two endpoints of the spring
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        double force = spring.stiffness * (norm(p1-p2)-spring.length);
        Vec<D> dir12 = 1.0/norm(p1-p2) * (p2-p1);

        if (c1.type == Connector::MASS)
          fmat.row(c1.nr) += force*dir12;
        if (c2.type == Connector::MASS)
          fmat.row(c2.nr) -= force*dir12;
      }
    
    for(auto chain : mss.chains()){
        auto [c1,c2] = chain.connectors;
        Vec<D> p1, p2; // positions of the two endpoints of the chain
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);
        
        double force = 0.0;
        if(norm(p1-p2) > chain.length){
            force = chain.stiffness * (norm(p1-p2)-chain.length);
            Vec<D> dir12 = 1.0/norm(p1-p2) * (p2-p1);
            if (c1.type == Connector::MASS)
              fmat.row(c1.nr) += force*dir12;
            if (c2.type == Connector::MASS)
              fmat.row(c2.nr) -= force*dir12;
        }
    }

    // constraints: add G^T * lambda
    for (size_t k = 0; k < nc; k++)
      {
        auto & C = mss.constraints()[k];
        auto [c1,c2] = C.connectors;

        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        double dist = norm(p1-p2);
        Vec<D> dC1 = 1.0/norm(p1-p2) * (p2-p1);
        Vec<D> dC2 = -dC1;

        double lam = lambda(k);

        if (c1.type == Connector::MASS)
          fmat.row(c1.nr) += lam * dC1;
        if (c2.type == Connector::MASS)
          fmat.row(c2.nr) += lam * dC2;
      }

    // constraint equations C(x)=0
    for (size_t k = 0; k < nc; k++)
      {
        auto & C = mss.constraints()[k];
        auto [c1,c2] = C.connectors;

        Vec<D> p1, p2;
        if (c1.type == Connector::FIX)
          p1 = mss.fixes()[c1.nr].pos;
        else
          p1 = xmat.row(c1.nr);
        if (c2.type == Connector::FIX)
          p2 = mss.fixes()[c2.nr].pos;
        else
          p2 = xmat.row(c2.nr);

        f(nm*D + k) = norm(p1-p2) - C.length;
      }
    // // divide by mass
    // for (size_t i = 0; i < mss.masses().size(); i++)
    //   fmat.row(i) *= 1.0/mss.masses()[i].mass;
  }
  
  virtual void evaluateDeriv (VectorView<double> x, MatrixView<double> df) const override
  {
      df = 0.0;

      size_t nm = mss.masses().size();
      size_t nc = mss.constraints().size();

      auto xmat = x.range(0, nm * D).asMatrix(nm, D);
      auto lambda = x.range(nm * D, nm * D + nc);

      // ----------------------------------------------------------------------
      // 1. SPRING JACOBIAN  (exact tangent stiffness)
      // ----------------------------------------------------------------------
      for (auto &spring : mss.springs())
      {
          auto [c1, c2] = spring.connectors;

          Vec<D> p1 = (c1.type == Connector::FIX)
                        ? mss.fixes()[c1.nr].pos
                        : xmat.row(c1.nr);
          Vec<D> p2 = (c2.type == Connector::FIX)
                        ? mss.fixes()[c2.nr].pos
                        : xmat.row(c2.nr);

          Vec<D> d  = p2 - p1;
          double L  = norm(d);
          if (L < 1e-12) continue;

          Vec<D> n = (1.0 / L) * d;  // FIXED: Use scalar * vector

          double k  = spring.stiffness;
          double L0 = spring.length;

          // K = k*n*n^T + (k*(L-L0)/L) * (I - n*n^T)
          Matrix<double> K(D, D);  // FIXED: Use Matrix<double>
          double force_over_L = k * (L - L0) / L;

          for (size_t i=0;i<D;i++)
              for (size_t j=0;j<D;j++)
              {
                  double ninj = n(i)*n(j);
                  K(i,j) = k*ninj + force_over_L*((i==j?1.0:0.0) - ninj);
              }

          // Insert into global Jacobian: df = dF/dx
          auto addBlock = [&](size_t r, size_t c, double val)
          {
              df(r, c) += val;
          };

          if (c1.type == Connector::MASS)
          {
              for (size_t i=0;i<D;i++)
                  for (size_t j=0;j<D;j++)
                      addBlock(c1.nr*D + i, c1.nr*D + j, -K(i,j));   // df1/dp1
          }

          if (c2.type == Connector::MASS)
          {
              for (size_t i=0;i<D;i++)
                  for (size_t j=0;j<D;j++)
                      addBlock(c2.nr*D + i, c2.nr*D + j, -K(i,j));   // df2/dp2
          }

          if (c1.type == Connector::MASS && c2.type == Connector::MASS)
          {
              for (size_t i=0;i<D;i++)
                  for (size_t j=0;j<D;j++)
                  {
                      addBlock(c1.nr*D + i, c2.nr*D + j,  K(i,j));   // df1/dp2
                      addBlock(c2.nr*D + i, c1.nr*D + j,  K(i,j));   // df2/dp1
                  }
          }
      }

      // ----------------------------------------------------------------------
      // 2. DISTANCE CONSTRAINT JACOBIAN (exact)
      // ----------------------------------------------------------------------
      for (size_t k = 0; k < nc; k++)
      {
          auto & Ck = mss.constraints()[k];
          auto [c1, c2] = Ck.connectors;

          Vec<D> p1 = (c1.type == Connector::FIX)
                        ? mss.fixes()[c1.nr].pos
                        : xmat.row(c1.nr);
          Vec<D> p2 = (c2.type == Connector::FIX)
                        ? mss.fixes()[c2.nr].pos
                        : xmat.row(c2.nr);

          Vec<D> d = p2 - p1;
          double L = norm(d);
          if (L < 1e-12) continue;

          Vec<D> n = (1.0 / L) * d;  // FIXED: Use scalar * vector
          double lam = lambda(k);

          // Matrix M = (I - n·n^T)/L
          Matrix<double> M(D, D);  // FIXED: Use Matrix<double>
          for (size_t i=0;i<D;i++)
              for (size_t j=0;j<D;j++)
                  M(i,j) = ((i==j)?1.0:0.0) - n(i)*n(j);
          M *= (1.0 / L);

          auto addBlock = [&](size_t r, size_t c, double val)
          {
              df(r, c) += val;
          };

          // ----------------------------------------------------------
          // 2A. dF/dλ: add n and -n into force rows
          // ----------------------------------------------------------
          if (c1.type == Connector::MASS)
              for (size_t i=0;i<D;i++)
                  addBlock(c1.nr*D + i, nm*D + k,  n(i));

          if (c2.type == Connector::MASS)
              for (size_t i=0;i<D;i++)
                  addBlock(c2.nr*D + i, nm*D + k, -n(i));

          // ----------------------------------------------------------
          // 2B. dF/dx: tangent stiffness from λ * n
          // ----------------------------------------------------------
          if (c1.type == Connector::MASS)
          {
              for (size_t i=0;i<D;i++)
                  for (size_t j=0;j<D;j++)
                      addBlock(c1.nr*D + i, c1.nr*D + j, -lam * M(i,j)); // df1/dp1

              if (c2.type == Connector::MASS)
                  for (size_t i=0;i<D;i++)
                      for (size_t j=0;j<D;j++)
                          addBlock(c1.nr*D + i, c2.nr*D + j, +lam * M(i,j)); // df1/dp2
          }

          if (c2.type == Connector::MASS)
          {
              if (c1.type == Connector::MASS)
                  for (size_t i=0;i<D;i++)
                      for (size_t j=0;j<D;j++)
                          addBlock(c2.nr*D + i, c1.nr*D + j, +lam * M(i,j)); // df2/dp1

              for (size_t i=0;i<D;i++)
                  for (size_t j=0;j<D;j++)
                      addBlock(c2.nr*D + i, c2.nr*D + j, -lam * M(i,j)); // df2/dp2
          }

          // ----------------------------------------------------------
          // 2C. dC/dx: bottom constraint row
          // ----------------------------------------------------------
          // FIXED: Signs were backwards!
          if (c1.type == Connector::MASS)
              for (size_t j=0;j<D;j++)
                  addBlock(nm*D + k, c1.nr*D + j, -n(j));   // dC/dp1 = -n

          if (c2.type == Connector::MASS)
              for (size_t j=0;j<D;j++)
                  addBlock(nm*D + k, c2.nr*D + j, +n(j));   // dC/dp2 = +n
      }
  }
  
};

template <int D>
class ConstrainedMassMatrixFunction : public NonlinearFunction
{
  MassSpringSystem<D> & mss;
public:
  ConstrainedMassMatrixFunction(MassSpringSystem<D> & _mss) : mss(_mss) {}

  virtual size_t dimX() const override {
    return D * mss.masses().size() + mss.constraints().size();
  }

  virtual size_t dimF() const override {
    return D * mss.masses().size() + mss.constraints().size();
  }

  virtual void evaluate(VectorView<double> x, VectorView<double> f) const override {
    size_t nm = mss.masses().size();
    size_t nc = mss.constraints().size();

    // Mass part
    for (size_t i = 0; i < nm; i++)
      for (size_t d = 0; d < D; d++)
        f(D * i + d) = mss.masses()[i].mass * x(D * i + d);

    // Constraint part
    for (size_t k = 0; k < nc; k++)
      f(nm * D + k) = 0.0;
  }

  virtual void evaluateDeriv(VectorView<double> x, MatrixView<double> df) const override {
    size_t nm = mss.masses().size();
    size_t nc = mss.constraints().size();
    size_t N = D * nm + nc;

    df = 0.0;
    // Diagonal for mass part
    for (size_t i = 0; i < nm; i++)
      for (size_t d = 0; d < D; d++)
        df(D * i + d, D * i + d) = mss.masses()[i].mass;
    // Constraints part remains zero
  }
};

#endif
