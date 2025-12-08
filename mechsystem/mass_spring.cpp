#include "mass_spring.hpp"
#include "Newmark.hpp"

int main()
{
  MassSpringSystem<2> mss;
  mss.setGravity( {0,-9.81} );
  auto fA = mss.addFix( { { 0.0, 0.0 } } );
  auto mA = mss.addMass( { 1, { 1.0, 0.0 } } );
  mss.addSpring ( { 1, 10, { fA, mA } }  );

  auto mB = mss.addMass( { 1, { 2.0, 0.0 } } );
  mss.addSpring ( { 1, 20, { mA, mB } } );
  mss.addConstraint( { 1.0, { mA, mB } } ); //TODO

  std::cout << "mss: " << std::endl << mss << std::endl;


  double tend = 10;
  double steps = 1000;
  size_t nc = mss.constraints().size();
  size_t nm = mss.masses().size();


  Vector<> x(2*mss.masses().size()+nc);
  Vector<> dx(2*mss.masses().size()+nc);
  Vector<> ddx(2*mss.masses().size()+nc);
  x = 0.0;
  dx = 0.0;
  ddx = 0.0;

  mss.getState(x.range(0, 2*nm), dx.range(0, 2*nm), ddx.range(0, 2*nm));
  //mss.getState (x, dx, ddx);
  auto mss_func = std::make_shared<MSS_Function<2>> (mss);

  //auto mass = std::make_shared<IdentityFunction> (x.size());
  auto mass = std::make_shared<ConstrainedMassMatrixFunction<2>>(mss);

  
  
  SolveODE_Newmark(tend, steps, x, dx,  mss_func, mass,
                   [](double t, VectorView<double> x) { std::cout << "t = " << t
                                                             << ", x = " << Vec<4>(x) << std::endl; });
}
