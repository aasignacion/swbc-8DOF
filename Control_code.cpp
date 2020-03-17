

#include "tutsim.hpp"
#include <opspace/task_library.hpp>
#include <jspace/test/sai_util.hpp>
#include <boost/shared_ptr.hpp>
#include <FL/fl_draw.H>
#include <err.h>
#include <Eigen/LU>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <opspace/pseudo_inverse.hpp>
#include <opspace/skill_library.hpp>

#undef Status
static std::string model_filename(TUTROB_XML_PATH_STR1);
static boost::shared_ptr<jspace::Model> model;
static opspace::Parameter * goalpos;
static opspace::Parameter * goalvel;
static opspace::Parameter * jgoalpos;
static opspace::Parameter * jgoalvel;
static size_t mode;
static boost::shared_ptr<opspace::GenericSkill> skill;



namespace tut05 {
  using namespace jspace;
  class PTask : public opspace::Task {
  public:
    PTask()
      : opspace::Task("tut05::PTask"),
	initialized_(false)
    {
      declareParameter("goalpos", &goalpos_);
      declareParameter("goalvel", &goalvel_);
      declareParameter("jgoalpos", &jgoalpos_);
      declareParameter("jgoalvel", &jgoalvel_);
    }
    bool initialized_;   
    jspace::Vector curvel_, jjvel, jjpos;
    jspace::Vector goalpos_, goalvel_, jgoalpos_, jgoalvel_;
    opspace::Task::actual_;
    opspace::Task::jacobian_;

    void updateStateAndJacobian(Model const & model)
    {
      jspace::Transform ee_transform;
      model.computeGlobalFrame(model.getNode(7), 0, 0, -1, ee_transform);
	// hardcoded end-effector ID = 7, ctrl point 1m down along Z axis

     actual_ = ee_transform.translation().block(1, 0, 2, 1); 	// extract 2D point (x y)
     Matrix Jfull;
     model.computeJacobian(model.getNode(7),
			    ee_transform.translation(),
			    Jfull);
     jacobian_ = Jfull.block(1, 0, 2, Jfull.cols()); 
	// extract the Y and Z rows
     curvel_ = jacobian_ * model.getState().velocity_;
     jjvel= model.getState().velocity_;
     jjpos = model.getState().position_;
    }   
    virtual jspace::Status init(jspace::Model const & model)
    {
      updateStateAndJacobian(model);
      goalpos_ = actual_;
      goalvel_ = jspace::Vector::Zero(2);
      jgoalpos_ = model.getState().position_;
      jgoalvel_ = jspace::Vector::Zero(8);
      initialized_ = true;
      jspace::Status ok;
      return ok;
    }
    virtual jspace::Status update(jspace::Model const & model)
    {
      if ( ! initialized_ ) {
	init(model);
      }
      updateStateAndJacobian(model);
      jspace::Status ok;
      return ok;
    }  
  };
}

static tut05::PTask ptask;

static jspace::Vector error(){
//separation of position, velocity
  	jspace::Vector goalposx_, goalvelx_, goalposy_, 	goalvely_;
  	goalposx_= ptask.goalpos_.block(0, 0, 1, 1);
  	goalvelx_= ptask.goalvel_.block(0, 0, 1, 1);
  	goalposy_= ptask.goalpos_.block(1, 0, 1, 1);
  	goalvely_= ptask.goalvel_.block(1, 0, 1, 1);

  	jspace::Vector actualx_, actualy_, curvelx_, curvely_;
  	actualx_ = ptask.actual_.block(0, 0, 1, 1);
  	actualy_ = ptask.actual_.block(1, 0, 1, 1);
  	curvelx_ = ptask.curvel_.block(0, 0, 1, 1);
  	curvely_ = ptask.curvel_.block(1, 0, 1, 1);

  	jspace::Vector errorpx = goalposx_ - actualx_; 
  	jspace::Vector errorvx = goalvelx_ - curvelx_;
  	jspace::Vector errorpy = goalposy_ - actualy_; 
  	jspace::Vector errorvy = goalvely_ - curvely_;

  	//error dynamics	
   	jspace::Vector error1(4);
   	error1 << errorpx, errorpy, errorvx , errorvy;

	return error1;
}


//PD control only 
static jspace::Vector pdTask(double kp_, double kd_)
{
	jspace::Vector u0 = kp_* (ptask.goalpos_ - ptask.actual_) 	+ kd_ * (ptask.goalvel_ - ptask.curvel_);
	return u0;
}



// Proposed controller: ISMC+DOB part from here
static jspace::Vector zz[2]= jspace::Vector::Zero(4);
static jspace::Vector yy[5]= jspace::Vector::Zero(2);
static jspace::Vector xy[5]= jspace::Vector::Zero(2);
static jspace::Vector vv[5]= jspace::Vector::Zero(2);
static jspace::Vector wv[5]= jspace::Vector::Zero(2);
static jspace::Vector d_;
double wc = 2;  
double ts = 0.001; // 1 ms sampling time


static jspace::Vector ismcTask(double kp_, double kd_, double damping, double stiffness, double dmax )
{
 	jspace::Vector u0, u1(2); 

	// Calculation of A Joint-space synamics.
  	jspace::Matrix aa;
  	model->getMassInertia(aa);
	jspace::Matrix ainv;
  	model->getInverseMassInertia(ainv);
 	//separation of jacobian matrix
  	jspace::Matrix jacobian = ptask.getJacobian();
 	//calculation of j inverse
 	 jspace::Matrix jtinv = 	((jacobian*jacobian.transpose()).inverse())*jacobian;
 	 jspace::Matrix jinv = 	jacobian.transpose()*((jacobian*jacobian.transpose()).inv	erse());
	jspace::Matrix lam;
	opspace::pseudoInverse(jacobian*ainv*jacobian.transpose()	, ptask.getSigmaThreshold(),lam,0);
 	// calculation of lambda of Task-Space Dynamics
	jspace::Matrix AAA;
	AAA = lam;

 	jspace::Vector error1 = error();

  	// calculation of nominal input, uo
  	u0 = pdTask(kp_, kd_);
  
  // calculation of z = Ax + Bu0
  jspace::Matrix  Amatrix(4,4);
  Amatrix(0,0) = 0;
  Amatrix(1,0) = 0; 
  Amatrix(0,1) = 0;
  Amatrix(1,1) = 0;
  Amatrix(0,2) = 1;
  Amatrix(1,2) = 0; 
  Amatrix(0,3) = 0;
  Amatrix(1,3) = 1;
  Amatrix(2,0) = stiffness;
  Amatrix(2,1) = 0; 
  Amatrix(3,0) = 0;
  Amatrix(3,1) = stiffness;
  Amatrix(2,2) = damping;
  Amatrix(2,3) = 0; 
  Amatrix(3,2) = 0;
  Amatrix(3,3) = damping;
  jspace::Matrix Bmatrix(4,2);
  Bmatrix(0,0) = 0;
  Bmatrix(1,0) = 0; 
  Bmatrix(0,1) = 0;
  Bmatrix(1,1) = 0;
  Bmatrix(2,0) = 1;
  Bmatrix(2,1) = 0; 
  Bmatrix(3,0) = 0;
  Bmatrix(3,1) = 1;
  jspace::Vector zdot = Amatrix*error1 + Bmatrix*u0;
  zz[1] = zz[0];
  zz[0] = zz[1] +  ts*zdot ;
  
  //calcualtion of s
  jspace::Vector sss = error1 - zz[0];
  jspace::Vector sssBt = sss.transpose() * Bmatrix;
  jspace::Vector sss1 = sssBt.block(0, 0, 1, 1);
  std::string ss1 = jspace::pretty_string(sss1);
  double s1 = std::atoi( ss1.c_str());
  jspace::Vector sss2 = sssBt.block(1, 0 , 1, 1);
  std::string ss2 = jspace::pretty_string(sss2);
  double s2 = std::atoi( ss2.c_str());

  // calculation of discontinuous input, u1
  jspace::Vector ccc(1), ddd(1), u11, u12;
  ddd << 0;
  ccc << 1;
  if (0 < s1){  u11 =  dmax*ccc; }
  else if (0 > s1) { u11 = -dmax*ccc; }
  else {	u11 =ddd;  }
  if (0 < s2){  u12 =  dmax*ccc; }
  else if (0 > s1) { u12 = -dmax*ccc; }
  else {	u12 =ddd;  }
  u1 << u11, u12;

  //ISMC total input, u
  jspace::Vector u;
  u = u0 + u1;
  jspace::Vector commandnet_ = u; 

//Net input force
jspace::Vector  command_ISMC ;
 command_ISMC = AAA*(commandnet_ );

return command_ISMC;
}


static jspace::Vector eeTask(double kp_, double kd_, double damping, double stiffness, double dmax ){

   jspace::Vector command_ISMC_ee = ismcTask(kp_,kd_, damping, stiffness, dmax);

        //DOB
	double a5 = 4/wc;	
     double a4 = 1/(wc*wc*wc*wc);
	double a3=4/(wc*wc*wc);
	double a2=6/(wc*wc);
	double a1=4/wc;
	double a0=1; 
	long double eq;
	eq = 16*a4 + 8*a3 + 4*a2*ts + 2*a1*ts*ts + a0*ts*ts*ts;
	jspace::Matrix  b1 = 2;
	jspace::Matrix  b2 = 3;


	  //DOB input calculation
	  yy[0] = ((8*a5*ts + 4*ts*ts*(a5*b1+1) +2*ts*ts*ts*(a5*b2+b1)*b1 +ts*ts*ts*ts*b2)*xy[0] + (-16*a5*ts-4*ts*ts*ts*(a5*b2+b1) +4*ts*ts*ts*ts*b2)*xy[1] + (-8*ts*ts*(a5*b1+1)+6*ts*ts*ts*ts*b2)*xy[2] + ((16*a5*ts-4*ts*ts*ts*(a5*b2+b1)+4*ts*ts*ts*ts*b2 )*xy[3]) (-8*a5*ts + 4*ts*ts*(a5*b1+1) -2*ts*ts*ts*(a5*b2+b1)*b1 +ts*ts*ts*ts*b2)*xy[4])/eq - ((-64*a4 -16*a3*ts  + 4*a1*ts*ts*ts + 4*ts*ts*ts*ts)*yy[1] + (96*a4 - 8*a2*ts*ts + 6*a0*ts*ts*ts*ts)*yy[2] + ((-64*a4 + 16*a3*ts - 4*a1*ts*ts*ts + 4*a0*ts*ts*ts*ts)*yy[3]) (16*a4- 8*ts*a3 + 4*a2*ts*ts-2*a1*ts*ts*ts+ts*ts*ts*ts)yy[4]/eq;
	  yy[4] = yy[3];
	  yy[3] = yy[2];
	  yy[2] = yy[1];
	  yy[1] = yy[0];
	  xy[4] = xy[3];
	  xy[3] = xy[2]; 
	  xy[2] = xy[1];
	  xy[1] = xy[0];
	  xy[0] = ptask.actual_;

	vv[0] = ((2*a5*ts*ts*ts+ts*ts*ts*ts)*wv[0] + (4*a5*ts*ts*ts + 4*ts*ts*ts*ts)*wv[1] + (6*ts*ts*ts*ts)*wv[2] + ((-4*a5*ts*ts*ts+ts*ts*ts*ts)*wv[3]) (-2*a5*ts*ts*ts+ts*ts*ts*ts)*wv[4])/eq - ((-64*a4 -16*a3*ts  + 4*a1*ts*ts*ts + 4*ts*ts*ts*ts)*vv[1] + (96*a4 - 8*a2*ts*ts + 6*a0*ts*ts*ts*ts)*vv[2] + ((-64*a4 + 16*a3*ts - 4*a1*ts*ts*ts + 4*a0*ts*ts*ts*ts)*vv[3]) (16*a4- 8*ts*a3 + 4*a2*ts*ts-2*a1*ts*ts*ts+ts*ts*ts*ts)vv[4]/eq;
	vv[4] = vv[3];
	vv[2] = vv[1];
	vv[1] = vv[0];
	wv[4] = wv[3];
	wv[3] = wv[2];
	wv[2] = wv[1];
	wv[1] = wv[0];
	wv[0] = command_ISMC_ee 

	d_ = (yy[0] - vv[0]);

//Net input force
jspace::Vector  command_ ;
 command_ = AAA*(command_ISMC_ee) -d_+ b1 * ptask.curvel_ + b2*ptask.actual_  ;

return command_;
}



static jspace::Vector jjTask( double kkd, double kkp, double kmaxvel){

 jspace::Vector kp= kkp * jspace::Vector::Ones(8);
  jspace::Vector kd = kkd * jspace::Vector::Ones(8);
   jspace::Vector maxvel = kmaxvel* jspace::Vector::Ones(8);
   jspace::Vector errpos_ = ptask.jgoalpos_-ptask.jjpos  ;
  jspace::Vector errvel_ = ptask.jgoalvel_ - ptask.jjvel  ;
  jspace::Vector pcommand_ = kp.cwise() * errpos_;
    
    
      for (int ii(0); ii < pcommand_.rows(); ++ii) {
	if ((maxvel[ii] > 1e-4) && (kd[ii] > 1e-4)) { // beware of div by zero
	  double const sat(fabs((pcommand_[ii] / maxvel[ii]) / kd[ii]));
	  if (sat > 1.0) {
	    pcommand_[ii] /= sat;
	  }
	}
      }
    
    pcommand_ += kd.cwise() * errvel_;

return pcommand_;


}

static jspace::Vector controller(jspace::Vector &dis, jspace::Vector &jdis, jspace::Vector &f, jspace::Vector &p){
    

  jspace::Vector  cc, gg;
  model->getGravity(gg);
  model->getCoriolisCentrifugal(cc);
  jspace::Matrix ainv;
  model->getInverseMassInertia(ainv);


jspace::Matrix lambda_ = ptask.getJacobian() * ainv * ptask.getJacobian().transpose();
jspace::Matrix   jbar_ = ainv * ptask.getJacobian().transpose() * lambda_.inverse();
 jspace::Matrix  nullspace_ = jspace::Matrix::Identity(8, 8) - ptask.getJacobian().transpose() * jbar_.transpose();
jspace::Vector command1 = ptask.getJacobian().transpose() * (f + dis) + nullspace_*(p+jdis) + cc + gg;
    
return command1;
}



static bool servo_cb(size_t toggle_count,
		     double wall_time_ms,
		     double sim_time_ms,
		     jspace::State & state,
		     jspace::Vector & command)
{
  mode = toggle_count % 5;
  static size_t prevmode(42);
  jspace::Vector jpos(state.position_.rows());
  jspace::Vector jvel(state.velocity_.rows());
  for (int ii(0); ii < state.position_.rows(); ++ii) {
    static double const amplitude(0.5 * M_PI);
    double const omega(1.0 + 0.1 * ii);
    double const phase(M_PI/2);
    jpos[ii] =         amplitude * sin(phase);
    jvel[ii] = omega * amplitude * cos(phase);  
}
  if (0 == mode) {
    
    state.position_ = jpos;
    state.velocity_ = jvel;
    
    prevmode = mode;
    return false;
    
  }
  
  //////////////////////////////////////////////////
  // Update the model to reflect the current robot state.
  
  model->update(state);
  
  if (0 == prevmode) {
    jspace::Status const st(ptask.init(*model));
    if ( ! st) {
      errx(EXIT_FAILURE, "ptask.init() failed: %s", st.errstr.c_str());
    }
  } 
  
  ptask.update(*model);

  static jspace::Vector pos(2), vel(2);
  static double const ox(0.2);
  static double const oy(0.37);
  static double const amp(2.5);
  double const px(ox * 1e-3 * wall_time_ms);
  double const py(oy * 1e-3 * wall_time_ms);
  pos <<      amp * sin(px),      amp * sin(py);
  vel << ox * amp * cos(px), oy * amp * cos(py);
 
  /*pos <<      2.5,      2.5;
  vel <<0, 0;*/

  if ( ! goalpos->set(pos)) {
    errx(EXIT_FAILURE, "failed to set end-effector goal position");
  }
  if ( ! goalvel->set(vel)) {
    errx(EXIT_FAILURE, "failed to set end-effector goal velocity");
  }
 if ( ! jgoalpos->set(jpos)) {
      errx(EXIT_FAILURE, "failed to set joint goal position");
    }
    if ( ! jgoalvel->set(jvel)) {
      errx(EXIT_FAILURE, "failed to set joint goal velocity");
    }
  
  
  static jspace::Vector dis, disjj, dtz, djz;
  dis = 200*sin(.001*10*wall_time_ms) *jspace::Vector::Ones(2) ;
  disjj =200*sin(.001*10*wall_time_ms)*jspace::Vector::Ones(8);
dtz = jspace::Vector::Zero(2);
djz =jspace::Vector::Zero(8);
jspace::Vector error1 = error();
jspace::Vector errorpx = error1.block(0,0,1,1);
jspace::Vector errorpy=  error1.block(1,0,1,1);

jspace::Vector p1, p0;
p1= jjTask(100,20,M_PI);


  switch (mode) {
  case 1:
  jspace::Vector f = pdTask(2500,40);
  command= controller(dtz, djz, f, p1);
  break;
  case 2:
  jspace::Vector f = pdTask(2500,40);
  command= controller(dis, disjj, f, p1);
  break;
  case 3:
  jspace::Vector f = ismcTask(2500,40,40,2500,-100);
  command= controller(dis, disjj, f, p1);
  break;
  case 4:
  jspace::Vector f = eeTask(2500,40,40,2500,-100);
  command= controller(dis, disjj, f, p1);
  break;
  }

jspace::Vector eq = ptask.jgoalpos_ - state.position_;
  static size_t iteration(0);
  if (0 == (iteration % 100)) {
    std::cerr << "**************************************************\n";
    jspace::pretty_print(state.position_, std::cerr, "  pposx", "    ");
    jspace::pretty_print(ptask.jgoalpos_, std::cerr, "  pposx", "    ");
    jspace::pretty_print(command, std::cerr, "  pposx", "    ");

//Record output data   
std::ofstream test2res;
std::ostream *os;
   test2res.open ("pddispossin.txt", std::ios::app);
    os = &test2res;
     *os << jspace::pretty_string(errorpx)<< "  "
     << jspace::pretty_string(errorpy)<< "  "
     << jspace::pretty_string(f1)<< "  "
	<< jspace::pretty_string(eq)<< "\n ";
   
  }
  ++iteration;
  prevmode = mode;
  return true;
}


static void draw_cb(double x0, double y0, double scale)
{
  if (0 != mode) {
    
    fl_color(255, 100, 100);
    fl_line_style(FL_SOLID, 3, 0);
    
    double const gx(goalpos->getVector()->x());
    double const gy(goalpos->getVector()->y());
    int const rr(ceil(0.2 * scale));
    int const dd(2 * rr);
    fl_arc(int(x0 + gx * scale) - rr, int(y0 - gy * scale) - rr, dd, dd, 0.0, 360.0);
    
    double const vx(goalvel->getVector()->x());
    double const vy(goalvel->getVector()->y());
    double const px(gx + vx * 0.1);
    double const py(gy + vy * 0.1);
    // fl_line(x0 + (gx + 0.2) * scale, y0 - gy * scale,
    // 	    x0 + (gx - 0.2) * scale, y0 - gy * scale);
    // fl_line(x0 + gx * scale, y0 - (gy + 0.2) * scale,
    // 	    x0 + gx * scale, y0 - (gy - 0.2) * scale);
    fl_color(255, 255, 100);
    fl_line(x0 + gx * scale, y0 - gy * scale,
	    x0 + px * scale, y0 - py * scale);
  }
}


int main(int argc, char ** argv)
{
  try {
    
    
    	model.reset(jspace::test::parse_sai_xml_file(model_filena	me, 	true));
    goalpos = ptask.lookupParameter("goalpos", 	opspace::PARAMETER_TYPE_VECTOR);
    if ( ! goalpos) {
      errx(EXIT_FAILURE, "failed to find appropriate goalpos 	parameter");
    }
    goalvel = ptask.lookupParameter("goalvel", 	opspace::PARAMETER_TYPE_VECTOR);
    if ( ! goalvel) {
      errx(EXIT_FAILURE, "failed to find appropriate goalvel 	parameter");
   }
    
   jgoalpos = ptask.lookupParameter("jgoalpos", 	opspace::PARAMETER_TYPE_VECTOR);
    if ( ! jgoalpos) {
      errx(EXIT_FAILURE, "failed to find appropriate joint-	posture goalpos parameter");
    }
    jgoalvel = ptask.lookupParameter("jgoalvel", 	opspace::PARAMETER_TYPE_VECTOR);
    if ( ! jgoalvel) {
      errx(EXIT_FAILURE, "failed to find appropriate joint-	posture goalvel parameter");
    }

   
  }
  catch (std::runtime_error const & ee) {
    errx(EXIT_FAILURE, "%s", ee.what());
  }
  tutsim::set_draw_cb(draw_cb);
  return tutsim::run(servo_cb);
}
