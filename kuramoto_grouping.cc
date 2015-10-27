#include <vector>

#include <boost/numeric/odeint.hpp>

typedef std::vector<float> state_type;


struct osc_ensemble {


	void operator()(state_type &x, state_type &dxdt, double t) {

	}

};

int main(int argc, char* argv[]) {

	osc_ensemble net;
	state_type x;


	boost::numeric::odeint::runge_kutta4<state_type> stepper;
	boost::numeric::odeint::integrate_const(stepper, boost::ref(net), x, 0.0, 100.0, 0.01);

	return 0;
}
