#include <vector>

#include <boost/numeric/odeint.hpp>
#include <boost/operators.hpp>


struct osc :
	boost::additive1<osc,
	boost::additive2<osc, double,
	boost::multiplicative2<osc, double> > > {
	osc() : phase(0.0), frequency(0) {}
	osc(double p, int f) :phase(p), frequency(f) {}

	double phase;
	int frequency;

	osc& operator+=(const osc &o) {
		this->phase+=o.phase;
		return *this;
	}
 	
	osc operator*(const double d) {
		osc tmp;
		tmp.phase = d*this->phase;
		return tmp;
  }
 
	osc operator*(const osc &o) {
		osc tmp;
		tmp.phase = this->phase*o.phase;
		return tmp;
  }
 
	osc& operator*=(const double d) {
		this->phase *= d;
		return *this;
  }
};

typedef std::vector<osc> state_type;


struct osc_ensemble {
	osc_ensemble(size_t num_osc, size_t num_layer, double K) :
	num(num_osc), layer(num_layer), K(K) {
		support.resize(layer);
	}

	~osc_ensemble() {}

	void operator()(state_type &x, state_type &dxdt, double t) {
		size_t n = x.size();

		for(size_t i=0; i<n; i++) {
			const osc &mx = x[i];
			osc &mdxdt = dxdt[i];
			double tmp = 0.0;
			std::fill(support.begin(), support.end(), 0.0);

			for(size_t j=0; j<n; j++) {
				tmp += f(i,j)*sin(x[j].phase-mx.phase);
			}
			mdxdt.phase = K*tmp; 
			int nf = std::distance(support.begin(), std::max_element(support.begin(), support.end()));
			x[i].frequency = nf;
		}
	}

	double f(size_t i, size_t j) {
		return 1;
	}

	size_t num;
	size_t layer;
	double K;
	std::vector<double> support;
};

int main(int argc, char* argv[]) {
	
	size_t units = 100;
	size_t layer = 10;
	double K = 2.0;

	osc_ensemble net(units, layer, K);
	state_type x;


	boost::numeric::odeint::runge_kutta4<state_type> stepper;
	boost::numeric::odeint::integrate_const(stepper, boost::ref(net), x, 0.0, 100.0, 0.01);

	return 0;
}
