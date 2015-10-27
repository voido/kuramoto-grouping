#include <vector>

#include <boost/numeric/odeint.hpp>
#include <boost/operators.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real.hpp>

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

double freq_from_id(int i, int layer) {
	return M_PI/2.0 + (double) i/(double) layer;
}


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
				support[x[j].frequency] += f(i,j) * 0.5 * (cos(x[j].phase-mx.phase)+1.0);
			}
			mdxdt.phase = freq_from_id(mx.frequency, layer) + K*tmp; 
			int nf = std::distance(support.begin(), std::max_element(support.begin(), support.end()));
			x[i].frequency = nf;
		}
	}

	double f(size_t i, size_t j) {
		if(i<num/2 && j<num/2)
			return 1.0;
		if(i>=num/2 && j>=num/2)
			return 1.0;

		return -1.0;
	}

	size_t num;
	size_t layer;
	double K;
	std::vector<double> support;
};


struct gnuplot_observer {
  gnuplot_observer(int l) : layer(l) {
		std::cout << "set term x11" << std::endl;
		std::cout << "set polar" << std::endl;
		std::cout << "set size square" << std::endl;
		std::cout << "set grid polar" << std::endl;
		std::cout << "set xrange [-4:4]" << std::endl;
		std::cout << "set yrange [-4:4]" << std::endl;
	}

  template<class State>
	void operator()(State &x, double t) {
		std::cout << "plot '-' u 1:2:3 lc var lt 7 t 'time " << t << "'" << std::endl;
		size_t n = x.size();
		for(size_t i=0; i<n; i++) {
			int label = (i<n/2) ? 1 : 3;
			std::cout << x[i].phase << " " << freq_from_id(x[i].frequency, layer) << " " << label << std::endl;
		}
		std::cout << "e" << std::endl;
  }
	int layer;
};


int main(int argc, char* argv[]) {
	
	size_t units = 100;
	size_t layer = 10;
	double K = 2.0;

	osc_ensemble net(units, layer, K);
	gnuplot_observer obs(layer);
	
	boost::random::mt19937 rng;
	rng.seed(static_cast<unsigned int>(std::time(0)));
	boost::uniform_real<> phase_dist(0.0, 2.0*M_PI);
	boost::random::uniform_int_distribution<> freq_dist(0,layer-1);

	state_type x;
	for(size_t i=0; i<units; i++) {
		osc o;
		o.phase = phase_dist(rng);
		o.frequency = freq_dist(rng);
		x.push_back(o);
	}

	boost::numeric::odeint::runge_kutta4<state_type> stepper;
	boost::numeric::odeint::integrate_const(stepper, boost::ref(net), x, 0.0, 100.0, 0.01, boost::ref(obs));

	return 0;
}
