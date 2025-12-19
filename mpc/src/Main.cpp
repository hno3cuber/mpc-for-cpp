#include "MPCCal.h"
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

int main() {
	unique_ptr<GetRefPath> ref_ptr = make_unique<GetRefPath>("ref.csv");

	VectorXd sta_local = VectorXd::Zero(STANUM);
	vector<double> x;
	vector<double> y;

	uint size = ref_ptr->path.size();
	x.reserve(size);
	y.reserve(size);

	unique_ptr<MPCCal> mpc_ptr = make_unique<MPCCal>(sta_local);

	if (mpc_ptr->MPCSolver(ref_ptr.get())) {
		cout << "solve complete" << endl;
	}

	for (int i = 0; i < size; i++) {
		x.push_back(ref_ptr->path[i][0]);
		y.push_back(ref_ptr->path[i][1]);
	}

	plt::figure_size(1200, 780);
	plt::plot(x, y, "k--");

	x.clear();
	y.clear();
	for (int i = 0; i < size; i++) {
		x.push_back(mpc_ptr->status_node[i][0]);
		y.push_back(mpc_ptr->status_node[i][1]);
	}

	plt::plot(x,y,"r-");
	
	plt::show();
}