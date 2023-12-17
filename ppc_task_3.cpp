#include <iostream>
#include <functional>
#include <vector>

double findM(const std::vector<double>& X, const std::function<double(double)>& f) {
	double M = 0;
	for (int i = 0; i < X.size() - 1; i++) {
		M = std::max(M, std::abs(f(X[i]) - f(X[i + 1])) / std::abs(X[i + 1] - X[i]));
	}
	return M;
}

double getm(double M, double r) {
	if (M) {
		return r * M;
	}
	return 1;
}

std::vector<double> getR(double m, const std::vector<double>& X, std::function<double(double)> f) {
	std::vector<double> R(X.size() - 1);
	for (int i = 0; i < X.size() - 1; i++) {
		R[i] = m * abs(X[i + 1] - X[i]) + (f(X[i]) - f(X[i + 1])) * (f(X[i]) - f(X[i + 1]))
			/ m / abs(X[i + 1] - X[i]) - 2 * (f(X[i]) + f(X[i + 1]));
	}
	return R;
}

double getXk_1(double xt_1, double xt, double r,
	double m, std::function<double(double)> f) {

	return (xt_1 + xt) / 2 - ((f(xt) - f(xt_1)) > 0 ? 1 : -1) * r / m * std::abs(f(xt) - f(xt_1)) / (2*r);
}

int getMax(const std::vector<double>& v) {
	int maxIndex = 0;
	for (int i = 0; i < v.size(); i++) {
		if (v[maxIndex] < v[i]) {
			maxIndex = i;
		}
	}
	return maxIndex;
}

double getMinSequential(std::function<double(double)> f, double leftBound,
	double rightBound, double eps, int maxIterations, double r) {

	std::vector<double> X = { leftBound, rightBound };

	double lastX = 0;

	for (int i = 0; i < maxIterations; i++) {
		double M = findM(X, f);
		double m = getm(M, r);
		std::vector<double> R = getR(m, X, f);
		int maxRindex = getMax(R);

		double xk_1 = getXk_1(X[maxRindex], X[maxRindex + 1], r, m, f);

		if (xk_1 < leftBound) {
			xk_1 = leftBound;
		}
		else if (xk_1 > rightBound) {
			xk_1 = rightBound;
		}

		lastX = xk_1;


		if (std::abs(xk_1 - X[maxRindex]) < eps) {
			return xk_1;
		}

		auto xk_1Place = std::upper_bound(X.begin(), X.end(), xk_1);

		X.insert(xk_1Place, xk_1);
	}

	return lastX;
}


int main() {
	auto f = [](double x) {return -x * x * x; };

	double val = getMinSequential(f, -5, 5, 0.0001, 10000, 2);

	std::cout << val << std::endl;
}

