package edu.harvard.sysbio.gavea.finitedifferencing;

import edu.harvard.sysbio.gavea.utils.MatrixMath;
import edu.harvard.sysbio.gavea.utils.ScalarMath;

/**
 * Implements finite difference schemes for first order gradients
 * 
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 23, 2007
 */
public class GradientHandler {
	private double[][] _levelSet;

	private double _dx;

	private int _iSize;

	private double[][] _temp;

	// TODO delete when done with testing
	private double[][] _temp2;

	private double[][] _original;

	public GradientHandler(int n, double dx) {
		_dx = dx;
		_iSize = n;
		_temp = new double[n][n];
		_original = new double[n][n];
		// TODO remove this variable
		_temp2 = new double[n][n];
	}

	/**
	 * set the Level-set matrix
	 * 
	 * @param ls
	 */
	public void setLevelSet(double[][] ls) {
		_levelSet = ls;
	}

	/**
	 * Compute the velocity using a central difference second order scheme with
	 * ghost node
	 * 
	 * @param p
	 *            the pressure matrix
	 * @param speedI
	 * @param speedJ
	 */
	public void computeSpeedCentralDifference(double[][] p, double[][] speedI,
			double[][] speedJ) {
		int i, j;
		double aux = 1 / _dx;
		int n = _iSize;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
				if (_levelSet[i][j] < 0) {
					int iPlus = ScalarMath.cyclicUpwards(i + 1, n);
					int iMinus = ScalarMath.cyclicDownwards(i - 1, n);
					int jPlus = ScalarMath.cyclicUpwards(j + 1, n);
					int jMinus = ScalarMath.cyclicDownwards(j - 1, n);
					// i direction
					speedI[i][j] = -(p[iPlus][j] - p[iMinus][j]) * 0.5 * aux;
					// j-direction
					speedJ[i][j] = -(p[i][jPlus] - p[i][jMinus]) * 0.5 * aux;
				} else {
					// outside biomass
					speedI[i][j] = 0;
					speedJ[i][j] = 0;
				}
			}
	}

	/**
	 * Compute the gradient of some data using central difference
	 * 
	 * @param p
	 * @param speedI
	 * @param speedJ
	 * @param n
	 * @param dx
	 */
	public static void computeSpeedCentralDifferenceGradient(double[][] p,
			double[][] speedI, double[][] speedJ, int n, double dx) {
		int i, j;
		double aux = 1 / dx;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
				int iPlus = ScalarMath.cyclicUpwards(i + 1, n);
				int iMinus = ScalarMath.cyclicDownwards(i - 1, n);
				int jPlus = ScalarMath.cyclicUpwards(j + 1, n);
				int jMinus = ScalarMath.cyclicDownwards(j - 1, n);
				// i direction
				speedI[i][j] = (p[iPlus][j] - p[iMinus][j]) * 0.5 * aux;
				// j-direction
				speedJ[i][j] = (p[i][jPlus] - p[i][jMinus]) * 0.5 * aux;
			}
	}

	/**
	 * Compute the velocity using a central difference second order scheme with
	 * ghost node
	 * 
	 * @param p
	 *            the pressure matrix
	 * @param speedI
	 * @param speedJ
	 */
	public void computeSpeedUpwindGhostNodeCentralDifference(double[][] p,
			double[][] speedI, double[][] speedJ) {
		int i, j;
		double aux = 1 / _dx;
		int n = _iSize;
		double speedIPlus = 0;
		double speedIMinus = 0;
		double speedJPlus = 0;
		double speedJMinus = 0;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
				if (_levelSet[i][j] < 0) {
					int iPlus = ScalarMath.cyclicUpwards(i + 1, n);
					int iMinus = ScalarMath.cyclicDownwards(i - 1, n);
					int jPlus = ScalarMath.cyclicUpwards(j + 1, n);
					int jMinus = ScalarMath.cyclicDownwards(j - 1, n);
					if (isBorder(i, j, iPlus, iMinus, jPlus, jMinus)) {
						// speedIPlus = ghostNodeDifference(p, i, j, iPlus, j);
						// speedIMinus = -ghostNodeDifference(p, i, j, iMinus,
						// j);
						// speedJPlus = ghostNodeDifference(p, i, j, i, jPlus);
						// speedJMinus = -ghostNodeDifference(p, i, j, i,
						// jMinus);
						// // i direction
						// speedI[i][j] = atributeUpwindSpeed(speedIPlus,
						// speedIMinus)
						// * aux;
						// // j-direction
						// speedJ[i][j] = atributeUpwindSpeed(speedJPlus,
						// speedJMinus)
						// * aux;
						// ghost node second order
						double thetaIPlus = ghostNodeTheta(p, i, j, iPlus, j);
						double thetaIMinus = ghostNodeTheta(p, i, j, iMinus, j);
						double thetaJPlus = ghostNodeTheta(p, i, j, i, jPlus);
						double thetaJMinus = ghostNodeTheta(p, i, j, i, jMinus);
						// i direction
						speedI[i][j] = -(p[iPlus][j] - p[iMinus][j]) * aux
								/ (thetaIPlus + thetaIMinus);
						// j-direction
						speedJ[i][j] = -(p[i][jPlus] - p[i][jMinus]) * aux
								/ (thetaJPlus + thetaJMinus);
					} else {
						// // double speedIPlus = -(u[iPlus][j] - u[i][j]);
						// // double speedIMinus = -(u[i][j] - u[iMinus][j]);
						// // double speedJPlus = -(u[i][jPlus] - u[i][j]);
						// // double speedJMinus = -(u[i][j] - u[i][jMinus]);
						// speedIPlus = -WenoScheme.wenoIPlus(p, n, i, j);
						// speedIMinus = -WenoScheme.wenoIMinus(p, n, i, j);
						// speedJPlus = -WenoScheme.wenoJPlus(p, n, i, j);
						// speedJMinus = -WenoScheme.wenoJMinus(p, n, i, j);
						// // i direction
						// speedI[i][j] = atributeUpwindSpeed(speedIPlus,
						// speedIMinus)
						// * aux;
						// // j-direction
						// speedJ[i][j] = atributeUpwindSpeed(speedJPlus,
						// speedJMinus)
						// * aux;
						// second order
						// i direction
						speedI[i][j] = -(p[iPlus][j] - p[iMinus][j]) * 0.5
								* aux;
						// j-direction
						speedJ[i][j] = -(p[i][jPlus] - p[i][jMinus]) * 0.5
								* aux;
					}
				} else {
					// outside biomass
					speedI[i][j] = 0;
					speedJ[i][j] = 0;
				}
			}
	}

	/**
	 * Compute the velocity using a central difference second order scheme with
	 * ghost node
	 * 
	 * @param p
	 *            the pressure matrix
	 * @param speedI
	 * @param speedJ
	 */
	public void computeSpeedUpwindGhostNodeFirstOrder(double[][] p,
			double[][] speedI, double[][] speedJ) {
		int i, j;
		double aux = 1 / _dx;
		int n = _iSize;
		double speedIPlus = 0;
		double speedIMinus = 0;
		double speedJPlus = 0;
		double speedJMinus = 0;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
				if (_levelSet[i][j] < 0) {
					int iPlus = ScalarMath.cyclicUpwards(i + 1, n);
					int iMinus = ScalarMath.cyclicDownwards(i - 1, n);
					int jPlus = ScalarMath.cyclicUpwards(j + 1, n);
					int jMinus = ScalarMath.cyclicDownwards(j - 1, n);
					if (isBorder(i, j, iPlus, iMinus, jPlus, jMinus)) {
						// Ghost node first order
						speedIPlus = ghostNodeDifference(p, i, j, iPlus, j);
						speedIMinus = -ghostNodeDifference(p, i, j, iMinus, j);
						speedJPlus = ghostNodeDifference(p, i, j, i, jPlus);
						speedJMinus = -ghostNodeDifference(p, i, j, i, jMinus);
						// i direction
						speedI[i][j] = atributeUpwindSpeed(speedIPlus,
								speedIMinus)
								* aux;
						// j-direction
						speedJ[i][j] = atributeUpwindSpeed(speedJPlus,
								speedJMinus)
								* aux;
						// // ghost node second order
						// double thetaIPlus = ghostNodeTheta(p, i, j, iPlus,
						// j);
						// double thetaIMinus = ghostNodeTheta(p, i, j, iMinus,
						// j);
						// double thetaJPlus = ghostNodeTheta(p, i, j, i,
						// jPlus);
						// double thetaJMinus = ghostNodeTheta(p, i, j, i,
						// jMinus);
						// // i direction
						// speedI[i][j] = -(p[iPlus][j] - p[iMinus][j]) * aux
						// / (thetaIPlus + thetaIMinus);
						// // j-direction
						// speedJ[i][j] = -(p[i][jPlus] - p[i][jMinus]) * aux
						// / (thetaJPlus + thetaJMinus);
					} else {
						// // double speedIPlus = -(u[iPlus][j] - u[i][j]);
						// // double speedIMinus = -(u[i][j] - u[iMinus][j]);
						// // double speedJPlus = -(u[i][jPlus] - u[i][j]);
						// // double speedJMinus = -(u[i][j] - u[i][jMinus]);
						// speedIPlus = -WenoScheme.wenoIPlus(p, n, i, j);
						// speedIMinus = -WenoScheme.wenoIMinus(p, n, i, j);
						// speedJPlus = -WenoScheme.wenoJPlus(p, n, i, j);
						// speedJMinus = -WenoScheme.wenoJMinus(p, n, i, j);
						// // i direction
						// speedI[i][j] = atributeUpwindSpeed(speedIPlus,
						// speedIMinus)
						// * aux;
						// // j-direction
						// speedJ[i][j] = atributeUpwindSpeed(speedJPlus,
						// speedJMinus)
						// * aux;
						// second order
						// i direction
						speedI[i][j] = -(p[iPlus][j] - p[iMinus][j]) * 0.5
								* aux;
						// j-direction
						speedJ[i][j] = -(p[i][jPlus] - p[i][jMinus]) * 0.5
								* aux;
					}
				} else {
					// outside biomass
					speedI[i][j] = 0;
					speedJ[i][j] = 0;
				}
			}
	}

	/**
	 * @param p
	 * @param i
	 * @param j
	 * @param iNeib
	 * @param jNeib
	 * @return the theta value
	 */
	private final double ghostNodeTheta(double[][] p, int i, int j, int iNeib,
			int jNeib) {
		if (_levelSet[iNeib][jNeib] > 0) {
			return _levelSet[i][j]
					/ (_levelSet[i][j] - _levelSet[iNeib][jNeib]);
		}
		return 1.0;
	}

	private double ghostNodeDifference(double[][] p, int i, int j, int iNeib,
			int jNeib) {
		if (_levelSet[iNeib][jNeib] > 0) {
			double theta = _levelSet[i][j]
					/ (_levelSet[i][j] - _levelSet[iNeib][jNeib]);
			return p[i][j] / theta;
		}
		return p[i][j] - p[iNeib][jNeib];
	}

	/**
	 * @param i
	 * @param j
	 * @param iPlus
	 * @param iMinus
	 * @param jPlus
	 * @param jMinus
	 * @return true if any of the neighboring points are outside
	 */
	private boolean isBorder(int i, int j, int iPlus, int iMinus, int jPlus,
			int jMinus) {
		return (_levelSet[iPlus][j] >= 0) | (_levelSet[iMinus][j] >= 0)
				| (_levelSet[i][jPlus] >= 0) | (_levelSet[i][jMinus] >= 0);
	}

	/**
	 * Compute the velocity using a central difference second order scheme with
	 * ghost node
	 * 
	 * @param speedI
	 * @param speedJ
	 */
	public void computeSpeedUpwind(double[][] u, double[][] speedI,
			double[][] speedJ) {
		int i, j;
		double aux = 1 / _dx;
		int n = _iSize;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++) {
				if (_levelSet[i][j] < 0) {
					// int iPlus = ScalarMath.cyclicUpwards(i + 1, n);
					// int iMinus = ScalarMath.cyclicDownwards(i - 1, n);
					// int jPlus = ScalarMath.cyclicUpwards(j + 1, n);
					// int jMinus = ScalarMath.cyclicDownwards(j - 1, n);
					// double speedIPlus = -(u[iPlus][j] - u[i][j]);
					// double speedIMinus = -(u[i][j] - u[iMinus][j]);
					// double speedJPlus = -(u[i][jPlus] - u[i][j]);
					// double speedJMinus = -(u[i][j] - u[i][jMinus]);
					double speedIPlus = -WenoScheme.wenoIPlus(u, n, i, j);
					double speedIMinus = -WenoScheme.wenoIMinus(u, n, i, j);
					double speedJPlus = -WenoScheme.wenoJPlus(u, n, i, j);
					double speedJMinus = -WenoScheme.wenoJMinus(u, n, i, j);

					speedI[i][j] = atributeUpwindSpeed(speedIPlus, speedIMinus)
							* aux;
					// j-direction
					speedJ[i][j] = atributeUpwindSpeed(speedJPlus, speedJMinus)
							* aux;
				} else {
					speedI[i][j] = 0;
					speedJ[i][j] = 0;
				}
			}
	}

	private static final double atributeUpwindSpeed(double speedPlus,
			double speedMinus) {
		// return Math.min(speedPlus, 0) + Math.max(speedMinus, 0); // new
		return Math.max(speedPlus, 0) + Math.min(speedMinus, 0); // new
	}

	/**
	 * smoothen using Runge-Kutta
	 */
	public double[][] extendFirstWithRungeKutta(double[][] u) {
		for (int c = 0; c < 1; c++) {
			// for (int c = 0; c < 30; c++) {
			MatrixMath.copyTo(u, _original);
			// advance twice, _levelset will hold n+2 temporary solution
			extendFirst(u);
			extendFirst(u);
			// linear combination, placing the result (n+3/2) in _levelSet
			MatrixMath.linearCombination(_original, 3.0 / 4.0, u, 1.0 / 4.0, u);
			// advance once more, _levelset will hold n+3/2 temporary solution
			extendFirst(u);
			// linear combination, placing the result (n+3/2) in _levelSet
			MatrixMath.linearCombination(_original, 1.0 / 3.0, u, 2.0 / 3.0, u);
		}
		return u;
	}

	/**
	 * Extend velocity to the outside of the biofilm using first order upwind
	 * gradient scheme
	 * 
	 * @param u
	 */
	public void extendFirst(double[][] u) {
		int i, j, c;
		int n = _iSize;
		double cfl = 0.25; // CFL is 0.25 since it combines 0.5 (the
		// recommended
		// value) with the central difference scheme
		// ith the fact that central difference is used for level set
		// derivatives.
		// dx is omitted since it disappears in the fraction
		double dudI, dudJ, dLsdI, dLsdJ;
		for (c = 0; c < 2; c++) {
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					if (_levelSet[i][j] >= 0) {
						int iPlus = ScalarMath.cyclicUpwards(i + 1, n);
						int iMinus = ScalarMath.cyclicDownwards(i - 1, n);
						int jPlus = ScalarMath.cyclicUpwards(j + 1, n);
						int jMinus = ScalarMath.cyclicDownwards(j - 1, n);

						// central differences
						dLsdI = (_levelSet[iPlus][j] - _levelSet[iMinus][j]);
						dLsdJ = (_levelSet[i][jPlus] - _levelSet[i][jMinus]);

						// first order (update CFL to 0.5)
						// double dLsdIPlus = (_levelSet[iPlus][j] -
						// _levelSet[i][j]);
						// double dLsdIMinus = (_levelSet[i][j] -
						// _levelSet[iMinus][j]);
						// double dLsdJPlus = (_levelSet[i][jPlus] -
						// _levelSet[i][j]);
						// double dLsdJMinus = (_levelSet[i][j] -
						// _levelSet[i][jMinus]);

						// WENO (update CFL to 0.5)
						// double dLsdIPlus = WenoScheme.wenoIPlus(_levelSet, n,
						// i, j);
						// double dLsdIMinus = WenoScheme.wenoIMinus(_levelSet,
						// n,
						// i, j);
						// double dLsdJPlus = WenoScheme.wenoJPlus(_levelSet, n,
						// i, j);
						// double dLsdJMinus = WenoScheme.wenoJMinus(_levelSet,
						// n,
						// i, j);

						if (_levelSet[iMinus][j] < _levelSet[iPlus][j]) {
							dudI = (u[i][j] - u[iMinus][j]);
						} else {
							dudI = (u[iPlus][j] - u[i][j]);
						}
						if (_levelSet[i][jMinus] < _levelSet[i][jPlus]) {
							dudJ = (u[i][j] - u[i][jMinus]);
						} else {
							dudJ = (u[i][jPlus] - u[i][j]);
						}
						// account for the zero gradients
						double norm = Math.sqrt(ScalarMath.square(dLsdI)
								+ ScalarMath.square(dLsdJ));
						_temp[i][j] = norm != 0 ? u[i][j]
								- cfl
								* (dudI * dLsdI + dudJ * dLsdJ)
								/ Math.sqrt(ScalarMath.square(dLsdI)
										+ ScalarMath.square(dLsdJ)) * 0.5 : 0;
					}
				}
			}
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					if (_levelSet[i][j] >= 0) {
						u[i][j] = _temp[i][j];
					}
				}
			}

		}
	}

	/**
	 * Advance the levelset in time accorind to the speed vector using RK4
	 * 
	 * @param speedI
	 * @param speedJ
	 * @param dt
	 */
	public void advanceWithRungeKutta(double[][] speedI, double[][] speedJ,
			double dt) {
		MatrixMath.copyTo(_levelSet, _original);
		// advance twice, _levelset will hold n+2 temporary solution
		advance(speedI, speedJ, dt);
		advance(speedI, speedJ, dt);
		// linear combination, placing the result (n+3/2) in _levelSet
		MatrixMath.linearCombination(_original, 3.0 / 4.0, _levelSet,
				1.0 / 4.0, _levelSet);
		// advance once mor4, _levelset will hold n+3/2 temporary solution
		advance(speedI, speedJ, dt);
		// linear combination, placing the result (n+3/2) in _levelSet
		MatrixMath.linearCombination(_original, 1.0 / 3.0, _levelSet,
				2.0 / 3.0, _levelSet);
	}

	/**
	 * Advance the levelset in time accorind to the speed vector
	 * 
	 * @param speedI
	 * @param speedJ
	 * @param dt
	 */
	public void advance(double[][] speedI, double[][] speedJ, double dt) {
		int i, j;
		int n = _iSize;
		double dLsdI, dLsdJ;
		double invDx = 1 / _dx;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				// UPWIND first order
				int iPlus = (i + 1 < n) ? (i + 1) : 0;
				int iMinus = (i - 1 < 0) ? (n - 1) : (i - 1);
				int jPlus = (j + 1 < n) ? (j + 1) : 0;
				int jMinus = (j - 1 < 0) ? (n - 1) : (j - 1);
				double dLsdIPlus = (_levelSet[iPlus][j] - _levelSet[i][j]);
				double dLsdIMinus = (_levelSet[i][j] - _levelSet[iMinus][j]);
				double dLsdJPlus = (_levelSet[i][jPlus] - _levelSet[i][j]);
				double dLsdJMinus = (_levelSet[i][j] - _levelSet[i][jMinus]);
				// // WENO scheme
				// double dLsdIPlus = WenoScheme.wenoIPlus(_levelSet, n, i, j);
				// double dLsdIMinus = WenoScheme.wenoIMinus(_levelSet, n, i,
				// j);
				// double dLsdJPlus = WenoScheme.wenoJPlus(_levelSet, n, i, j);
				// double dLsdJMinus = WenoScheme.wenoJMinus(_levelSet, n, i,
				// j);

				// if (speedI[i][j] > 0) {
				// dLsdI = dLsdIMinus;
				// } else {
				// dLsdI = dLsdIPlus;
				// }
				// if (speedJ[i][j] > 0) {
				// dLsdJ = dLsdJMinus;
				// } else {
				// dLsdJ = dLsdJPlus;
				// }
				// // TODO there is a dx missing somwhere
				// _temp[i][j] = dt
				// * (speedI[i][j] * dLsdI + speedJ[i][j] * dLsdJ)
				// / Math.sqrt(ScalarMath.square(dLsdI)
				// + ScalarMath.square(dLsdJ));
				// //
				_temp[i][j] = dt
						* (Math.max(speedI[i][j], 0) * dLsdIMinus
								+ Math.min(speedI[i][j], 0) * dLsdIPlus
								+ Math.max(speedJ[i][j], 0) * dLsdJMinus + Math
								.min(speedJ[i][j], 0)
								* dLsdJPlus) * invDx;
			}
		}

		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				_levelSet[i][j] -= _temp[i][j];
			}
		}
	}

	public void reinitializeWithRungeKutta() {
		MatrixMath.copyTo(_levelSet, _original);
		// advance twice, _levelset will hold n+2 temporary solution
		reinitialize();
		reinitialize();
		// linear combination, placing the result (n+1/2) in _levelSet
		MatrixMath.linearCombination(_original, 3.0 / 4.0, _levelSet,
				1.0 / 4.0, _levelSet);
		// advance once mor4, _levelset will hold n+3/2 temporary solution
		reinitialize();
		// linear combination, placing the result (n+3/2) in _levelSet
		MatrixMath.linearCombination(_original, 1.0 / 3.0, _levelSet,
				2.0 / 3.0, _levelSet);
	}

	/**
	 * Advance the levelset in time accorind to the speed vector
	 * 
	 * @param speedI
	 * @param speedJ
	 * @param dt
	 */
	public void reinitialize() {
		int i, j;
		int n = _iSize;
		int c;
		// auciliary
		double a1 = 1 / _dx;
		double a2 = 1 / (2 * _dx);
		double a3 = ScalarMath.square(_dx);
		double cfl;
		for (c = 0; c < 1; c++) {
			cfl = 0;
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					double v = _levelSet[i][j];
					double iPlus = _levelSet[(i + 1 < n) ? (i + 1) : 0][j];
					double iMinus = _levelSet[(i - 1 < 0) ? (n - 1) : (i - 1)][j];
					double jPlus = _levelSet[i][(j + 1 < n) ? (j + 1) : 0];
					double jMinus = _levelSet[i][(j - 1 < 0) ? (n - 1)
							: (j - 1)];

					// no need to divide by 2*h since it is only used to
					// 
					double dLsdICentral = (iPlus - iMinus) * a2;
					double dLsdJCentral = (jPlus - jMinus) * a2;
					// double dLsdIPlus = (_levelSet[(i + 1 < n) ? (i + 1) :
					// 0][j] - _levelSet[i][j])
					// * a1;
					// double dLsdIMinus = (_levelSet[i][j] - _levelSet[(i - 1 <
					// 0) ? (n - 1)
					// : (i - 1)][j])
					// * a1;
					// double dLsdJPlus = (_levelSet[i][(j + 1 < n) ? (j + 1) :
					// 0] - _levelSet[i][j])
					// * a1;
					// double dLsdJMinus = (_levelSet[i][j] - _levelSet[i][(j -
					// 1 < 0) ? (n - 1)
					// : (j - 1)])
					// * a1;
					double dLsdIPlus = WenoScheme.wenoIPlus(_levelSet, n, i, j)
							* a1;
					double dLsdIMinus = WenoScheme.wenoIMinus(_levelSet, n, i,
							j)
							* a1;

					double dLsdJPlus = WenoScheme.wenoJPlus(_levelSet, n, i, j)
							* a1;
					double dLsdJMinus = WenoScheme.wenoJMinus(_levelSet, n, i,
							j)
							* a1;
					// equation 3.12 b from Alpkvist Licensiat thesis
					// to approximate sign
					double sign = v
							/ Math
									.sqrt(ScalarMath.square(v)
											+ (ScalarMath.square(dLsdICentral) + ScalarMath
													.square(dLsdJCentral)) * a3);
					// store the maximum sign value temporarily in cfl;
					cfl = Math.max(cfl, sign);
					// the increment
					_temp[i][j] = chooseReinitializationDerivative(v, sign,
							dLsdIMinus, dLsdIPlus, dLsdJMinus, dLsdJPlus);
					_temp2[i][j] = sign;
				}
			}
			cfl = _dx / (2.0 * cfl);
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					_levelSet[i][j] -= cfl * _temp[i][j];
				}
			}
		}
	}

	/**
	 * From Alpkvist licenciat thesis
	 * 
	 * @param v
	 * @param s
	 * @param dIMinus
	 * @param dIPlus
	 * @param dJMinus
	 * @param dJPlus
	 * @return
	 */
	private static double chooseReinitializationDerivative(double v, double s,
			double dIMinus, double dIPlus, double dJMinus, double dJPlus) {
		double t1 = ScalarMath.square(Math.max(dIMinus, 0))
				+ ScalarMath.square(Math.min(dIPlus, 0));
		double t2 = ScalarMath.square(Math.max(dIPlus, 0))
				+ ScalarMath.square(Math.min(dIMinus, 0));
		double t3 = ScalarMath.square(Math.max(dJMinus, 0))
				+ ScalarMath.square(Math.min(dJPlus, 0));
		double t4 = ScalarMath.square(Math.max(dJPlus, 0))
				+ ScalarMath.square(Math.min(dJMinus, 0));
		return Math.max(s, 0) * (Math.sqrt(t1 + t3) - 1) + Math.min(s, 0)
				* (Math.sqrt(t2 + t4) - 1);
	}

	// TODO delete after testing
	public double[][] extendWithFirstOrder(double[][] u) {
		extendFirst(u);
		return u;
	}

	public double[][] getLevelSet() {
		return _levelSet;
	}

	public double[][] getSign() {
		return _temp2;
	}
}
