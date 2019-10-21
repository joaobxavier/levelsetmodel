package edu.harvard.sysbio.gavea.finitedifferencing;

import edu.harvard.sysbio.gavea.utils.ScalarMath;

/**
 * Implements numeric differentiation using WENO scheme (see Osher and Fedwik
 * level set book)
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 25, 2007
 */
public class WenoScheme {
	private double[][] _gradI;

	private double[][] _gradJ;

	public WenoScheme(double[][] u, double dx) {
		int n = u.length;
		_gradI = new double[n][n];
		_gradJ = new double[n][n];
		gradient(u, dx, _gradI, _gradJ);
	}

	public static void centralGradient(double[][] u, double dx,
			double[][] gradI, double[][] gradJ) {
		int n = u.length;
		int i, j;
		int iPlus, iMinus, jPlus, jMinus;
		double aux = 1 / (2 * dx);
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				iPlus = ScalarMath.cyclicUpwards(i + 1, n);
				iMinus = ScalarMath.cyclicDownwards(i - 1, n);
				jPlus = ScalarMath.cyclicUpwards(j + 1, n);
				jMinus = ScalarMath.cyclicDownwards(j - 1, n);
				gradI[i][j] = (u[iPlus][j] - u[iMinus][j]) * aux;
				gradJ[i][j] = (u[i][jPlus] - u[i][jMinus]) * aux;
			}
		}
	}

	public static void gradient(double[][] u, double dx, double[][] gradI,
			double[][] gradJ) {
		int n = u.length;
		int i, j;
		double iPlus3, iPlus2, iPlus1, iMinus1, iMinus2, iMinus3;
		double jPlus3, jPlus2, jPlus1, jMinus1, jMinus2, jMinus3;
		double v;
		double invDx = 1 / dx;
		double tmpGradMinus, tmpGradPlus;
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				v = u[i][j];
				// i's
				iPlus3 = u[ScalarMath.cyclicUpwards(i + 3, n)][j];
				iPlus2 = u[ScalarMath.cyclicUpwards(i + 2, n)][j];
				iPlus1 = u[ScalarMath.cyclicUpwards(i + 1, n)][j];
				iMinus1 = u[ScalarMath.cyclicDownwards(i - 1, n)][j];
				iMinus2 = u[ScalarMath.cyclicDownwards(i - 2, n)][j];
				iMinus3 = u[ScalarMath.cyclicDownwards(i - 3, n)][j];
				//
				tmpGradMinus = weno(iMinus3, iMinus2, iMinus1, v, iPlus1,
						iPlus2)
						* invDx;

				tmpGradPlus = -weno(iPlus3, iPlus2, iPlus1, v, iMinus1, iMinus2)
						* invDx;
				gradI[i][j] = -(tmpGradPlus + tmpGradMinus) * 0.5;
				// j's
				jPlus3 = u[i][ScalarMath.cyclicUpwards(j + 3, n)];
				jPlus2 = u[i][ScalarMath.cyclicUpwards(j + 2, n)];
				jPlus1 = u[i][ScalarMath.cyclicUpwards(j + 1, n)];
				jMinus1 = u[i][ScalarMath.cyclicDownwards(j - 1, n)];
				jMinus2 = u[i][ScalarMath.cyclicDownwards(j - 2, n)];
				jMinus3 = u[i][ScalarMath.cyclicDownwards(j - 3, n)];
				//
				tmpGradMinus = weno(jMinus3, jMinus2, jMinus1, v, jPlus1,
						jPlus2)
						* invDx;
				tmpGradPlus = -weno(jPlus3, jPlus2, jPlus1, v, jMinus1, jMinus2)
						* invDx;
				gradJ[i][j] = -(tmpGradPlus + tmpGradMinus) * 0.5;
			}
		}
	}

	private static double v1, v2, v3, v4, v5;

	private static double diff1, diff2, diff3;

	private static double w1, w2, w3;

	private static double s1, s2, s3;

	private static double a1, a2, a3;

	private static double eps;

	private final static double ONETHIRD = 1.0 / 3.0;

	private final static double ONESIXTH = 1.0 / 6.0;

	private final static double FIVESIXTHS = 5.0 / 6.0;

	private final static double SEVENSIXTHS = 7.0 / 6.0;

	private final static double ELEVENSIXTHS = 11.0 / 6.0;

	private final static double ONEFORTH = 1.0 / 4.0;

	private final static double THIRTEENTWELVTHS = 13.0 / 12.0;

	private static double weno(double i1, double i2, double i3, double i4,
			double i5, double i6) {
		// v's
		v1 = (i2 - i1);
		v2 = (i3 - i2);
		v3 = (i4 - i3);
		v4 = (i5 - i4);
		v5 = (i6 - i5);
		// diff's
		diff1 = ONETHIRD * v1 - SEVENSIXTHS * v2 + ELEVENSIXTHS * v3;
		diff2 = -ONESIXTH * v2 + FIVESIXTHS * v3 + ONETHIRD * v4;
		diff3 = ONETHIRD * v3 + FIVESIXTHS * v4 - ONESIXTH * v5;
		// s's
		s1 = THIRTEENTWELVTHS * ScalarMath.square(v1 - 2 * v2 + v3) + ONEFORTH
				* ScalarMath.square(v1 - 4 * v2 + 3 * v3);
		s2 = THIRTEENTWELVTHS * ScalarMath.square(v2 - 2 * v3 + v4) + ONEFORTH
				* ScalarMath.square(v2 - v4);
		s3 = THIRTEENTWELVTHS * ScalarMath.square(v3 - 2 * v4 + v5) + ONEFORTH
				* ScalarMath.square(3 * v3 - 4 * v4 + 3 * v5);
		// alphas's
		eps = 1e-6*maxVSquared(v1, v2, v3, v4, v5) + 1e-99;
		a1 = 0.1 / ScalarMath.square(s1 + eps);
		a2 = 0.6 / ScalarMath.square(s2 + eps);
		a3 = 0.3 / ScalarMath.square(s3 + eps);
		// weights
		double aux = a1 + a2 + a3;
		w1 = a1 / aux;
		w2 = a2 / aux;
		w3 = a3 / aux;
		// finally:
		return w1 * diff1 + w2 * diff2 + w3 * diff3;
	}

	private static double maxVSquared(double i1, double i2, double i3,
			double i4, double i5) {
		i1 = ScalarMath.square(i1);
		i2 = ScalarMath.square(i2);
		i3 = ScalarMath.square(i3);
		i4 = ScalarMath.square(i4);
		i5 = ScalarMath.square(i5);
		i1 = (i1 > i2 ? i1 : i2);
		i1 = (i1 > i3 ? i1 : i3);
		i1 = (i1 > i4 ? i1 : i4);
		i1 = (i1 > i5 ? i1 : i5);
		return i1;
	}

	// ///////////////////////////////////////////////////////////
	// Static methods for WENO scheme
	private static double i1, i2, i3, i4, i5, i6;

	public static double wenoIMinus(double[][] u, int n, int i, int j) {
		i1 = u[ScalarMath.cyclicDownwards(i - 3, n)][j];
		i2 = u[ScalarMath.cyclicDownwards(i - 2, n)][j];
		i3 = u[ScalarMath.cyclicDownwards(i - 1, n)][j];
		i4 = u[i][j];
		i5 = u[ScalarMath.cyclicUpwards(i + 1, n)][j];
		i6 = u[ScalarMath.cyclicUpwards(i + 2, n)][j];
		return weno(i1, i2, i3, i4, i5, i6);
	}

	public static double wenoIPlus(double[][] u, int n, int i, int j) {
		i6 = u[ScalarMath.cyclicDownwards(i - 2, n)][j];
		i5 = u[ScalarMath.cyclicDownwards(i - 1, n)][j];
		i4 = u[i][j];
		i3 = u[ScalarMath.cyclicUpwards(i + 1, n)][j];
		i2 = u[ScalarMath.cyclicUpwards(i + 2, n)][j];
		i1 = u[ScalarMath.cyclicUpwards(i + 3, n)][j];
		return -weno(i1, i2, i3, i4, i5, i6);
	}

	public static double wenoJMinus(double[][] u, int n, int i, int j) {
		i1 = u[i][ScalarMath.cyclicDownwards(j - 3, n)];
		i2 = u[i][ScalarMath.cyclicDownwards(j - 2, n)];
		i3 = u[i][ScalarMath.cyclicDownwards(j - 1, n)];
		i4 = u[i][j];
		i5 = u[i][ScalarMath.cyclicUpwards(j + 1, n)];
		i6 = u[i][ScalarMath.cyclicUpwards(j + 2, n)];
		return weno(i1, i2, i3, i4, i5, i6);
	}

	public static double wenoJPlus(double[][] u, int n, int i, int j) {
		i6 = u[i][ScalarMath.cyclicDownwards(j - 2, n)];
		i5 = u[i][ScalarMath.cyclicDownwards(j - 1, n)];
		i4 = u[i][j];
		i3 = u[i][ScalarMath.cyclicUpwards(j + 1, n)];
		i2 = u[i][ScalarMath.cyclicUpwards(j + 2, n)];
		i1 = u[i][ScalarMath.cyclicUpwards(j + 3, n)];
		return -weno(i1, i2, i3, i4, i5, i6);
	}

	public double[][] getGradI() {
		return _gradI;
	}

	public double[][] getGradJ() {
		return _gradJ;
	}
}
