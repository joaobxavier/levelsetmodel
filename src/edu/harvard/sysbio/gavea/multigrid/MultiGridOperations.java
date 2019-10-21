package edu.harvard.sysbio.gavea.multigrid;

/**
 * Implements multigrid matrix operations such as restriction and interpolations
 * 
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 15, 2007
 */
public abstract class MultiGridOperations {

	/**
	 * @param n
	 * @return ((n - 1) / 2) + 1
	 */
	public static final int coarserSize(int n) {
		return ((n - 1) / 2) + 1;
	}

	/**
	 * @param n
	 * @return ((n - 1) * 2) + 1
	 */
	public static final int finnerSize(int n) {
		return ((n - 1) * 2) + 1;
	}

	/**
	 * @param n
	 * @return true is n is of the form 2^n + 1
	 */
	public static final boolean validMultigridSize(int n) {
		double n2 = 2;
		while (n2 < n) {
			if (n2 + 1 == n)
				return true;
			n2 *= 2;
		}
		return false;
	}

	/**
	 * Restricts the data in matrix uf to a grid one order coarser.
	 * 
	 * @param uc
	 *            coarser grid
	 * @param uf
	 *            finner grid
	 */
	public static void restrictCyclic8Point(double[][] uc, double[][] uf) {
		int nf = uf.length;
		int nc = uc.length;

		int ic, iif, jc, jf;
		// Interior points.
		for (jf = 2, jc = 1; jc < nc - 1; jc++, jf += 2) {
			for (iif = 2, ic = 1; ic < nc - 1; ic++, iif += 2) {
				uc[ic][jc] = 0.25
						* uf[iif][jf]
						+ 0.125
						* (uf[iif + 1][jf] + uf[iif - 1][jf] + uf[iif][jf + 1] + uf[iif][jf - 1])
						+ 0.0625
						* (uf[iif + 1][jf + 1] + uf[iif - 1][jf - 1]
								+ uf[iif - 1][jf + 1] + uf[iif + 1][jf - 1]);

			}
		}

		//
		// edges
		for (jf = 0, jc = 0; jc < nc; jc++, jf += 2) { // Boundary points.
			uc[jc][0] = interpolate8PointWithCyclic(uf, jf, 0, nf);
			uc[jc][nc - 1] = interpolate8PointWithCyclic(uf, jf, nf - 1, nf);
			// top and bottom edges
			uc[0][jc] = interpolate8PointWithCyclic(uf, 0, jf, nf);
			uc[nc - 1][jc] = interpolate8PointWithCyclic(uf, nf - 1, jf, nf);
		}
	}

	/**
	 * @param uf
	 * @param iif
	 * @param jf
	 * @param nf
	 * @return the 8 point interpolation applying cyclic borders
	 */
	private static final double interpolate8PointWithCyclic(double[][] uf,
			int iif, int jf, int nf) {
		int ifPlus = (iif + 1 < nf) ? (iif + 1) : 0;
		int ifMinus = (iif - 1 < 0) ? (nf - 1) : (iif - 1);
		int jfPlus = (jf + 1 < nf) ? (jf + 1) : 0;
		int jfMinus = (jf - 1 < 0) ? (nf - 1) : (jf - 1);
		return 0.25
				* uf[iif][jf]
				+ 0.125
				* (uf[ifPlus][jf] + uf[ifMinus][jf] + uf[iif][jfPlus] + uf[iif][jfMinus])
				+ 0.0625
				* (uf[ifPlus][jfPlus] + uf[ifMinus][jfMinus]
						+ uf[ifMinus][jfPlus] + uf[ifPlus][jfMinus]);
	}

	/**
	 * Restricts the data in matrix uf to a grid one order coarser.
	 * 
	 * @param uc
	 *            coarser grid
	 * @param uf
	 *            finner grid
	 */
	public static void restrictCyclic4Point(double[][] uc, double[][] uf) {
		int nf = uf.length;
		int nc = uc.length;

		int ic, iif, jc, jf;
		// Interior points.
		for (jf = 2, jc = 1; jc < nc - 1; jc++, jf += 2) {
			for (iif = 2, ic = 1; ic < nc - 1; ic++, iif += 2) {
				uc[ic][jc] = 0.5
						* uf[iif][jf]
						+ 0.125
						* (uf[iif + 1][jf] + uf[iif - 1][jf] + uf[iif][jf + 1] + uf[iif][jf - 1]);

			}
		}

		//
		// edges
		for (jf = 0, jc = 0; jc < nc; jc++, jf += 2) { // Boundary points.
			uc[jc][0] = interpolate4PointWithCyclic(uf, jf, 0, nf);
			uc[jc][nc - 1] = interpolate4PointWithCyclic(uf, jf, nf - 1, nf);
			// top and bottom edges
			uc[0][jc] = interpolate4PointWithCyclic(uf, 0, jf, nf);
			uc[nc - 1][jc] = interpolate4PointWithCyclic(uf, nf - 1, jf, nf);
		}
	}

	/**
	 * @param uf
	 * @param iif
	 * @param jf
	 * @param nf
	 * @return the 8 point interpolation applying cyclic borders
	 */
	private static final double interpolate4PointWithCyclic(double[][] uf,
			int iif, int jf, int nf) {
		int ifPlus = (iif + 1 < nf) ? (iif + 1) : 0;
		int ifMinus = (iif - 1 < 0) ? (nf - 1) : (iif - 1);
		int jfPlus = (jf + 1 < nf) ? (jf + 1) : 0;
		int jfMinus = (jf - 1 < 0) ? (nf - 1) : (jf - 1);
		return 0.5
				* uf[iif][jf]
				+ 0.125
				* (uf[ifPlus][jf] + uf[ifMinus][jf] + uf[iif][jfPlus] + uf[iif][jfMinus]);
	}

	/**
	 * Restricts the data in matrix uf to a grid one order coarser by copying
	 * the same location points (no neighborhood averaging)
	 * 
	 * @param uc
	 *            coarser grid
	 * @param uf
	 *            finner grid
	 */
	public static void restrictByDirectCopy(double[][] uc, double[][] uf) {
		int nc = uc.length;

		int ic, iif, jc, jf;
		// Interior points.
		for (jf = 0, jc = 0; jc < nc; jc++, jf += 2) {
			for (iif = 0, ic = 0; ic < nc; ic++, iif += 2) {
				uc[ic][jc] = uf[iif][jf];
			}
		}
	}

	public static void restrictByExpanding(boolean[][] uc, boolean[][] uf) {
		int nc = uc.length;
		int nf = uf.length;

		int ic, iif, jc, jf;
		// Interior points.
		for (jf = 2, jc = 1; jc < nc - 1; jc++, jf += 2) {
			for (iif = 2, ic = 1; ic < nc - 1; ic++, iif += 2) {
				uc[ic][jc] = uf[iif][jf] | uf[iif + 1][jf] | uf[iif - 1][jf]
						| uf[iif][jf + 1] | uf[iif][jf - 1];
			}
		}

		//
		// edges
		for (jf = 0, jc = 0; jc < nc; jc++, jf += 2) { // Boundary points.
			uc[jc][0] = expand4PointWithCyclic(uf, jf, 0, nf);
			uc[jc][nc - 1] = expand4PointWithCyclic(uf, jf, nf - 1, nf);
			// top and bottom edges
			uc[0][jc] = expand4PointWithCyclic(uf, 0, jf, nf);
			uc[nc - 1][jc] = expand4PointWithCyclic(uf, nf - 1, jf, nf);
		}
	}

	/**
	 * @param uf
	 * @param iif
	 * @param jf
	 * @param nf
	 * @return the 8 point interpolation applying cyclic borders
	 */
	private static final boolean expand4PointWithCyclic(boolean[][] uf,
			int iif, int jf, int nf) {
		int ifPlus = (iif + 1 < nf) ? (iif + 1) : 0;
		int ifMinus = (iif - 1 < 0) ? (nf - 1) : (iif - 1);
		int jfPlus = (jf + 1 < nf) ? (jf + 1) : 0;
		int jfMinus = (jf - 1 < 0) ? (nf - 1) : (jf - 1);
		return uf[iif][jf] | uf[ifPlus][jf] | uf[ifMinus][jf] | uf[iif][jfPlus]
				| uf[iif][jfMinus];
	}

	/**
	 * Restricts the data in matrix uf to a grid one order coarser.
	 * 
	 * @param uc
	 *            coarser grid
	 * @param uf
	 *            finner grid
	 * @param lsf
	 *            level set distance function at finner grid
	 */
	public static void restrictOnlyInside(double[][] uc, double[][] uf,
			double[][] lsf) {
		int nf = uf.length;
		int nc = uc.length;

		int ic, iif, jc, jf;
		// Interior points.
		for (jf = 2, jc = 1; jc < nc - 1; jc++, jf += 2) {
			for (iif = 2, ic = 1; ic < nc - 1; ic++, iif += 2) {
				double v = uf[iif][jf];
				if (lsf[iif][jf] < 0) {
					uc[ic][jc] = 0.5
							* v
							+ 0.125
							* (valueIfInside(v, uf, lsf, iif + 1, jf)
									+ valueIfInside(v, uf, lsf, iif - 1, jf)
									+ valueIfInside(v, uf, lsf, iif, jf + 1) + valueIfInside(
									v, uf, lsf, iif, jf - 1));
				} else
					uc[ic][jc] = v;
			}
		}
		//
		// edges
		for (jf = 0, jc = 0; jc < nc; jc++, jf += 2) { // Boundary points.
			uc[jc][0] = interpolate8PointWithCyclicIfInside(uf, lsf, jf, 0, nf);
			uc[jc][nc - 1] = interpolate8PointWithCyclicIfInside(uf, lsf, jf,
					nf - 1, nf);
			// top and bottom edges
			uc[0][jc] = interpolate8PointWithCyclicIfInside(uf, lsf, 0, jf, nf);
			uc[nc - 1][jc] = interpolate8PointWithCyclicIfInside(uf, lsf,
					nf - 1, jf, nf);
		}
	}

	/**
	 * @param v
	 * @param u
	 * @param ls
	 * @param i
	 * @param j
	 * @return the value of u if inside (ls < 0); v otherwise
	 */
	private static final double valueIfInside(double v, double[][] u,
			double[][] ls, int i, int j) {
		return (ls[i][j] < 0 ? u[i][j] : v);
	}

	/**
	 * @param uf
	 * @param iif
	 * @param jf
	 * @param nf
	 * @return the 8 point interpolation applying cyclic borders
	 */
	private static final double interpolate8PointWithCyclicIfInside(
			double[][] uf, double[][] ls, int iif, int jf, int nf) {
		double v = uf[iif][jf];
		if (ls[iif][jf] >= 0)
			return v;
		int ifPlus = (iif + 1 < nf) ? (iif + 1) : 0;
		int ifMinus = (iif - 1 < 0) ? (nf - 1) : (iif - 1);
		int jfPlus = (jf + 1 < nf) ? (jf + 1) : 0;
		int jfMinus = (jf - 1 < 0) ? (nf - 1) : (jf - 1);
		return 0.5
				* v
				+ 0.125
				* (valueIfInside(v, uf, ls, ifPlus, jf)
						+ valueIfInside(v, uf, ls, ifMinus, jf)
						+ valueIfInside(v, uf, ls, iif, jfPlus) + valueIfInside(
						v, uf, ls, iif, jfMinus));
	}

	/**
	 * Restricts the data in matrix uf to a grid one order coarser.
	 * 
	 * @param uc
	 *            coarser grid
	 * @param uf
	 *            finner grid
	 */
	public static void restrictAllBordersCyclicWithTempPaddedMatrix8Point(
			double[][] uc, double[][] uf) {
		int nf = uf.length;
		int nc = uc.length;

		// create new matrix with borders to enforce cyclic boundarie conditions
		double[][] ufBorders = new double[nf + 2][nf + 2];
		int i, j;
		for (i = 0; i < nf; i++) {
			ufBorders[0][i + 1] = uf[nf - 1][i];
			ufBorders[nf + 1][i + 1] = uf[0][i];
			ufBorders[i + 1][0] = uf[i][nf - 1];
			ufBorders[i + 1][nf + 1] = uf[i][0];
			for (j = 0; j < nf; j++) {
				ufBorders[i + 1][j + 1] = uf[i][j];
			}
		}
		ufBorders[0][0] = ufBorders[nf][nf];
		ufBorders[0][nf + 1] = ufBorders[nf][1];
		ufBorders[nf + 1][0] = ufBorders[1][nf];
		ufBorders[nf + 1][nf + 1] = ufBorders[1][1];

		int ic, iif, jc, jf;
		// Interior points.
		for (jf = 0, jc = 0; jc < nc; jc++, jf += 2) {
			for (iif = 0, ic = 0; ic < nc; ic++, iif += 2) {
				uc[ic][jc] = 0.25
						* ufBorders[iif + 1][jf + 1]
						+ 0.125
						* (ufBorders[iif + 2][jf + 1] + ufBorders[iif][jf + 1]
								+ ufBorders[iif + 1][jf + 2] + ufBorders[iif + 1][jf])
						+ 0.0625
						* (ufBorders[iif + 2][jf + 2] + ufBorders[iif][jf]
								+ ufBorders[iif][jf + 2] + ufBorders[iif + 2][jf]);
			}
		}
	}

	/**
	 * Interpolares the data in matrix uc to a grid one order finner.
	 * 
	 * @param uf
	 *            finner grid
	 * @param uc
	 *            coarser grid
	 */
	public static void interpolateAllBordersCyclicTo(double[][] uf,
			double[][] uc) {
		int nc = uc.length;
		int nf = uf.length;

		int ic, iif, jc, jf;
		for (jc = 0, jf = 0; jc < nc; jc++, jf += 2)
			// Do elements that are copies.
			for (ic = 0; ic < nc; ic++)
				uf[2 * ic][jf] = uc[ic][jc];
		for (jf = 0; jf < nf; jf += 2)
			// Do odd-numbered columns, interpolating vertically.
			for (iif = 1; iif < nf; iif += 2)
				uf[iif][jf] = 0.5 * (uf[iif + 1][jf] + uf[iif - 1][jf]);
		for (jf = 1; jf < nf; jf += 2)
			// Do even-numbered columns, interpolating horizontally.
			for (iif = 0; iif < nf; iif++)
				uf[iif][jf] = 0.5 * (uf[iif][jf + 1] + uf[iif][jf - 1]);
	}

	/**
	 * subtract all entries of matrix b to a and place the result in c
	 * 
	 * @param a
	 * @param b
	 * @param c
	 */
	public static void subtractTo(double[][] a, double[][] b, double[][] c) {
		int n = a.length;
		// TODO: remove after testing phase
		if ((a[0].length != n) || (b.length != n) || (b[0].length != n))
			throw new RuntimeException(
					"matrix a and b must be square and have same size");
		int i, j;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				c[i][j] = a[i][j] - b[i][j];
	}

	/**
	 * add all entries of matrix b to a
	 * 
	 * @param a
	 * @param b
	 */
	public static void addTo(double[][] a, double[][] b) {
		int n = a.length;
		// TODO: remove after testing phase
		if ((a[0].length != n) || (b.length != n) || (b[0].length != n))
			throw new RuntimeException(
					"matrix a and b must be square and have same size");
		int i, j;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				a[i][j] += b[i][j];
	}

	/**
	 * Set all the values in square matrix to 0
	 * 
	 * @param u
	 * @param n
	 *            size of side of matrix
	 */
	public static void resetMatrix(double[][] u, int n) {
		int i, j;
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
				u[i][j] = 0;
	}

}
