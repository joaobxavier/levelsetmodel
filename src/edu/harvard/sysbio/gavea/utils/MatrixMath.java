package edu.harvard.sysbio.gavea.utils;

/**
 * Scalar math operations
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 12, 2007
 */
public abstract class MatrixMath {
	/**
	 * @param aI
	 * @param aJ
	 * @return the maximum euclidean distance
	 */
	public static final double maximumNorm(double[][] aI, double[][] aJ) {
		int n = aI.length;
		int i, j;
		double maxNorm = 0.0;
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++) {
				// maxNorm = Math.max((Math.abs(aI[i][j]) + Math.abs(aJ[i][j])),
				// maxNorm);
				maxNorm = Math.max(Math.sqrt(ScalarMath.square(aI[i][j])
						+ ScalarMath.square(aJ[i][j])), maxNorm);
			}
		return maxNorm;
	}

	/**
	 * Copy all values from source to destination
	 * 
	 * @param source
	 * @param destination
	 */
	public static final void copyTo(double[][] source, double[][] destination) {
		int n = source.length;
		int i, j;
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
				destination[i][j] = source[i][j];
	}

	/**
	 * performs destination = aMultiplier*a + bMultiplier*b
	 * 
	 * @param a
	 * @param aMultiplier
	 * @param b
	 * @param bMultiplier
	 * @param destination
	 */
	public static final void linearCombination(double[][] a,
			double aMultiplier, double[][] b, double bMultiplier,
			double[][] destination) {
		int n = a.length;
		int i, j;
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
				destination[i][j] = aMultiplier * a[i][j] + bMultiplier
						* b[i][j];
	}

	/**
	 * Use to detect if an outside point is located at the border
	 * 
	 * @param i
	 * @param j
	 * @param iPlus
	 * @param iMinus
	 * @param jPlus
	 * @param jMinus
	 * @return true if any of the neighboring points are inside
	 */
	public static final boolean anyNeighborIsInside(double[][] levelSet, int i,
			int j, int iPlus, int iMinus, int jPlus, int jMinus) {
		return (levelSet[iPlus][j] <= 0) | (levelSet[iMinus][j] <= 0)
				| (levelSet[i][jPlus] <= 0) | (levelSet[i][jMinus] <= 0);
	}

	/**
	 * Use to detect if an inside point is located at the border
	 * 
	 * @param i
	 * @param j
	 * @param iPlus
	 * @param iMinus
	 * @param jPlus
	 * @param jMinus
	 * @return true if any of the neighboring points are outside
	 */
	public static final boolean anyNeighborIsOutside(double[][] levelSet,
			int i, int j, int iPlus, int iMinus, int jPlus, int jMinus) {
		return (levelSet[iPlus][j] > 0) | (levelSet[iMinus][j] > 0)
				| (levelSet[i][jPlus] > 0) | (levelSet[i][jMinus] > 0);
	}

	/**
	 * Returns true is a point is in the border, false otherwise
	 * 
	 * @param levelSet
	 * @param i
	 * @param j
	 * @param iPlus
	 * @param iMinus
	 * @param jPlus
	 * @param jMinus
	 * @return true is a point is in the border, false otherwise
	 */
	public static final boolean isBorderPoint(double[][] levelSet, int i,
			int j, int iPlus, int iMinus, int jPlus, int jMinus) {
		if (levelSet[i][j] * levelSet[iPlus][j] <= 0)
			return true;
		if (levelSet[i][j] * levelSet[iMinus][j] <= 0)
			return true;
		if (levelSet[i][j] * levelSet[i][jPlus] <= 0)
			return true;
		if (levelSet[i][j] * levelSet[i][jMinus] <= 0)
			return true;
		return false;
	}

	/**
	 * This method counts the number of imidiate neighbors that belong to the
	 * biofilm. To be used in the initiation of the distance matrix.
	 * 
	 * @param m
	 * @param i
	 * @param j
	 * @param iPlus
	 * @param iMinus
	 * @param jPlus
	 * @param jMinus
	 * @return the number of neighbors that have a value greater than 0
	 */
	public static final int numberOfPositiveValueNeighbors(double[][] m, int i,
			int j, int iPlus, int iMinus, int jPlus, int jMinus) {
		int n = 0;
		n += (m[iPlus][j] > 0 ? 1 : 0);
		n += (m[iMinus][j] > 0 ? 1 : 0);
		n += (m[i][jPlus] > 0 ? 1 : 0);
		n += (m[i][jMinus] > 0 ? 1 : 0);
		return n;
	}

	/**
	 * This method counts the number of imidiate neighbors that do not belong to
	 * the biofilm. To be used in the initiation of the distance matrix.
	 * 
	 * @param m
	 * @param i
	 * @param j
	 * @param iPlus
	 * @param iMinus
	 * @param jPlus
	 * @param jMinus
	 * @return the number of neighbors that have a value equal to 0
	 */
	public static final int numberOfZeroValueNeighbors(double[][] m, int i,
			int j, int iPlus, int iMinus, int jPlus, int jMinus) {
		int n = 0;
		n += (m[iPlus][j] == 0 ? 1 : 0);
		n += (m[iMinus][j] == 0 ? 1 : 0);
		n += (m[i][jPlus] == 0 ? 1 : 0);
		n += (m[i][jMinus] == 0 ? 1 : 0);
		return n;
	}

	/**
	 * @param m
	 * @return hte maximum value in the matrix
	 */
	public static double maxValue(double[][] m) {
		int n = m.length;
		double maxVal = Double.NEGATIVE_INFINITY;
		int i, j;
		for (j = 0; j < n; j++)
			for (i = 0; i < n; i++)
				maxVal = Math.max(m[i][j], maxVal);
		return maxVal;
	}
}
