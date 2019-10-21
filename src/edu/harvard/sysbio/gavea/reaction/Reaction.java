package edu.harvard.sysbio.gavea.reaction;

/**
 * The abstract class defining all reactions
 * 
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Feb 4, 2007
 */
public abstract class Reaction {
	/**
	 * @return the substrate rate (positive if produced, negative if consumed)
	 */
	public abstract double getSubstrateRate(double s, double b);

	/**
	 * @return the substrate rate derivative
	 */
	public abstract double getSubstrateRateDerivative(double s, double b);

	/**
	 * @return the specific growth ratte of biomass
	 */
	public abstract double getSpecificGrowthRate(double s);

	/**
	 * Creates a new matrix and fills it with the substrate rate
	 * 
	 * @param s
	 *            the substrate concentration
	 * @param b
	 *            the biomass concentration
	 * @return the substrate rate
	 */
	public double[][] computeSubstrateRate(double[][] s, double[][] b) {
		int n = s.length;
		double[][] r = new double[n][n];
		for (int i = 0; i < n; i++)
			for (int j = 0; j < n; j++)
				r[i][j] = getSubstrateRate(s[i][j], b[i][j]);
		return r;
	}
}
