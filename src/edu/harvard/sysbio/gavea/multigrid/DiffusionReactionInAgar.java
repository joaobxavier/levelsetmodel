package edu.harvard.sysbio.gavea.multigrid;

import edu.harvard.sysbio.gavea.reaction.Reaction;
import edu.harvard.sysbio.gavea.utils.MatrixMath;
import edu.harvard.sysbio.gavea.utils.ScalarMath;

/**
 * Implements the diffusion-reaction in a 2D plane with the nutrient comming
 * from above the plane. This class initializes the instance of biomass
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 12, 2007
 */
public class DiffusionReactionInAgar extends PoissonProblem {
	private Reaction _reaction;

	private double _diffusivity;

	private double _dz;

	private double _cBulk;

	private double _density;

	private double[][] _biomass;

	private double aux1;

	private double aux2;

	private double aux3;

	private double aux4;

	private double aux5;

	public DiffusionReactionInAgar(double diffusivity, double size, double dz,
			Reaction r, double cBulk, double density, int nGrid) {
		super(size);
		_diffusivity = diffusivity;
		_dz = dz;
		_reaction = r;
		_cBulk = cBulk;
		_density = density;
		// auxiliary variable (to speedup computation
		aux3 = _diffusivity / ScalarMath.square(_dz);
		// initializes the instance to biomass as well as all the coarse grids
		setSize(nGrid);
	}

	@Override
	protected void setSize(int n) {
		super.setSize(n);
		_biomass = new double[n][n];
		// auxiliary variable (to speedup computation
		aux2 = -4.0 * _diffusivity / ScalarMath.square(_dx) - _diffusivity
				/ ScalarMath.square(_dz);
		aux5 = _diffusivity / ScalarMath.square(_dx);
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
		// TODO Auto-generated method stub
		return super.clone();
	}

	@Override
	public void setOriginalRightHandSide() {
		MultiGridOperations.resetMatrix(rhs, _iSize);
	}

	@Override
	public void createInitialSolution() {
		MultiGridOperations.resetMatrix(s, _iSize);
	}

	@Override
	public double getLOperatorDerivative(int i, int j) {
		double s = this.s[i][j];
		double b = _biomass[i][j];
		double dr = _reaction.getSubstrateRateDerivative(s, b);
		return aux2 + dr;
	}

	@Override
	public double getLOperator(int i, int j) {
		double s = this.s[i][j];
		double b = _biomass[i][j];
		double soluteConvulution = soluteConvolution(this.s, i, j);
		// compute the rate
		double r = _reaction.getSubstrateRate(s, b);
		return aux5 * soluteConvulution + r + aux3 * (_cBulk - s);
	}

	@Override
	public double getLOperatorWithGhostNode(int i, int j) {
		// Ghost node is not needed for the solutes
		return getLOperator(i, j);
	}

	@Override
	public double getLOperatorDerivativeWithGhostNode(int i, int j) {
		return getLOperatorDerivative(i, j);
	}

	/**
	 * @param u
	 * @param i
	 * @param j
	 * @return
	 */
	private double soluteConvolution(double[][] u, int i, int j) {
		return u[i > 0 ? i - 1 : _iSize - 1][j]
				+ u[i < _iSize - 1 ? i + 1 : 0][j]
				+ u[i][j > 0 ? j - 1 : _iSize - 1]
				+ u[i][j < _iSize - 1 ? j + 1 : 0] - 4 * u[i][j];
	}

	@Override
	public void interpolateTo(double[][] uf, double[][] uc) {
		MultiGridOperations.interpolateAllBordersCyclicTo(uf, uc);
	}

	@Override
	public void restrictTo(double[][] uc, double[][] uf) {
		MultiGridOperations.restrictCyclic8Point(uc, uf);
	}

	@Override
	public void restrict() {
		DiffusionReactionInAgar pc = (DiffusionReactionInAgar) pCoarser;
		// restrict the the biomass
		restrictTo(pc._biomass, _biomass);
	}

	/**
	 * Computes the specific biomass growth rate and places the result in the
	 * matrix a.
	 * 
	 * @param a
	 *            matrix to place the result in
	 */
	public void updateSpecificGrowthRate(double[][] a) {
		int i, j;
		for (i = 0; i < _iSize; i++) {
			for (j = 0; j < _iSize; j++) {
				if (_biomass[i][j] > 0)
					// only assign values inside biofilm
					a[i][j] = _reaction.getSpecificGrowthRate(s[i][j]);
				else
					a[i][j] = 0;
			}
		}
	}

	
	/**
	 * Set values of biomass to the maximum biomass density if inside (distance <
	 * 0) or to 0 otherwise
	 * 
	 * @param setDistanceMatrix
	 */
	public void updateValuesFromDistanceMatrix(double[][] distanceMatrix) {
		for (int i = 0; i < _iSize; i++)
			for (int j = 0; j < _iSize; j++) {
				int iPlus = ScalarMath.cyclicUpwards(i + 1, _iSize);
				int iMinus = ScalarMath.cyclicDownwards(i - 1, _iSize);
				int jPlus = ScalarMath.cyclicUpwards(j + 1, _iSize);
				int jMinus = ScalarMath.cyclicDownwards(j - 1, _iSize);
				if (MatrixMath.isBorderPoint(distanceMatrix, i, j, iPlus,
						iMinus, jPlus, jMinus))
				// compute the fraction of volume that is filled with
				// biomass
				{
					if (distanceMatrix[i][j] < 0) {
						// otherwise compute ghostnode
						double thetaIPlus = ghostNodeThetaInsidePoints(
								distanceMatrix, i, j, iPlus, j);
						double thetaIMinus = ghostNodeThetaInsidePoints(
								distanceMatrix, i, j, iMinus, j);
						double thetaJPlus = ghostNodeThetaInsidePoints(
								distanceMatrix, i, j, i, jPlus);
						double thetaJMinus = ghostNodeThetaInsidePoints(
								distanceMatrix, i, j, i, jMinus);
						double f = (thetaIPlus + thetaIMinus)
								* (thetaJPlus + thetaJMinus);
						// different fraction for points outside and inside
						_biomass[i][j] = _density * f;
					} else {
						// otherwise compute ghostnode
						double thetaIPlus = ghostNodeThetaOutsidePoints(
								distanceMatrix, i, j, iPlus, j);
						double thetaIMinus = ghostNodeThetaOutsidePoints(
								distanceMatrix, i, j, iMinus, j);
						double thetaJPlus = ghostNodeThetaOutsidePoints(
								distanceMatrix, i, j, i, jPlus);
						double thetaJMinus = ghostNodeThetaOutsidePoints(
								distanceMatrix, i, j, i, jMinus);
						double f = (thetaIPlus + thetaIMinus)
								* (thetaJPlus + thetaJMinus);
						// different fraction for points outside and inside
						_biomass[i][j] = _density * (1 - f);
					}
				} else if (distanceMatrix[i][j] < 0) {
					_biomass[i][j] = _density;
				} else
					_biomass[i][j] = 0;
			}
	}

	/**
	 * @param levelSet
	 * @param i
	 * @param j
	 * @param iNeib
	 * @param jNeib
	 * @return the theta value
	 */
	private final double ghostNodeThetaInsidePoints(double[][] distanceMatrix,
			int i, int j, int iNeib, int jNeib) {
		if ((distanceMatrix[iNeib][jNeib] > 0) & (distanceMatrix[i][j] < _dx)) {
			double val = distanceMatrix[i][j]
					/ (distanceMatrix[i][j] - distanceMatrix[iNeib][jNeib]);
			return (val > 0.5 ? 0.5 : val);
		}
		return 0.5;
	}

	/**
	 * @param levelSet
	 * @param i
	 * @param j
	 * @param iNeib
	 * @param jNeib
	 * @return the theta value
	 */
	private final double ghostNodeThetaOutsidePoints(double[][] distanceMatrix,
			int i, int j, int iNeib, int jNeib) {
		if ((distanceMatrix[iNeib][jNeib] < 0)
				& ((-distanceMatrix[i][j]) < _dx)) {
			double val = distanceMatrix[i][j]
					/ (distanceMatrix[i][j] - distanceMatrix[iNeib][jNeib]);
			return (val > 0.5 ? 0.5 : val);
		}
		return 0.5;
	}

	/**
	 * @return the instance of biomass (there should be only one!)
	 */
	public double[][] getBiomass() {
		return _biomass;
	}
	
	/**
	 * Set values of biomass to the maximum biomass density if inside (distance <
	 * 0) or to 0 otherwise
	 * 
	 * @param setDistanceMatrix
	 */
	public void setBiomassValues(double[][] biomass) {
		for (int i = 0; i < _iSize; i++) {
			for (int j = 0; j < _iSize; j++) {
				// different fraction for points outside and inside
				_biomass[i][j] = biomass[i][j];
			}
		}
	}

}
