package edu.harvard.sysbio.gavea.multigrid;

import edu.harvard.sysbio.gavea.utils.MatrixMath;
import edu.harvard.sysbio.gavea.utils.ScalarMath;

/**
 * Implements the diffusion-reaction in a 2D plane with the nutrient comming
 * from above the plane
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 12, 2007
 */
public class PressureField extends PoissonProblem {
	private double[][] _specificGrowthRate;

	private double[][] _levelSetDistanceMatrix;

	private boolean[][] _isBorder;

	private double aux1;

	private double aux2;

	private double aux3;

	public PressureField(double size, int n) {
		super(size);
		setSize(n);
	}

	@Override
	protected void setSize(int n) {
		super.setSize(n);
		aux1 = 4.0 / ScalarMath.square(_dx);
		aux2 = 1.0 / ScalarMath.square(_dx);
		aux3 = 2.0 / ScalarMath.square(_dx);
		_levelSetDistanceMatrix = new double[n][n];
		_specificGrowthRate = new double[n][n];
		_isBorder = new boolean[n][n];
	}

	@Override
	public boolean pointIsInsideDomain(int i, int j) {
		// only points with biomass are inside domain
		// return (_levelSetDistanceMatrix[i][j] < 0);
		return !_isBorder[i][j];
	}

	@Override
	protected Object clone() throws CloneNotSupportedException {
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
		// return getLOperatorDerivativeNoGhostNode(i, j);
		return getLOperatorDerivativeWithGhostNode(i, j);
	}

	public double getLOperatorDerivativeNoGhostNode(int i, int j) {
		return aux1;
	}

	@Override
	public double getLOperator(int i, int j) {
		// return getLOperatorNoGhostNode(i, j);
		return getLOperatorWithGhostNode(i, j);
	}

	public double getLOperatorNoGhostNode(int i, int j) {
		double g = _specificGrowthRate[i][j];
		double soluteConvulution = soluteConvolution(this.s, i, j);
		return -(soluteConvulution * aux2 + g);
	}

	/**
	 * @param u
	 * @param i
	 * @param j
	 * @return the difference operator in usual (non-ghost node) condiitons
	 */
	private double soluteConvolution(double[][] u, int i, int j) {
		return u[i > 0 ? i - 1 : _iSize - 1][j]
				+ u[i < _iSize - 1 ? i + 1 : 0][j]
				+ u[i][j > 0 ? j - 1 : _iSize - 1]
				+ u[i][j < _iSize - 1 ? j + 1 : 0] - 4 * u[i][j];
	}

	@Override
	public double getLOperatorWithGhostNode(int i, int j) {
		int iPlus = ScalarMath.cyclicUpwards(i + 1, _iSize);
		int iMinus = ScalarMath.cyclicDownwards(i - 1, _iSize);
		int jPlus = ScalarMath.cyclicUpwards(j + 1, _iSize);
		int jMinus = ScalarMath.cyclicDownwards(j - 1, _iSize);
		// if point is not in border, use the regurlar function
		if (_levelSetDistanceMatrix[iPlus][j] < 0
				& _levelSetDistanceMatrix[iMinus][j] < 0
				& _levelSetDistanceMatrix[i][jPlus] < 0
				& _levelSetDistanceMatrix[i][jMinus] < 0)
			return getLOperatorNoGhostNode(i, j);
		// otherwise compute ghostnode
		double g = _specificGrowthRate[i][j];
		double thetaIPlus = ghostNodeTheta(i, j, iPlus, j);
		double thetaIMinus = ghostNodeTheta(i, j, iMinus, j);
		double thetaJPlus = ghostNodeTheta(i, j, i, jPlus);
		double thetaJMinus = ghostNodeTheta(i, j, i, jMinus);

		double d2sdi2 = ((s[iPlus][j] - s[i][j]) / thetaIPlus - (s[i][j] - s[iMinus][j])
				/ thetaIMinus)
				/ (thetaIPlus + thetaIMinus) * aux3;
		double d2sdj2 = ((s[i][jPlus] - s[i][j]) / thetaJPlus - (s[i][j] - s[i][jMinus])
				/ thetaJMinus)
				/ (thetaJPlus + thetaJMinus) * aux3;
		// TODO optimize this code
		double aux = -(d2sdi2 + d2sdj2 + g);

		return aux;
	}

	@Override
	public double getLOperatorDerivativeWithGhostNode(int i, int j) {
		int iPlus = ScalarMath.cyclicUpwards(i + 1, _iSize);
		int iMinus = ScalarMath.cyclicDownwards(i - 1, _iSize);
		int jPlus = ScalarMath.cyclicUpwards(j + 1, _iSize);
		int jMinus = ScalarMath.cyclicDownwards(j - 1, _iSize);
		// if point is not in border, use the regurlar function
		if (_levelSetDistanceMatrix[iPlus][j] < 0
				& _levelSetDistanceMatrix[iMinus][j] < 0
				& _levelSetDistanceMatrix[i][jPlus] < 0
				& _levelSetDistanceMatrix[i][jMinus] < 0)
			return getLOperatorDerivativeNoGhostNode(i, j);
		// otherwise compute ghostnode
		double thetaIPlus = ghostNodeTheta(i, j, iPlus, j);
		double thetaIMinus = ghostNodeTheta(i, j, iMinus, j);
		double thetaJPlus = ghostNodeTheta(i, j, i, jPlus);
		double thetaJMinus = ghostNodeTheta(i, j, i, jMinus);
		return (1.0 / thetaIPlus + 1.0 / thetaIMinus)
				/ (thetaIPlus + thetaIMinus) * aux3
				+ (1.0 / thetaJPlus + 1.0 / thetaJMinus)
				/ (thetaJPlus + thetaJMinus) * aux3;
	}

	/**
	 * @param p
	 * @param i
	 * @param j
	 * @param iNeib
	 * @param jNeib
	 * @return the theta value
	 */
	private final double ghostNodeTheta(int i, int j, int iNeib, int jNeib) {
		if ((_levelSetDistanceMatrix[iNeib][jNeib] > 0)
				& (_levelSetDistanceMatrix[i][j] < _dx)) {
			return _levelSetDistanceMatrix[i][j]
					/ (_levelSetDistanceMatrix[i][j] - _levelSetDistanceMatrix[iNeib][jNeib]);
			// TODO check which one is better
			// return -_levelSetDistanceMatrix[i][j] / _dx;
		}
		return 1.0;
	}

	@Override
	public void interpolateTo(double[][] uf, double[][] uc) {
		// MultiGridOperations.interpolateAllBordersCyclicTo(uf, uc);

		int nc = uc.length;
		int nf = uf.length;

		int ic, iif, jc, jf;
		for (jc = 0, jf = 0; jc < nc; jc++, jf += 2)
			// Do elements that are copies.
			for (ic = 0; ic < nc; ic++) {
				if (pointIsInsideDomain(2 * ic, jf))
					uf[2 * ic][jf] = uc[ic][jc];
				else
					uf[2 * ic][jf] = 0;
			}
		for (jf = 0; jf < nf; jf += 2)
			// Do odd-numbered columns, interpolating vertically.
			for (iif = 1; iif < nf; iif += 2)
				if (pointIsInsideDomain(iif, jf)) {
					uf[iif][jf] = 0.5 * (uf[iif + 1][jf] + uf[iif - 1][jf]);
					// uf[iif][jf] = average2(uf, iif + 1, jf, iif - 1, jf);
				} else
					uf[iif][jf] = 0;
		for (jf = 1; jf < nf; jf += 2)
			// Do even-numbered columns, interpolating horizontally.
			for (iif = 0; iif < nf; iif++)
				if (pointIsInsideDomain(iif, jf))
					uf[iif][jf] = 0.5 * (uf[iif][jf + 1] + uf[iif][jf - 1]);
				else
					uf[iif][jf] = 0;
		// uf[iif][jf] = average2(uf, iif, jf + 1, iif, jf - 1);
	}

	@Override
	public void restrictTo(double[][] uc, double[][] uf) {
		// MultiGridOperations.restrictOnlyInside(uc, uf,
		// _levelSetDistanceMatrix);
		MultiGridOperations.restrictCyclic4Point(uc, uf);
	}

	@Override
	public void restrict() {
		PressureField pRestricted = (PressureField) pCoarser;
		// restric the the specific growth rate
		restrictTo(pRestricted._specificGrowthRate, _specificGrowthRate);

		// restrict the isborder
		MultiGridOperations.restrictByExpanding(pRestricted._isBorder,
				_isBorder);
		// restrict the levelset
		restrictTo(pRestricted._levelSetDistanceMatrix, _levelSetDistanceMatrix);
	}

	/**
	 * TODO check if it is realy necessary to use ghost node here. If ghost node
	 * is not necessary then the instance to levelset does not have to be used
	 * here
	 * 
	 * @return the instance of the levelset matrix
	 */
	public double[][] getLevelSetDistanceMatrix() {
		return _levelSetDistanceMatrix;
	}

	/**
	 * @return the only instance of specific growth rate
	 */
	public double[][] getSpecificGrowthRate() {
		return _specificGrowthRate;
	}

	public double[][] getDistanceMatrix() {
		return _levelSetDistanceMatrix;
	}

	/**
	 * Copy the values to the instance of level set and update the _isBiomass
	 * 
	 * @param setDistanceMatrix
	 */
	public void copyValuesToDistanceMatrix(double[][] setDistanceMatrix) {
		MatrixMath.copyTo(setDistanceMatrix, _levelSetDistanceMatrix);
		updateIsBorder();
	}

	/**
	 * Must be ran each time before solving with multigrid
	 */
	public void updateIsBorder() {
		for (int i = 0; i < _iSize; i++)
			for (int j = 0; j < _iSize; j++)
				_isBorder[i][j] = (_levelSetDistanceMatrix[i][j] > 0);
	}

	/**
	 * Set a single value in the instance of level set and update the same
	 * location in the _isBiomass
	 * 
	 * @param setDistanceMatrix
	 */
	public void setValueInDistanceMatrix(double v, int i, int j) {
		_levelSetDistanceMatrix[i][j] = v;
		_isBorder[i][j] = (v > 0);
	}

}