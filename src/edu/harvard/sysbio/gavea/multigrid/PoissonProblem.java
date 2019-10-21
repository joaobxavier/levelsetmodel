package edu.harvard.sysbio.gavea.multigrid;

import edu.harvard.sysbio.gavea.utils.ScalarMath;

/**
 * The mother class of all poisson problems
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 12, 2007
 */
public abstract class PoissonProblem implements Cloneable {
	protected int _iSize;

	private double _size;

	protected double _dx;

	protected double _dx2;

	protected double residualGoal;

	private double _coarsestGridToUseInMultigridSolver = 3; // default 3x3

	private double _finnesttGridToUseInMultigridSolver;

	protected double[][] temp; // temp matrix used in multigrid

	public double res, trerr;

	protected double[][] s; // temp matrix used in multigrid

	protected double[][] rhs; // temp matrix used in multigrid

	protected PoissonProblem pCoarser;

	public PoissonProblem(double size) {
		_size = size;
	}

	public boolean pointIsInsideDomain(int i, int j) {
		return true;
	}

	/**
	 * Set the residual goal of the problem
	 * 
	 * @param rg
	 */
	public void setResidualGoal(double rg) {
		residualGoal = rg;
	}

	/**
	 * This method sets the ,atrix size and, IMPORTANTLY, initializes all the
	 * memory demanding data structures such as the s and rhs matrices. It also
	 * initializes the pCoarser attribute, triggering a recursive creation of
	 * coarser grid sizes until the coarses grid is reached. Classes that extend
	 * PoissonProblem (e.g. PressureField) must take care to implement all the
	 * necessary initializations
	 * 
	 * @param n
	 */
	protected void setSize(int n) {
		// declare finnest grid if this is first call
		if (_finnesttGridToUseInMultigridSolver == 0)
			_finnesttGridToUseInMultigridSolver = n;
		_iSize = n;
		temp = new double[n][n];
		s = new double[n][n];
		rhs = new double[n][n];
		_dx = _size / _iSize;
		_dx2 = ScalarMath.square(_dx);
		// create the corser instance if this is not the coarset instance
		if (!isCoarset()) {
			try {
				pCoarser = (PoissonProblem) this.clone();
				int coarserSize = MultiGridOperations.coarserSize(_iSize);
				pCoarser.setSize(coarserSize);
			} catch (CloneNotSupportedException e) {
				throw new RuntimeException(e);
			}
		}
	}

	/**
	 * @return the size of grid
	 */
	public int getSize() {
		return _iSize;
	}

	/**
	 * @return the reference to the solution
	 */
	public double[][] getSolution() {
		return s;
	}

	/**
	 * @return the initial right-hand-side of the poisson problem
	 */
	public abstract void setOriginalRightHandSide();

	/**
	 * Implements the L-operator for the poisson problem (may be non-linear)
	 * 
	 * @param i
	 * @param j
	 * @return the result of applying L-operator at location i,j
	 */
	public abstract double getLOperator(int i, int j);

	/**
	 * Implements the L-operator for the poisson problem (may be non-linear)
	 * with ghost node if needed
	 * 
	 * @param i
	 * @param j
	 * @return the result of applying L-operator at location i,j with ghostnode
	 */
	public abstract double getLOperatorWithGhostNode(int i, int j);

	/**
	 * The derivative of the L-operator
	 * 
	 * @param i
	 * @param j
	 * @return
	 */
	public abstract double getLOperatorDerivative(int i, int j);

	/**
	 * The derivative of the L-operator
	 * 
	 * @param i
	 * @param j
	 * @return
	 */
	public abstract double getLOperatorDerivativeWithGhostNode(int i, int j);

	/**
	 * create a new matrix with the initial solution
	 */
	public abstract void createInitialSolution();

	/**
	 * To assess if problem is defined on the coarsest grid. Throws a run-time
	 * exception if it detects that problem is not defined in a 2^n+1 side
	 * square grid.
	 * 
	 * @return true if problem is defined in the coarset grid (3x3) false
	 *         otherwise.
	 */
	public boolean isCoarset() {
		return _iSize == _coarsestGridToUseInMultigridSolver;
	}

	/**
	 * To assess if problem is defined on the finnest grid.
	 * 
	 * @return true if problem is defined in the finnest grid
	 */
	public boolean isFinnest() {
		return _iSize == _finnesttGridToUseInMultigridSolver;
	}

	/**
	 * Update the coraser version of this problem (restrict the problem)
	 * 
	 * @return the same problem, restricted to a coarse grid
	 */
	public abstract void restrict();

	/**
	 * @return the matrix interpolated to a finner grid
	 */
	public abstract void interpolateTo(double[][] uf, double[][] uc);

	/**
	 * @return the matrix restricted to a coarser grid
	 */
	public abstract void restrictTo(double[][] uc, double[][] uf);

	/**
	 * Create the L-Operator matrix
	 * 
	 * @param u
	 * @param rhs
	 * @return the L-Operator matrix
	 */
	public void setToLOperatorMatrix(double[][] a) {
		for (int i = 0; i < _iSize; i++) {
			for (int j = 0; j < _iSize; j++) {
				if (this.pointIsInsideDomain(i, j))
					a[i][j] = this.getLOperator(i, j);
			}
		}
	}

	/**
	 * Compute all the downward V operations and return the trerr
	 * 
	 * @param temp
	 * @param rhs
	 * @return the trerr
	 */
	public double downwardVOperations() {
		double tau;
		double trerr = 0;
		double count = 0.0;
		for (int i = 0; i < _iSize; i++) {
			for (int j = 0; j < _iSize; j++) {
				if (this.pointIsInsideDomain(i, j)) {
					// truncation error
					tau = getLOperator(i, j) - temp[i][j];
					rhs[i][j] += tau;
					trerr += ScalarMath.square(tau);
					count += 1.0;
				}
			}
		}
		return Math.sqrt(trerr / count);
	}

	/**
	 * @param u
	 * @return eucledian norm of matrix u
	 */
	public double anorm2Temp() {
		int i, j;
		double sum = 0.0;
		double count = 0.0;
		for (j = 0; j < _iSize; j++)
			for (i = 0; i < _iSize; i++)
				if (this.pointIsInsideDomain(i, j)) {
					sum += ScalarMath.square(temp[i][j]);
					count += 1.0;
				}
		return Math.sqrt(sum / count);
	}

	/**
	 * Stores the present grid size as the finest grid size. Initializes the
	 * matrices at coarser levels and calls solve with multigrid from the
	 * MultigridSolver
	 * 
	 * @return the solution
	 */
	public double[][] solveWithMultigrid(double alpha, int maxCycles) {
		// initialize multigrid copies
		PoissonProblem p = this;
		while (!p.isCoarset()) {
			p.restrict();
			p = p.pCoarser;
		}
		MultigridSolver.maxVCylces = maxCycles;
		MultigridSolver.solveWithMultigrid(this, residualGoal, alpha);
		return s;
	}

	/**
	 * @return the residual goeal
	 */
	public double getResidualGoal() {
		return residualGoal;
	}

	/**
	 * @param g
	 *            grid number to set
	 */
	public void setCoarsestGridToUseInMultigridSolver(int g) {
		if (!MultiGridOperations.validMultigridSize(g))
			throw new RuntimeException("trying to set coarsest grid to " + g);
		_coarsestGridToUseInMultigridSolver = g;
	}
}
