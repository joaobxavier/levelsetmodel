package edu.harvard.sysbio.gavea.multigrid;

import edu.harvard.sysbio.gavea.utils.ScalarMath;

/**
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 12, 2007
 */
public abstract class MultigridSolver {
	public static int maxVCylces = 1;

	public static double ALPHA = 0.3;

	public static int maxIterations = 1000;

	// number of pre-processing steps
	public static int NPREPROCESSING = 5;

	// number of post-processing steps
	public static int NPOSTPROCESSING = 5;

	// number of steps on coarsest grid
	private static int NCOARSEST = 5000;

	/**
	 * Solve the a poisson problem (independent of its resolution)
	 * 
	 * @param biomass
	 * @return the solved solute matrix
	 */
	public static void solveWithMultigrid(PoissonProblem p,
			double residualGoal, double alpha) {
		ALPHA = alpha;
		p.setOriginalRightHandSide();
		if (!p.isCoarset()) {
			// also initializes the solute and rhs matrices
			interpolateSolutionFromCoarserGridSolution(p, residualGoal);
			// corse grid correction
			int vCycleCounter = 0;
			while ((vCycleCounter++ < maxVCylces)
					& coarseGridCorrection(p, residualGoal)) {
			}
			// the solveByRelax execution ensures that solution is achieved
			// at al grid levels
			solveUntilResidualGoal(p, residualGoal, false);
		} else {
			p.createInitialSolution();
			relaxOnCoarsetGrid(p, false);
		}
	}

	/**
	 * Restrict the problem to a coarser grid, solve the problem in that grid,
	 * then prolong the obtained solution to a finner grid
	 * 
	 * @param p
	 * @return
	 */
	private static final void interpolateSolutionFromCoarserGridSolution(
			PoissonProblem p, double residualGoal) {
		solveWithMultigrid(p.pCoarser, residualGoal, ALPHA);
		p.interpolateTo(p.s, p.pCoarser.s);
	}

	/**
	 * Recursive function to solve the coarse grid correction
	 * 
	 * @param p
	 * @param s
	 * @param rhs
	 * @return true if coarse grid correction should be conducted again
	 */
	private static final boolean coarseGridCorrection(PoissonProblem p,
			double residualGoal) {
		if (p.isCoarset()) {
			// if coarser grid is reached, solve equation directly
			relaxOnCoarsetGrid(p, false);
			return false;
		}
		// pre-processing steps
		preProcess(p);
		// get the problem restricted to a corser grid
		PoissonProblem pCoarse = p.pCoarser;
		//
		p.setToLOperatorMatrix(p.temp);
		p.restrictTo(pCoarse.temp, p.temp);
		p.restrictTo(pCoarse.s, p.s);
		// optimization
		p.restrictTo(pCoarse.rhs, p.rhs);
		p.trerr = ALPHA * pCoarse.downwardVOperations();
		// solve
		// corse grid correction
		int vCycleCounter = 0;
		while (vCycleCounter++ < maxVCylces
				&& coarseGridCorrection(pCoarse, residualGoal))
			;
		// correct at finner grid
		p.restrictTo(pCoarse.temp, p.s);
		MultiGridOperations.subtractTo(pCoarse.s, pCoarse.temp, pCoarse.temp);
		p.interpolateTo(p.temp, pCoarse.temp);
		MultiGridOperations.addTo(p.s, p.temp); // update solution
		// pos-processing steps
		postProcess(p);
		// compute residual
		p.setToLOperatorMatrix(p.temp);
		MultiGridOperations.subtractTo(p.temp, p.rhs, p.temp);
		p.res = p.anorm2Temp();
		return (p.res >= p.trerr);
	}

	/**
	 * Solve using Red-Black relaxation steps
	 * 
	 * @param log
	 *            TODO
	 * 
	 */
	public static double[][] solveUntilResidualGoal(PoissonProblem p,
			double residualGoal, boolean log) {
		return performRelaxationSteps(p, residualGoal, maxIterations, log);
	}

	/**
	 * Solve using Red-Black relaxation steps
	 * 
	 */
	public static double[][] relaxOnCoarsetGrid(PoissonProblem p, boolean solve) {
		return performRelaxationSteps(p, 0, NCOARSEST, solve);
	}

	/**
	 * Perform the preprocessing iterations
	 * 
	 * @param p
	 * @param u
	 * @return the updated solution matrix
	 */
	public static double[][] preProcess(PoissonProblem p) {
		return performRelaxationSteps(p, 0, NPREPROCESSING, false);
	}

	/**
	 * Perform the postprocessing iterations
	 * 
	 * @param p
	 * @param u
	 * @return the updated solution matrix
	 */
	public static double[][] postProcess(PoissonProblem p) {
		return performRelaxationSteps(p, 0, NPOSTPROCESSING, false);
	}

	/**
	 * Performs a number of relaxation steps (nSteps). if a residual goal
	 * greater than zero is set, then realaxation iterations are carried out
	 * only until the goal is reached.
	 * 
	 * @param p
	 * @param u
	 * @param residualGoal
	 * @param nSteps
	 *            number of steps
	 * @return the reference to the solution u (the same as the input)
	 */
	private static double[][] performRelaxationSteps(PoissonProblem p,
			double residualGoal, int nSteps, boolean solving) {
		double[][] u = p.s;
		double[][] rhs = p.rhs;
		double sumResidual = 0;
		double lop; // the L-operator
		double dLop; // the derivative of the L-operator
		double r;
		int isw, jsw;
		double count = 0;
		// iterate for a given until maximum number
		int it;
		for (it = 0; it < nSteps; it++) {
			sumResidual = 0;
			count = 0;
			// alternate passes between red and black
			isw = 0;
			for (int pass = 0; pass <= 1; pass++, isw = 1 - isw) {
				jsw = isw;
				// execute relaxation step for each grid node
				for (int i = 0; i < u.length; i++, jsw = 1 - jsw) {
					for (int j = jsw; j < u[i].length; j += 2) {
						// compute only for points inside domain
						if (p.pointIsInsideDomain(i, j)) {
							lop = p.getLOperator(i, j);
							dLop = p.getLOperatorDerivative(i, j);
							r = (lop - rhs[i][j]) / dLop;
							// add the square of the residual to the cumulative
							// sum
							sumResidual += ScalarMath.square(r);
							count += 1.0;
							u[i][j] -= r;
						}
					}
				}
			}
			if (Math.sqrt(sumResidual / count) < residualGoal) {
				break;
			}
		}
		if (solving) {
			// print solution
			System.out.println("n = " + u.length + " solved in " + it
					+ " iterations, r = " + Math.sqrt(sumResidual / count)
					+ " of " + residualGoal);
		}
		return u;
	}
}
