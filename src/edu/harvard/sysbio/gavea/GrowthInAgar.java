package edu.harvard.sysbio.gavea;

import java.util.Random;

import edu.harvard.sysbio.gavea.fastmarching.*;
import edu.harvard.sysbio.gavea.finitedifferencing.*;
import edu.harvard.sysbio.gavea.multigrid.*;
import edu.harvard.sysbio.gavea.output.*;
import edu.harvard.sysbio.gavea.reaction.*;
import edu.harvard.sysbio.gavea.utils.*;

;

/**
 * Model colony growth on an agar plate
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 12, 2007
 */
public class GrowthInAgar {
	/**
	 * residual goal for solving the pressure field
	 */
	private static final double PRESID = 1e-3;

	private static long SEED = 10;

	private DiffusionReactionInAgar _diffusionReactionProblem;

	private double[][] _speedI;

	private double[][] _speedJ;

	private PressureField _pressureField;

	private double _time;

	private double _diffusivity;

	private double _cBulk;

	private double _sizeSystem;

	private double _dz;

	private int _n;

	private double _dx;

	private double _xMax;

	private Random _random;

	private GradientHandler _gradientHandler;

	private FastMarchingVelocityExtension _fastMarchingHandler;

	private Reaction _reaction;

	private static final boolean TIMMING = true;

	private long _timeBefore;

	public GrowthInAgar() {
	}

	public GrowthInAgar(double diffusivity, double yield, double qS, double kS,
			double cBulk, double xMax, double maintenanceFraction,
			double sizeSystem, double dz, int n) {
		setParameters(diffusivity, yield, qS, kS, cBulk, xMax,
				maintenanceFraction, sizeSystem, dz, n);
	}

	/**
	 * @param diffusivity
	 * @param yield
	 * @param qS
	 * @param kS
	 * @param cBulk
	 * @param sizeSystem
	 * @param dz
	 *            vertical distance of nutrient source
	 * @param n
	 */
	public void setParameters(double diffusivity, double yield, double qS,
			double kS, double cBulk, double xMax, double maintenanceFraction,
			double sizeSystem, double dz, int n) {
		_time = 0;
		_diffusivity = diffusivity;
		_cBulk = cBulk;
		_xMax = xMax;
		_sizeSystem = sizeSystem;
		_dz = dz;
		_n = n;
		_dx = _sizeSystem / _n;
		// create matrices
		_speedI = new double[n][n];
		_speedJ = new double[n][n];
		_random = new Random(SEED);
		// create the solvers and get the references of some of the large
		// matrices which are initialized withing the solvers:
		// Reaction r = new MonodGrowth(qS, yield, kS);
		double kd = qS * maintenanceFraction;
		_reaction = new MonodGrowthWithDecay(qS, yield, kS, kd);
		_diffusionReactionProblem = new DiffusionReactionInAgar(_diffusivity,
				_sizeSystem, _dz, _reaction, _cBulk, _xMax, n);
		_diffusionReactionProblem.setResidualGoal(cBulk * 1e-3);
		//
		_pressureField = new PressureField(sizeSystem, n);
		_pressureField.setResidualGoal(PRESID);
		//
		_gradientHandler = new GradientHandler(_n, _dx);
		_gradientHandler.setLevelSet(_pressureField.getDistanceMatrix());
		//
		_fastMarchingHandler = new FastMarchingVelocityExtension(_n, _dx);
		_fastMarchingHandler.setLevelSetAndSpeeds(_pressureField
				.getDistanceMatrix(), _speedI, _speedJ);
	}

	/**
	 * Initialize the system with same size circles placed at random locations
	 * 
	 * 
	 * @param nCircles
	 * @param radius
	 */
	public void initializeSystem(double[][] ls) {
		// copy, reference should not be overwriten
		_pressureField.copyValuesToDistanceMatrix(ls);
	}

	/**
	 * Initialize the system with same size circles placed at random locations
	 * 
	 * 
	 * @param nCircles
	 * @param radius
	 */
	public void initializeSystem(int nCircles, double radius) {
		double[][] circles = new double[nCircles][2];
		// create random locations for circles
		for (int i = 0; i < nCircles; i++) {
			// i coordinate
			circles[i][0] = _random.nextDouble() * ((double) _n);
			// j coordinate
			circles[i][1] = _random.nextDouble() * ((double) _n);
		}
		// initialize the levelset function with the distance to the closest
		// circle border
		double distance;
		for (int i = 0; i < _n; i++)
			for (int j = 0; j < _n; j++) {
				// initialize with the maximum distance possible
				distance = ScalarMath.square(_n);
				// iterate through circles to find the closest
				for (int c = 0; c < nCircles; c++) {
					double distI = circles[c][0] - i;
					distI = Math.abs(distI);
					distI = Math.min(distI, _n - distI);
					double distJ = circles[c][1] - j;
					distJ = Math.abs(distJ);
					distJ = Math.min(distJ, _n - distJ);
					distance = Math
							.min(distance, ScalarMath.norm(distI, distJ));
				}
				// assign distance to border (negative values for inside,
				// positive for outside)
				_pressureField.setValueInDistanceMatrix((distance * _dx)
						- radius, i, j);
			}
	}

	public void doNIterations(int n) {
		for (int i = 0; i < n; i++)
			do1Iteration();
	}

	/**
	 * Perform 1 set of the iterative cycle
	 */
	public void do1Iteration() {
		// update the biomass matrix
		_diffusionReactionProblem.updateValuesFromDistanceMatrix(_pressureField
				.getDistanceMatrix());
		// solve solute gradient
		beginTiming();
		_diffusionReactionProblem.solveWithMultigrid(0.33, 0);
		long drTime = endTiming();
		// Update the growth rate
		_diffusionReactionProblem.updateSpecificGrowthRate(_pressureField
				.getSpecificGrowthRate());
		// solve for pressure field
		_pressureField.updateIsBorder();
		beginTiming();
		_pressureField.solveWithMultigrid(5e-2, 5);
		long pressureTime = endTiming();
		// get the speed from gradient of pressure
		// alternative
		// _gradientHandler.computeSpeedUpwind(_pressureField.getSolution(),
		// _speedI, _speedJ);
		beginTiming();
		_gradientHandler.computeSpeedUpwindGhostNodeFirstOrder(_pressureField
				.getSolution(), _speedI, _speedJ);
		long speedTime = endTiming();
		// Extending using relaxation steps
		// extending and reinitializing simultanewously
		beginTiming();
		_fastMarchingHandler.reinitiallizeDistanceFunctionBothWays();
		long reinitTime = endTiming();
		// determine the time step through CFL criterium
		double velocityNorm = MatrixMath.maximumNorm(_speedI, _speedJ);
		double dt = (_dx / velocityNorm) * 1.0;
		// advance time
		_time += dt;
		System.out.println("dt = " + dt + "; t = " + _time);
		if (Double.isInfinite(dt))
			throw new RuntimeException(
					"Infinite time step (biofilm not growing)");
		// advance the levelset
		beginTiming();
		_gradientHandler.advanceWithRungeKutta(_speedI, _speedJ, dt);
		long advanceTime = endTiming();
		long totalTime = drTime + pressureTime + speedTime + reinitTime
				+ advanceTime;
		System.out.println("DR " + percentage(drTime, totalTime) + ", P "
				+ percentage(pressureTime, totalTime) + ", speed "
				+ percentage(speedTime, totalTime) + ", reinit "
				+ percentage(reinitTime, totalTime) + ", advance "
				+ percentage(advanceTime, totalTime));
		// System.out.println("DR " + drTime + ", P " + pressureTime + ", speed
		// "
		// + speedTime + ", reinit " + reinitTime + ", advance "
		// + advanceTime);
	}

	/**
	 * External interface for solving reaction-diffusion in matlab
	 * 
	 * @param biomass
	 * @return
	 */
	public double[][] solveSoluteField(double[][] biomass) {
		// update the biomass matrix
		_diffusionReactionProblem.setBiomassValues(biomass);
		// solve solute gradient
		beginTiming();
		_diffusionReactionProblem.solveWithMultigrid(0.33, 0);
		return getSolutes();
	}

	/**
	 * External interface for matlab
	 * 
	 * @param pressure
	 * @param levelset
	 */
	public void computeVelocity(double[][] pressure, double[][] levelset) {
		int n = pressure.length;
		_speedI = new double[n][n];
		_speedJ = new double[n][n];
		// create a new instance of the handler
		_gradientHandler = new GradientHandler(n, 1 / ((double) (n)));
		_gradientHandler.setLevelSet(levelset);
		_gradientHandler.computeSpeedUpwindGhostNodeFirstOrder(pressure,
				_speedI, _speedJ);
	}

	/**
	 * @param x
	 * @param levelset
	 */
	public void computeCentralDifference(double[][] x) {
		int n = x.length;
		double dx = 1 / ((double) (n));
		_speedI = new double[n][n];
		_speedJ = new double[n][n];
		// create a new instance of the handler
		GradientHandler.computeSpeedCentralDifferenceGradient(x, _speedI,
				_speedJ, n, dx);
	}

	private String percentage(long t, long total) {
		double percentage = 100.0 * ((double) t) / ((double) total);
		return "" + percentage + " %";
	}

	/**
	 * save timing
	 */
	private void beginTiming() {
		_timeBefore = System.currentTimeMillis();
	}

	/**
	 * return elapsed time since last begin timing
	 */
	private long endTiming() {
		return System.currentTimeMillis() - _timeBefore;
	}

	/**
	 * @return the specific growth rate matrix
	 */
	public double[][] getSpecificGrowthRate() {
		return _pressureField.getSpecificGrowthRate();
	}

	/**
	 * @return the distance matrix used in levelset method
	 */
	public double[][] getDistanceMatrix() {
		// return _pressureField.getDistanceMatrix();
		return _gradientHandler.getLevelSet();
	}

	/**
	 * @return the distance matrix used in levelset method
	 */
	public double[][] getBiomass() {
		return _diffusionReactionProblem.getBiomass();
	}

	/**
	 * @return the solute concentration field
	 */
	public double[][] getSolutes() {
		return _diffusionReactionProblem.getSolution();
	}

	/**
	 * @return the pressure field
	 */
	public double[][] getPressure() {
		return _pressureField.getSolution();
	}

	/**
	 * @return the speed in the i direction
	 */
	public double[][] getSpeedI() {
		return _speedI;
	}

	/**
	 * @return the speed in the j direction
	 */
	public double[][] getSpeedJ() {
		return _speedJ;
	}

	/**
	 * @return the present time
	 */
	public double getTime() {
		return _time;
	}

	public static long getSEED() {
		return SEED;
	}

	public static void setRandomNumberGeneratorSeed(long seed) {
		SEED = seed;
	}

	/**
	 * @param args
	 *            the output directory
	 */
	public static void main(String[] args) {
		if (args.length == 0)
			throw new RuntimeException("Output directory must be supplied");

		// create the output directory
		OutputDirectory outDir = new OutputDirectory(args[0]);

		// create a system
		int n = 129;

		// parameters
		double xMax = 200;
		double cBulk = 3e-3;
		double diffusivity = 8.33e6; // um2/h
		double yield = 0.5;
		double qS = 5;
		double kS = 1e-5;
		// double L = 100;
		double L = 50;
		// double dz = 40;
		double dz = 5;
		double maintenanceFraction = 0.9962;
		// 0.98 is wholes
		// 0.99 is wiggly wholes
		// 0.992 is long wiggly wholes
		// 0.994 is labyrintish
		// 0.995 worms and labyrinth
		// 0.996 spots and worms
		// 0.997 is small spots

		GrowthInAgar s = new GrowthInAgar(diffusivity, yield, qS, kS, cBulk,
				xMax, maintenanceFraction, L, dz, n);

		// initialize state matrix randonmly
		int nCircles = 70;
		double circleRadius = 0.5;
		s.initializeSystem(nCircles, circleRadius);
		int i = 0;
		// while (s.getTime() < 146.4) {
		boolean stopSimulation = false;
		while (true) {
			// while (s.getTime() < 48) {
			// while (i < 20) {
			i++;
			// for (int i = 0; i < 20; i++) {
			System.out.println("iteration " + i);
			try {
				long tBefore = System.currentTimeMillis();
				s.doNIterations(5);
				System.out.println("iteration took "
						+ (System.currentTimeMillis() - tBefore) + " mSecs");
			} catch (RuntimeException e) {
				stopSimulation = true;
			}
			outDir.createNewThreadToWriteMatrix(s.getBiomass(), "biomass" + i
					+ ".txt");
			outDir.createNewThreadToWriteMatrix(s.getPressure(), "pressure" + i
					+ ".txt");
			outDir.createNewThreadToWriteMatrix(s.getSpecificGrowthRate(),
					"specificGrowthRate" + i + ".txt");
			outDir.createNewThreadToWriteMatrix(s.getSolutes(), "solutes" + i
					+ ".txt");
			outDir.createNewThreadToWriteMatrix(s.getDistanceMatrix(),
					"levelset" + i + ".txt");
			outDir.createNewThreadToWriteMatrix(s.getSpeedI(), "speedI" + i
					+ ".txt");
			outDir.createNewThreadToWriteMatrix(s.getSpeedJ(), "speedJ" + i
					+ ".txt");
			if (stopSimulation) {
				System.out.println("max solute concentration = "
						+ MatrixMath.maxValue(s.getSolutes()));
				System.out.println("min for growth = "
						+ ((MonodGrowthWithDecay) s._reaction)
								.getMinimumSoluteConcentration());
				break;
			}
		}
		// i = 1;
		// outDir.createNewThreadToWriteMatrix(s.getBiomass(), "biomass" + i
		// + ".txt");
		// outDir.createNewThreadToWriteMatrix(s.getPressure(), "pressure" + i
		// + ".txt");
		// outDir.createNewThreadToWriteMatrix(s.getSpecificGrowthRate(),
		// "specificGrowthRate" + i + ".txt");
		// outDir.createNewThreadToWriteMatrix(s.getSolutes(), "solutes" + i
		// + ".txt");
		// outDir.createNewThreadToWriteMatrix(s.getDistanceMatrix(), "levelset"
		// + i + ".txt");
		// outDir.createNewThreadToWriteMatrix(s.getSpeedI(), "speedI" + i
		// + ".txt");
		// outDir.createNewThreadToWriteMatrix(s.getSpeedJ(), "speedJ" + i
		// + ".txt");
		//
		//
		System.out.println("done!");
	}
}
