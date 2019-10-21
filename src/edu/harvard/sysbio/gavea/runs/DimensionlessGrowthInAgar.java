package edu.harvard.sysbio.gavea.runs;

import edu.harvard.sysbio.gavea.GrowthInAgar;
import edu.harvard.sysbio.gavea.output.*;
import edu.harvard.sysbio.gavea.utils.*;

;

/**
 * Model colony growth on an agar plate
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 12, 2007
 */
public class DimensionlessGrowthInAgar extends GrowthInAgar {
	// Dimensionless numbers used in simulation
	private double _thiele; // the thiele modulus

	private double _ksDimless; // the dimensionless kS

	private double _maintenanceFraction; // the fraction of uptake rate that

	// goes to maintenance

	private double _dimLessL;

	private OutputDirectory _outDir;

	/**
	 * @param thiele
	 * @param ksDimless
	 * @param fraction
	 * @param dimLessL
	 * @param n
	 * @param outDir
	 */
	public DimensionlessGrowthInAgar(double thiele, double ksDimless,
			double fraction, double dimLessL, int n, String outDir) {
		super();
		_thiele = thiele;
		_ksDimless = ksDimless;
		_maintenanceFraction = fraction;
		_dimLessL = dimLessL;
		setParameters(1.0, 1.0, _thiele, _ksDimless, 1.0, 1.0,
				_maintenanceFraction, _dimLessL, 1.0, n);
		_outDir = new OutputDirectory(outDir);
	}

	public void run() {
		System.out.println("started");
		// initialize state matrix randonmly
		int nCircles = 70;
		double circleRadius = 0.1; // in dimensionless units
		initializeSystem(nCircles, circleRadius);
		int i = 0;
		// while (s.getTime() < 146.4) {
		boolean stopSimulation = false;
		// while (s.getTime() < 48) {
		while (i < 20) {
			i++;
			// for (int i = 0; i < 20; i++) {
			System.out.println("iteration " + i);
			try {
				long tBefore = System.currentTimeMillis();
				doNIterations(5);
				System.out.println("iteration took "
						+ (System.currentTimeMillis() - tBefore) + " mSecs");
			} catch (RuntimeException e) {
				stopSimulation = true;
			}
			_outDir.createNewThreadToWriteMatrix(getBiomass(), "biomass" + i
					+ ".txt");
			_outDir.createNewThreadToWriteMatrix(getPressure(), "pressure" + i
					+ ".txt");
			_outDir.createNewThreadToWriteMatrix(getSpecificGrowthRate(),
					"specificGrowthRate" + i + ".txt");
			_outDir.createNewThreadToWriteMatrix(getSolutes(), "solutes" + i
					+ ".txt");
			_outDir.createNewThreadToWriteMatrix(getDistanceMatrix(),
					"levelset" + i + ".txt");
			_outDir.createNewThreadToWriteMatrix(getSpeedI(), "speedI" + i
					+ ".txt");
			_outDir.createNewThreadToWriteMatrix(getSpeedJ(), "speedJ" + i
					+ ".txt");
			// allow to finish writing beofre shutting down
			if (stopSimulation) {
				break;
			}
		}
		//
		System.out.println("done!");
	}

	/**
	 * @param args
	 *            the output directory
	 */
	public static void main(String[] args) {
		if (args.length == 0)
			throw new RuntimeException("Output directory must be supplied");
		// create a system
		int n = 129;

		// parameters
		double xMax = 200;
		//double cBulk = 8e-3; // no patterns (completely occupied surface)
		double cBulk = 1e-4;
		//double cBulk = 1e-5; // dots
		double diffusivity = 8.33e-6;
		double qS = 0.47*0.505;
		double kS = 3.5e-4;
		double L = 100*1e-6;
		double dz = 2.8810e-005;
		double maintenanceFraction = 0.0307;

		double thiele = qS * xMax / (diffusivity / ScalarMath.square(dz))
				/ cBulk;
		DimensionlessGrowthInAgar s = new DimensionlessGrowthInAgar(thiele, kS
				/ cBulk, maintenanceFraction, L / dz, n, args[0]);
		s.run();
	}
}
