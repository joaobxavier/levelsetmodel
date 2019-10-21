package edu.harvard.sysbio.gavea.fastmarching;

import java.util.Iterator;
import java.util.TreeSet;

import edu.harvard.sysbio.gavea.utils.MatrixMath;
import edu.harvard.sysbio.gavea.utils.ScalarMath;

/**
 * Handles execution of the fast-marching levelset method for reinitialization
 * of the distance function and extension of the velocity components outside the
 * biomass, along the normal direction
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Feb 4, 2007
 */
public class FastMarchingVelocityExtension {
	private static final double F = 1.0;

	private double[][] _levelSet;

	private double[][] _speedI;

	private double[][] _speedJ;

	private Location[][] _locations;

	private double _dx;

	private int _iSize;

	private TreeSet<Location> _narrowBandOutside = new TreeSet<Location>();

	private TreeSet<Location> _narrowBandInside = new TreeSet<Location>();

	/**
	 * Stores a location and its level set value
	 * 
	 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Feb 3, 2007
	 */
	private class Location implements Comparable {

		int i;

		int j;

		double value;

		public int compareTo(Object o) {
			Location loc = (Location) o;
			if (this.equals(o))
				return 0;
			if (value == loc.value) {
				if (i == loc.i)
					return (j > loc.j ? 1 : -1);
				return (i > loc.i ? 1 : -1);
			}
			return (value > loc.value ? 1 : -1);
		}

		/**
		 * For comparison purposes only
		 * 
		 * @param i
		 * @param j
		 */
		public Location(int i, int j) {
			this.i = i;
			this.j = j;
		}

		public Location(int i, int j, double levelSetValue) {
			this.i = i;
			this.j = j;
			this.value = levelSetValue;
		}

		/**
		 * Set values
		 * 
		 * @param i
		 * @param j
		 */
		public void setValue(double value) {
			this.value = value;
		}

		@Override
		public boolean equals(Object obj) {
			return (this.i == ((Location) obj).i)
					& (this.j == ((Location) obj).j);
		}

		@Override
		public String toString() {
			return "i = " + i + ", j = " + j + ", val = " + value;
		}
	}

	/**
	 * Initialize the handler
	 * 
	 * @param n
	 * @param dx
	 */
	public FastMarchingVelocityExtension(int n, double dx) {
		_iSize = n;
		_dx = dx;
		_locations = new Location[n][n];
		for (int i = 0; i < _locations.length; i++) {
			for (int j = 0; j < _locations[i].length; j++) {
				_locations[i][j] = new Location(i, j);
			}
		}
	}

	/**
	 * set the Level-set matrix and the speed matrices
	 * 
	 * @param ls
	 * @param speedI
	 * @param speedJ
	 */
	public void setLevelSetAndSpeeds(double[][] ls, double[][] speedI,
			double[][] speedJ) {
		_levelSet = ls;
		_speedI = speedI;
		_speedJ = speedJ;
	}

	/**
	 * Get the Level-set matrix
	 * 
	 * @return the instance to the level set matrix
	 */
	public double[][] getLevelSet() {
		return _levelSet;
	}

	/**
	 * @return the speedI
	 */
	public double[][] getSpeedI() {
		return _speedI;
	}

	/**
	 * @return the speedJ
	 */
	public double[][] getSpeedJ() {
		return _speedJ;
	}

	private class OutsideThread extends Thread {
		FastMarchingVelocityExtension h;

		/**
		 * Copies the matrix to a new instance.
		 * 
		 * @param m
		 * @param name
		 */
		public OutsideThread(FastMarchingVelocityExtension h) {
			this.h = h;
		}

		public void run() {
			h.reinitiallizeDistanceAndExtendVelocitiesOutside();
		}
	}

	private class InsideThread extends Thread {
		FastMarchingVelocityExtension h;

		/**
		 * Copies the matrix to a new instance.
		 * 
		 * @param m
		 * @param name
		 */
		public InsideThread(FastMarchingVelocityExtension h) {
			this.h = h;
		}

		public void run() {
			h.reinitiallizeDistanceFunctionInside();
		}
	}

	/**
	 * Perform distance function resinitialization using the fast marching level
	 * set method. Points at the border are kept and the distance function is
	 * computed for all other points. This is executed for both the points
	 * inside and outside the border, which means that the fast marching
	 * algorithm is executed twice
	 * 
	 * TODO paralelize execution using multithreading
	 */
	public void reinitiallizeDistanceFunctionBothWays() {
		// // Use this for parallel execution of both sides
		// Thread tout = new OutsideThread(this);
		// tout.start();
		// Thread tin = new InsideThread(this);
		// tin.start();
		// try {
		// tout.join();
		// tin.join();
		// } catch (InterruptedException e) {
		// throw new RuntimeException(e);
		// }
		reinitiallizeDistanceFunctionInside();
		reinitiallizeDistanceAndExtendVelocitiesOutside();
	}

	public void reinitiallizeDistanceAndExtendVelocitiesOutside() {
		// intialize narrow band and far-away sets
		for (int i = 0; i < _levelSet.length; i++) {
			for (int j = 0; j < _levelSet[i].length; j++) {
				int iPlus = ScalarMath.cyclicUpwards(i + 1, _iSize);
				int iMinus = ScalarMath.cyclicDownwards(i - 1, _iSize);
				int jPlus = ScalarMath.cyclicUpwards(j + 1, _iSize);
				int jMinus = ScalarMath.cyclicDownwards(j - 1, _iSize);
				if (_levelSet[i][j] > 0) // // the check for points in region
					// // the check for points in border
					if (MatrixMath.anyNeighborIsInside(_levelSet, i, j, iPlus,
							iMinus, jPlus, jMinus)) {
						// points in border are placed in narrowband
						_locations[i][j].setValue(_levelSet[i][j]);
						_narrowBandOutside.add(_locations[i][j]);
					} else {
						// and location in distance matrix is set to infinite
						_levelSet[i][j] = Double.POSITIVE_INFINITY;
					}
			}
		}
		for (Iterator iter = _narrowBandOutside.iterator(); iter.hasNext();) {
			Location loc = (Location) iter.next();
			//
			int i = loc.i;
			int j = loc.j;
			int iPlus = ScalarMath.cyclicUpwards(i + 1, _iSize);
			int iMinus = ScalarMath.cyclicDownwards(i - 1, _iSize);
			int jPlus = ScalarMath.cyclicUpwards(j + 1, _iSize);
			int jMinus = ScalarMath.cyclicDownwards(j - 1, _iSize);
			// update the speed values
			_speedI[i][j] = updateFromUpwindGradient(_speedI, i, j, iMinus,
					iPlus, jMinus, jPlus);
			_speedJ[i][j] = updateFromUpwindGradient(_speedJ, i, j, iMinus,
					iPlus, jMinus, jPlus);
		}
		// iterate until narrow band is empty
		// OrderByLevelSetValue comparator = new OrderByLevelSetValue();
		// int lixo = 0;
		while (_narrowBandOutside.size() > 0) {
			// if (lixo++ > 5000) {
			// System.out.println("done " + lixo);
			// break;
			// }
			// get the minimum value
			Location minVal = _narrowBandOutside.first();
			// get the neighbors that are still in _faraway
			int iPlus = ScalarMath.cyclicUpwards(minVal.i + 1, _iSize);
			int iMinus = ScalarMath.cyclicDownwards(minVal.i - 1, _iSize);
			int jPlus = ScalarMath.cyclicUpwards(minVal.j + 1, _iSize);
			int jMinus = ScalarMath.cyclicDownwards(minVal.j - 1, _iSize);
			transferNeighborOutside(iPlus, minVal.j);
			transferNeighborOutside(iMinus, minVal.j);
			transferNeighborOutside(minVal.i, jPlus);
			transferNeighborOutside(minVal.i, jMinus);
			// remove the value from the narrowband
			_narrowBandOutside.remove(minVal);
			// update the speed values
			_speedI[minVal.i][minVal.j] = updateFromUpwindGradient(_speedI,
					minVal.i, minVal.j, iMinus, iPlus, jMinus, jPlus);
			_speedJ[minVal.i][minVal.j] = updateFromUpwindGradient(_speedJ,
					minVal.i, minVal.j, iMinus, iPlus, jMinus, jPlus);
		}
		// done
	}

	/**
	 * @param narrowBand
	 * @param faraway
	 * @param i
	 * @param j
	 */
	private void transferNeighborOutside(int i, int j) {
		if (Double.isInfinite(_levelSet[i][j])) {
			Location location = _locations[i][j];
			// solve the quadratic for the location
			location.value = solveQuadraticOutside(location);
			// add to the narrow band
			_narrowBandOutside.add(location);
			// update value in distance matrix
			_levelSet[i][j] = location.value;
		}
	}

	private double solveQuadraticOutside(Location location) {
		// get the neighbors that are still in _faraway
		int iPlus = ScalarMath.cyclicUpwards(location.i + 1, _iSize);
		int iMinus = ScalarMath.cyclicDownwards(location.i - 1, _iSize);
		int jPlus = ScalarMath.cyclicUpwards(location.j + 1, _iSize);
		int jMinus = ScalarMath.cyclicDownwards(location.j - 1, _iSize);
		//
		double ux1 = _levelSet[iPlus][location.j];
		double ux2 = _levelSet[iMinus][location.j];
		double uy1 = _levelSet[location.i][jPlus];
		double uy2 = _levelSet[location.i][jMinus];
		//
		double ux = Math.min(ux1, ux2);
		double uy = Math.min(uy1, uy2);
		// TODO pre-compute
		double f = _dx / F;
		//
		// solve the quadratic equation
		double a = 2;
		double b = -2 * (ux + uy);
		double c = ScalarMath.square(ux) + ScalarMath.square(uy)
				- ScalarMath.square(f);
		double delta = ScalarMath.square(b) - 4 * a * c;
		// extraordinary solutions
		double u;
		if (Double.isInfinite(ux) | Double.isInfinite(uy) | delta < 0) {
			if (Double.isInfinite(ux))
				u = uy;
			else if (Double.isInfinite(uy))
				u = ux;
			else
				// the case where (delta < 0), which produces a
				// imaginary root
				u = Math.min(ux, uy);
			return u + f;
		}
		double sol1 = (-b + Math.sqrt(delta)) / (2 * a);
		double sol2 = (-b - Math.sqrt(delta)) / (2 * a);
		return Math.max(sol1, sol2);
	}

	private double updateFromUpwindGradient(double[][] u, int i, int j,
			int iMinus, int iPlus, int jMinus, int jPlus) {
		double dLsdIPlus = -Math.min(_levelSet[iPlus][j] - _levelSet[i][j], 0);
		double dLsdIMinus = Math.max(_levelSet[i][j] - _levelSet[iMinus][j], 0);
		double dLsdJPlus = -Math.min(_levelSet[i][jPlus] - _levelSet[i][j], 0);
		double dLsdJMinus = Math.max(_levelSet[i][j] - _levelSet[i][jMinus], 0);
		// upwind scheme
		return (dLsdIMinus * u[iMinus][j] + dLsdIPlus * u[iPlus][j]
				+ dLsdJMinus * u[i][jMinus] + dLsdJPlus * u[i][jPlus])
				/ (dLsdIMinus + dLsdIPlus + dLsdJMinus + dLsdJPlus);
	}

	public void reinitiallizeDistanceFunctionInside() {
		// intialize narrow band and far-away sets
		int i, j;
		for (i = 0; i < _levelSet.length; i++) {
			for (j = 0; j < _levelSet[i].length; j++) {
				int iPlus = ScalarMath.cyclicUpwards(i + 1, _iSize);
				int iMinus = ScalarMath.cyclicDownwards(i - 1, _iSize);
				int jPlus = ScalarMath.cyclicUpwards(j + 1, _iSize);
				int jMinus = ScalarMath.cyclicDownwards(j - 1, _iSize);
				if (_levelSet[i][j] <= 0) // // the check for points in region
					// // the check for points in border
					if (MatrixMath.anyNeighborIsOutside(_levelSet, i, j, iPlus,
							iMinus, jPlus, jMinus)) {
						// points in border are placed in narrowband
						_locations[i][j].setValue(_levelSet[i][j]);
						_narrowBandInside.add(_locations[i][j]);
					} else {
						// and location in distance matrix is set to infinite
						_levelSet[i][j] = Double.NEGATIVE_INFINITY;
					}
			}
		}
		// iterate until narrow band is empty
		// OrderByLevelSetValue comparator = new OrderByLevelSetValue();
		// while (_narrowBandInside.size() > 0) {
		while (_narrowBandInside.size() > 0) {
			// get the maximum value
			// Location maxVal = _narrowBandInside.last();
			Location maxVal = _narrowBandInside.last();
			// System.out.println("maxVal " + maxVal);
			// System.out.println("_narrowBandInside.contains(maxVal) "
			// + _narrowBandInside.contains(maxVal));
			// System.out.println("maxVal.compareTo(maxVal) "
			// + maxVal.compareTo(maxVal));

			// get the neighbors that are still in _faraway
			int iPlus = ScalarMath.cyclicUpwards(maxVal.i + 1, _iSize);
			int iMinus = ScalarMath.cyclicDownwards(maxVal.i - 1, _iSize);
			int jPlus = ScalarMath.cyclicUpwards(maxVal.j + 1, _iSize);
			int jMinus = ScalarMath.cyclicDownwards(maxVal.j - 1, _iSize);
			transferNeighborInside(iPlus, maxVal.j);
			transferNeighborInside(iMinus, maxVal.j);
			transferNeighborInside(maxVal.i, jPlus);
			transferNeighborInside(maxVal.i, jMinus);
			// remove the value from the narrowband
			// _narrowBandInside.remove(maxVal);
			_narrowBandInside.remove(maxVal);
		}
		// done
	}

	/**
	 * @param narrowBand
	 * @param faraway
	 * @param i
	 * @param j
	 */
	private void transferNeighborInside(int i, int j) {
		if (Double.isInfinite(-_levelSet[i][j])) {
			Location location = _locations[i][j];
			// solve the quadratic for the location
			location.value = solveQuadraticInside(location);
			// add to the narrow band
			_narrowBandInside.add(location);
			// update value in distance matrix
			_levelSet[i][j] = location.value;
		}
	}

	private double solveQuadraticInside(Location location) {
		// get the neighbors that are still in _faraway
		int iPlus = ScalarMath.cyclicUpwards(location.i + 1, _iSize);
		int iMinus = ScalarMath.cyclicDownwards(location.i - 1, _iSize);
		int jPlus = ScalarMath.cyclicUpwards(location.j + 1, _iSize);
		int jMinus = ScalarMath.cyclicDownwards(location.j - 1, _iSize);
		//
		double ux1 = -_levelSet[iPlus][location.j];
		double ux2 = -_levelSet[iMinus][location.j];
		double uy1 = -_levelSet[location.i][jPlus];
		double uy2 = -_levelSet[location.i][jMinus];
		//
		double ux = Math.min(ux1, ux2);
		double uy = Math.min(uy1, uy2);
		// TODO pre-compute
		double f = _dx / F;
		//
		// solve the quadratic equation
		double a = 2;
		double b = -2 * (ux + uy);
		double c = ScalarMath.square(ux) + ScalarMath.square(uy)
				- ScalarMath.square(f);
		double delta = ScalarMath.square(b) - 4 * a * c;
		// extraordinary solutions
		double u;
		if (Double.isInfinite(ux) | Double.isInfinite(uy) | delta < 0) {
			if (Double.isInfinite(ux))
				u = uy;
			else if (Double.isInfinite(uy))
				u = ux;
			else
				// the case where (delta < 0), which produces a
				// imaginary root
				u = Math.min(ux, uy);
			return -(u + f);
		}
		double sol1 = (-b + Math.sqrt(delta)) / (2 * a);
		double sol2 = (-b - Math.sqrt(delta)) / (2 * a);
		return -(Math.max(sol1, sol2));
	}

}
