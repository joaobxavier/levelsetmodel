package edu.harvard.sysbio.gavea.fastmarching;

import java.util.Comparator;
import java.util.TreeSet;

import edu.harvard.sysbio.gavea.utils.MatrixMath;
import edu.harvard.sysbio.gavea.utils.ScalarMath;

/**
 * A genereal purpose fast-marching level set algorithm. Use for initialization
 * of distance functions.
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - apr 24, 2007
 */
public class GeneralFastMarchingLevelSetHandler {
	private static final double F = 1.0;

	private double[][] _levelSet;

	private double _dx;

	private int _iSize;

	private int _jSize;

	private TreeSet<Location> _narrowBandOutside = new TreeSet<Location>();

	private TreeSet<Location> _narrowBandInside = new TreeSet<Location>();

	/**
	 * Compares location based on its level set value
	 * 
	 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Feb 3, 2007
	 */
	public static class OrderByMatrixLocation implements Comparator<Location> {
		public int compare(Location b1, Location b2) {
			if (b1.equals(b2))
				return 0;
			if (b1.i == b2.i)
				return (b1.j > b2.j ? 1 : -1);
			return (b1.i > b2.i ? 1 : -1);
		}
	}

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
	 * Initializes the handler with a biomass matrix that must be non-zero
	 * everywhere outside the biofilm
	 * 
	 * @param m
	 *            biomass matrix
	 * @param dx
	 *            pixelsize
	 */
	public GeneralFastMarchingLevelSetHandler(double[][] m, double dx) {
		_iSize = m.length;
		_jSize = m[0].length;
		_dx = dx;
		_levelSet = new double[_iSize][_jSize];
		initializeBorderPointFromBiomassMatrix(m);
		reinitiallizeDistanceFunctionBothWays();
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
	 * Initializes the points located at each side of the border
	 * 
	 * @param biomass
	 */
	private void initializeBorderPointFromBiomassMatrix(double biomass[][]) {
		for (int i = 0; i < biomass.length; i++) {
			for (int j = 0; j < biomass[i].length; j++) {
				// get the neighbors that are still in _faraway
				int iPlus = zeroFluxUpwards(i + 1, _iSize);
				int iMinus = zeroFluxDownwards(i - 1, _iSize);
				int jPlus = zeroFluxUpwards(j + 1, _jSize);
				int jMinus = zeroFluxDownwards(j - 1, _jSize);
				if (biomass[i][j] > 0) {
					// inside points
					int n = MatrixMath.numberOfZeroValueNeighbors(biomass, i,
							j, iPlus, iMinus, jPlus, jMinus);
					if (n == 1) {
						// one neighbor, direct distance to surface
						_levelSet[i][j] = -_dx * 0.5;
					} else if (n > 1) {
						// more than one neighbor, diagonal distance to surface
						_levelSet[i][j] = -Math.sqrt(_dx * 2) * 0.25;
					} else {
						// other points, any negative value
						_levelSet[i][j] = Double.NEGATIVE_INFINITY;
					}
				} else {
					// outside points
					int n = MatrixMath.numberOfPositiveValueNeighbors(biomass,
							i, j, iPlus, iMinus, jPlus, jMinus);
					if (n == 1) {
						// one neighbor, direct distance to surface
						_levelSet[i][j] = _dx * 0.5;
					} else if (n > 1) {
						// more than one neighbor, diagonal distance to surface
						_levelSet[i][j] = Math.sqrt(_dx * 2) * 0.25;
					} else {
						// other points, any negative value
						_levelSet[i][j] = Double.POSITIVE_INFINITY;
					}
				}
			}
		}
	}

	/**
	 * Impose zeroflux borders for downstream points
	 * 
	 * @param i
	 * @param n
	 * @return (n+i)% n
	 */
	private static final int zeroFluxDownwards(int i, int n) {
		return (i < 0 ? 0 : i);
	}

	/**
	 * Impose zeroflux borders for uptream points
	 * 
	 * @param i
	 * @param n
	 * @return i% n
	 */
	private static final int zeroFluxUpwards(int i, int n) {
		return (i >= n ? n - 1 : i);
	}

	/**
	 * Perform distance function resinitialization using the fast marching level
	 * set method. Points at the border are kept and the distance function is
	 * computed for all other points. This is executed for both the points
	 * inside and outside the border, which means that the fast marching
	 * algorithm is executed twice
	 */
	public void reinitiallizeDistanceFunctionBothWays() {
		// // Use this for parallel execution of both sides
		reinitiallizeDistanceFunctionInside();
		reinitiallizeDistanceFunctionOutside();
	}

	private void reinitiallizeDistanceFunctionOutside() {
		// intialize narrow band and far-away sets
		int i, j;
		for (i = 0; i < _levelSet.length; i++) {
			for (j = 0; j < _levelSet[i].length; j++) {
				int iPlus = zeroFluxUpwards(i + 1, _iSize);
				int iMinus = zeroFluxDownwards(i - 1, _iSize);
				int jPlus = zeroFluxUpwards(j + 1, _jSize);
				int jMinus = zeroFluxDownwards(j - 1, _jSize);
				if (_levelSet[i][j] > 0) // // the check for points in region
					// // the check for points in border
					if (MatrixMath.anyNeighborIsInside(_levelSet, i, j, iPlus,
							iMinus, jPlus, jMinus)) {
						// points in border are placed in narrowband
						_narrowBandOutside.add(new Location(i, j,
								_levelSet[i][j]));
					} else {
						// and location in distance matrix is set to infinite
						_levelSet[i][j] = Double.POSITIVE_INFINITY;
					}
			}
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
			int iPlus = zeroFluxUpwards(minVal.i + 1, _iSize);
			int iMinus = zeroFluxDownwards(minVal.i - 1, _iSize);
			int jPlus = zeroFluxUpwards(minVal.j + 1, _jSize);
			int jMinus = zeroFluxDownwards(minVal.j - 1, _jSize);
			transferNeighborOutside(iPlus, minVal.j);
			transferNeighborOutside(iMinus, minVal.j);
			transferNeighborOutside(minVal.i, jPlus);
			transferNeighborOutside(minVal.i, jMinus);
			// remove the value from the narrowband
			_narrowBandOutside.remove(minVal);
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
			Location location = new Location(i, j);
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
		int iPlus = zeroFluxUpwards(location.i + 1, _iSize);
		int iMinus = zeroFluxDownwards(location.i - 1, _iSize);
		int jPlus = zeroFluxUpwards(location.j + 1, _jSize);
		int jMinus = zeroFluxDownwards(location.j - 1, _jSize);
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

	private void reinitiallizeDistanceFunctionInside() {
		// intialize narrow band and far-away sets
		int i, j;
		for (i = 0; i < _levelSet.length; i++) {
			for (j = 0; j < _levelSet[i].length; j++) {
				int iPlus = zeroFluxUpwards(i + 1, _iSize);
				int iMinus = zeroFluxDownwards(i - 1, _iSize);
				int jPlus = zeroFluxUpwards(j + 1, _jSize);
				int jMinus = zeroFluxDownwards(j - 1, _jSize);
				if (_levelSet[i][j] <= 0) // // the check for points in region
					// // the check for points in border
					if (MatrixMath.anyNeighborIsOutside(_levelSet, i, j, iPlus,
							iMinus, jPlus, jMinus)) {
						// points in border are placed in narrowband
						// _narrowBandInside.add(new Location(i, j));
						_narrowBandInside.add(new Location(i, j));
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
			int iPlus = zeroFluxUpwards(maxVal.i + 1, _iSize);
			int iMinus = zeroFluxDownwards(maxVal.i - 1, _iSize);
			int jPlus = zeroFluxUpwards(maxVal.j + 1, _jSize);
			int jMinus = zeroFluxDownwards(maxVal.j - 1, _jSize);
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
			Location location = new Location(i, j);
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
		int iPlus = zeroFluxUpwards(location.i + 1, _iSize);
		int iMinus = zeroFluxDownwards(location.i - 1, _iSize);
		int jPlus = zeroFluxUpwards(location.j + 1, _jSize);
		int jMinus = zeroFluxDownwards(location.j - 1, _jSize);
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
