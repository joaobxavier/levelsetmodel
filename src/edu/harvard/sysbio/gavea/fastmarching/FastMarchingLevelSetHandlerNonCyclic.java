package edu.harvard.sysbio.gavea.fastmarching;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.TreeSet;

import edu.harvard.sysbio.gavea.utils.MatrixMath;
import edu.harvard.sysbio.gavea.utils.ScalarMath;

/**
 * A fast marching algorithm that works for non-cyclic borders
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Feb 3, 2007
 */
public class FastMarchingLevelSetHandlerNonCyclic {
	private static final double F = 1.0;

	private double[][] _levelSet;

	private double _dx;

	private int _iSize;

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
	 * Initialize the handler
	 * 
	 * @param n
	 * @param dx
	 */
	public FastMarchingLevelSetHandlerNonCyclic(int n, double dx) {
		_iSize = n;
		_dx = dx;
	}

	/**
	 * set the Level-set matrix
	 * 
	 * @param ls
	 */
	public void setLevelSet(double[][] ls) {
		_levelSet = ls;
	}

	/**
	 * Get the Level-set matrix
	 * 
	 * @return the instance to the level set matrix
	 */
	public double[][] getLevelSet() {
		return _levelSet;
	}

	private class OutsideThread extends Thread {
		FastMarchingLevelSetHandlerNonCyclic h;

		/**
		 * Copies the matrix to a new instance.
		 * 
		 * @param m
		 * @param name
		 */
		public OutsideThread(FastMarchingLevelSetHandlerNonCyclic h) {
			this.h = h;
		}

		public void run() {
			h.reinitiallizeDistanceFunctionOutside();
		}
	}

	private class InsideThread extends Thread {
		FastMarchingLevelSetHandlerNonCyclic h;

		/**
		 * Copies the matrix to a new instance.
		 * 
		 * @param m
		 * @param name
		 */
		public InsideThread(FastMarchingLevelSetHandlerNonCyclic h) {
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
		// Use this for parallel execution of both sides
		reinitiallizeDistanceFunctionInside();
		reinitiallizeDistanceFunctionOutside();
	}

	public void reinitiallizeDistanceFunctionOutside() {
		// intialize narrow band and far-away sets
		int i, j;
		for (i = 0; i < _levelSet.length; i++) {
			for (j = 0; j < _levelSet[i].length; j++) {
				int iPlus = imposeBorder(i + 1, _iSize);
				int iMinus = imposeBorder(i - 1, _iSize);
				int jPlus = imposeBorder(j + 1, _iSize);
				int jMinus = imposeBorder(j - 1, _iSize);
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
			int iPlus = imposeBorder(minVal.i + 1, _iSize);
			int iMinus = imposeBorder(minVal.i - 1, _iSize);
			int jPlus = imposeBorder(minVal.j + 1, _iSize);
			int jMinus = imposeBorder(minVal.j - 1, _iSize);
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
		int iPlus = imposeBorder(location.i + 1, _iSize);
		int iMinus = imposeBorder(location.i - 1, _iSize);
		int jPlus = imposeBorder(location.j + 1, _iSize);
		int jMinus = imposeBorder(location.j - 1, _iSize);
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

	public void reinitiallizeDistanceFunctionInside() {
		// intialize narrow band and far-away sets
		int i, j;
		for (i = 0; i < _levelSet.length; i++) {
			for (j = 0; j < _levelSet[i].length; j++) {
				int iPlus = imposeBorder(i + 1, _iSize);
				int iMinus = imposeBorder(i - 1, _iSize);
				int jPlus = imposeBorder(j + 1, _iSize);
				int jMinus = imposeBorder(j - 1, _iSize);
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
			int iPlus = imposeBorder(maxVal.i + 1, _iSize);
			int iMinus = imposeBorder(maxVal.i - 1, _iSize);
			int jPlus = imposeBorder(maxVal.j + 1, _iSize);
			int jMinus = imposeBorder(maxVal.j - 1, _iSize);
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
		int iPlus = imposeBorder(location.i + 1, _iSize);
		int iMinus = imposeBorder(location.i - 1, _iSize);
		int jPlus = imposeBorder(location.j + 1, _iSize);
		int jMinus = imposeBorder(location.j - 1, _iSize);
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

	/**
	 * Impose cyclic borders for downstream points
	 * 
	 * @param i
	 * @param n
	 * @return (n+i)% n
	 */
	public static final int imposeBorder(int i, int n) {
		if (i > n - 1)
			return n - 2;
		if (i < 0)
			return 1;
		return i;
	}
}
