package edu.harvard.sysbio.gavea.utils;

/**
 * Scalar math operations
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Jan 12, 2007
 */
public abstract class ScalarMath {
	/**
	 * @param x
	 * @return x squared
	 */
	public static final double square(double x) {
		return x * x;
	}

	public static final double norm(double x, double y) {
		return Math.sqrt(square(x) + square(y));
	}

	/**
	 * Impose cyclic borders for downstream points
	 * 
	 * @param i
	 * @param n
	 * @return (n+i)% n
	 */
	public static final int cyclicDownwards(int i, int n) {
		return (n + i) % n;
	}

	/**
	 * Impose cyclic borders for uptream points
	 * 
	 * @param i
	 * @param n
	 * @return i% n
	 */
	public static final int cyclicUpwards(int i, int n) {
		return i % n;
	}

	/**
	 * @param a
	 * @param b
	 * @param c
	 * @return the minimum of the modulus of the 3 values
	 */
	public static final int minMod(int a, int b, int c) {
		a = (a < 0 ? -a : a);
		b = (b < 0 ? -b : b);
		c = (c < 0 ? -c : c);
		return (a < b ? Math.min(a, c) : Math.min(b, c));
	}

}
