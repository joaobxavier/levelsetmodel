package edu.harvard.sysbio.gavea.reaction;

import edu.harvard.sysbio.gavea.utils.ScalarMath;

/**
 * Implements growth with Monod kinetics
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Feb 4, 2007
 */
public class MonodGrowth extends Reaction {
	private double _yield;

	private double _qS;

	private double _kS;

	private double aux1;

	private double aux2;

	/**
	 * @param qS
	 *            the substrate uptake rate
	 * @param yield
	 *            the growth yield
	 * @param kS
	 *            the saturation constant
	 */
	public MonodGrowth(double qS, double yield, double kS) {
		_qS = qS;
		_yield = yield;
		_kS = kS;
		aux1 = -_qS * _kS;
		aux2 = -_qS;
	}

	@Override
	public double getSubstrateRateDerivative(double s, double b) {
		return (s <= 0 ? aux2 : aux1 / ScalarMath.square(s + _kS) * b);
	}

	@Override
	public double getSubstrateRate(double s, double b) {
		return (s <= 0 ? 0 : aux2 * s / (s + _kS) * b);
	}

	@Override
	public double getSpecificGrowthRate(double s) {
		return (s <= 0 ? 0 : _qS * _yield * s / (s + _kS));
	}

}
