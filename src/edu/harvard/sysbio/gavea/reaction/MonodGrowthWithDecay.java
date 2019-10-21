package edu.harvard.sysbio.gavea.reaction;

import edu.harvard.sysbio.gavea.utils.ScalarMath;

/**
 * Implements growth with Monod kinetics and first order decay
 * 
 * @author Joao Xavier (jxavier@cgr.harvard.edu) - Feb 4, 2007
 */
public class MonodGrowthWithDecay extends Reaction {
	private double _yield;

	private double _qS;

	private double _kS;

	private double _kDecay;

	private double _uMax;

	private double _minimumS;

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
	public MonodGrowthWithDecay(double qS, double yield, double kS,
			double kDecay) {
		_qS = qS;
		_yield = yield;
		_kS = kS;
		_kDecay = kDecay;
		aux1 = -_qS * _kS;
		aux2 = -_qS;
		_uMax = (_qS - _kDecay) * _yield;
		_minimumS = _kS * _kDecay * yield / _uMax;
	}

	@Override
	public double getSubstrateRateDerivative(double s, double b) {
		if (b == 0)
			return 0;
		return (s <= 0 ? aux2 : aux1 / ScalarMath.square(s + _kS) * b);
	}

	@Override
	public double getSubstrateRate(double s, double b) {
		if (b == 0)
			return 0;
		return (s <= 0 ? 0 : aux2 * s / (s + _kS) * b);
	}

	@Override
	public double getSpecificGrowthRate(double s) {
		// return (s <= _minimumS ? 0 : _uMax * s / (s + _kS) - _kDecay * _yield
		// * _kS / (s + _kS));
		return (_uMax * s / (s + _kS) - _kDecay * _yield * _kS / (s + _kS));
	}

	/**
	 * @return the value of the minimum solute concentration required for growth
	 */
	public double getMinimumSoluteConcentration() {
		return _minimumS;
	}

}
