#include "CnSwitchingFunction.h"
#include "mathutilities.h"

using std::vector;
using std::string;

namespace ProtoMol
{
	//  ____________________________________________________ CnSwitchingFunction

	Real CnSwitchingFunction::swcoef[][MAXEQNN] = {
		{10., -15., 6., 0, 0, 0, 0},
		{35., -84., 70., -20., 0, 0, 0},
		{126., -420., 540., -315., 70., 0, 0},
		{0, 0, 0, 0, 0, 0, 0},
		{1716., -9009., 20020., -24024., 16380., -6006., 924.0}};

	Real CnSwitchingFunction::dswcoef[][MAXEQNN] = {
		{30., -60., 30., 0, 0, 0, 0},
		{140., -420., 420., -140., 0, 0, 0},
		{630., -2520., 3780., -2520., 630., 0, 0},
		{0, 0, 0, 0, 0, 0, 0},
		{12012., -72072., 180180., -240240., 180180., -72072., 12012.0}};

	Real CnSwitchingFunction::d2swcoef[][MAXEQNN] = {
		{60., -180., 120., 0, 0, 0, 0},
		{420., -1680., 2100., -840., 0, 0, 0},
		{2520., -12600., 22680., -17640., 5040., 0, 0},
		{0, 0, 0, 0, 0, 0, 0},
		{72072., -504504., 1441440., -2162160., 1801800., -792792., 144144.0}};

	//  --------------------------------------------------------------------  //

	CnSwitchingFunction::CnSwitchingFunction(): mySwitchon(-1.0), myCutoff(0.0), myOrder(0.0), mySwitchoff(0.0)
	{
	}

	//  --------------------------------------------------------------------  //

	CnSwitchingFunction::CnSwitchingFunction(Real switchon, Real cutoff, Real
	                                         order, Real switchoff) : mySwitchon(switchon),
	                                                                  mySwitchon2(switchon * switchon),
	                                                                  myOrder(order), mySwitchoff(switchoff)
	{
		if (switchoff)
			myCutoff = switchoff; //use switchoff if defined 
		else
			myCutoff = cutoff;

		myCutoff2 = myCutoff * myCutoff;

		ordInt = (int)myOrder;
		ordIdx = ordInt - 2;

		if (ordInt == 5)
			ordInt = 2; // fix as 5 dosn't exist

		myIRange[0] = pow(switchon - cutoff, -(ordInt + 1));

		for (int i = 1; i < ordInt + 1; i++)
			myIRange[i] = myIRange[i - 1] / (switchon - cutoff);
	}


	//  --------------------------------------------------------------------  //

	void CnSwitchingFunction::getParameters(vector<Parameter>& parameters) const
	{
		parameters.push_back(Parameter("-switchon", Value(mySwitchon, ConstraintValueType::NotNegative()), Text("Cn swf switch on")));
		parameters.push_back(Parameter("-switchoff", Value(myCutoff, ConstraintValueType::NotNegative()), 0.0, Text("Cn swf switchoff")));
		parameters.push_back(Parameter("-cutoff", Value(myCutoff, ConstraintValueType::Positive()), Text("Cn swf cutoff")));
		parameters.push_back(Parameter("-n", Value(myOrder, ConstraintValueType::Positive()), Text("Cn swf smoothness")));
	}

	//  --------------------------------------------------------------------  //

	CnSwitchingFunction CnSwitchingFunction::make(string& errMsg, vector<Value> values)
	{
		Real switchon, cutoff, order, switchoff;

		values[0].get(switchon);
		values[1].get(switchoff);
		values[2].get(cutoff);
		values[3].get(order);

		if (!values[0].valid() || !values[1].valid() || !values[2].valid()
			|| switchon < 0.0 || cutoff <= 0.0 || switchon >= cutoff
			|| order < 2.0 || order > 6.0)
		{
			if (switchoff)
				errMsg += keyword + " switching function: 0 <= switchon (="
					+ values[0].getString() + ") < switchoff (="
					+ values[1].getString() + ") <= cutoff (="
					+ values[2].getString() + ").";

			else
				errMsg += keyword + " switching function: 0 <= switchon (="
					+ values[0].getString() + ") <= cutoff (="
					+ values[2].getString() + ").";

			errMsg += ", order 2, 3, 4, 6 (=" + values[3].getString() + ").";

			return CnSwitchingFunction();
		}

		return CnSwitchingFunction(switchon, cutoff, order, switchoff);
	}

	//  --------------------------------------------------------------------  //

	Matrix3by3 CnSwitchingFunction::hessian(const Vector3D& rij, Real a) const
	{
		Real sqrta = sqrt(a),
			invA = 1 / a;
		Real c[MAXEQNN + 1];
		Real swDiff;

		swDiff = sqrta - myCutoff;

		c[0] = pow(swDiff, ordInt - 1);

		for (int i = 1; i < ordInt + 2; i++)
			c[i] = c[i - 1] * swDiff;

		Real tm3 = 0.,
			tm4 = 0.;

		for (int i = 0; i < ordInt + 1; i++)
		{
			tm3 += dswcoef[ordIdx][i] * c[i + 1] * myIRange[i];
			tm4 += d2swcoef[ordIdx][i] * c[i] * myIRange[i];
		}

		tm3 /= sqrta;

		tm4 *= invA;

		tm4 -= tm3 * invA;

		return (Matrix3by3(tm3, 0, 0, 0, tm3, 0, 0, 0, tm3) +
			Matrix3by3(rij, rij) * tm4);
	}


	//  --------------------------------------------------------------------  //

	void CnSwitchingFunction::operator()(Real& value, Real& deriv,
	                                     Real distSquared) const
	{
		deriv = 0.0;
		value = 0.0;

		if (distSquared > myCutoff2)
			return;

		else if (distSquared < mySwitchon2)
		{
			value = 1.0;
			return;
		}

		Real c[MAXEQNN + 1],
		     swDiff,
		     sqrtd = sqrt(distSquared);

		swDiff = sqrtd - myCutoff;

		c[0] = pow(swDiff, ordInt);

		for (int i = 1; i < ordInt + 2; i++)
			c[i] = c[i - 1] * swDiff;

		for (int i = 0; i < ordInt + 1; i++)
		{
			value += swcoef[ordIdx][i] * c[i + 1] * myIRange[i];
			deriv += dswcoef[ordIdx][i] * c[i] * myIRange[i];
		}

		deriv /= sqrtd;
	}

	//  --------------------------------------------------------------------  //
}
