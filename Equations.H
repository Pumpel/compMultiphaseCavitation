/*
 * PEqn.H
 *
 *  Created on: May 8, 2018
 *      Author: timo
 */

#ifndef EQUATIONS_H_
#define EQUATIONS_H_

#include "fvCFD.H"

namespace Foam
{
class Equations
{
public:

	Equations();
	Equations(const Foam::Equations& eqns);

	~Equations();

	void pEqn();

};
}





#endif /* EQUATIONS_H_ */
