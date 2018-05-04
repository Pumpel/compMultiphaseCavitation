/*
 * utilies.H
 *
 *  Created on: May 4, 2018
 *      Author: timo
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include "fvCFD.H"
#include "multiphaseCavitationMixture.H"

namespace Foam
{

class Utilities
{
public:
	Utilities();
	Utilities(const Foam::Utilities& aUtilClass);

	~Utilities();

	scalar readAdjustTimeStep(Foam::Time& runTime);
	scalar readMaxCo(Foam::Time& runTime);
	scalar readMaxDeltaT(Foam::Time& runTime);
	scalar readMaxAcousticCo(Foam::Time& runTime);
	scalar calculateCoNumber(
								Time& runTime,
								fvMesh& mesh,
								surfaceScalarField& phi);
	scalar calculateMeanCoNumber(
									Time& runTime,
									fvMesh& mesh,
									surfaceScalarField& phi);
	scalar calculateAcoustincCoNumber(
										Time& runTime,
										multiphaseCavitationMixture& mixture,
										fvMesh& mesh,
										surfaceScalarField& phi);
	void printCoNumbers(
						scalar CoNumber,
						scalar meanCoNumber,
						scalar acousticCoNumber);

	void setInitialExtendedDeltaT(
									bool adjustTimeStep,
									scalar CoNum,
									scalar maxCo,
									scalar maxAcousticCo,
									scalar acousticCoNum,
									Time& runTime,
									scalar maxDeltaT);

	scalar calculateAlphaCoNum(
								fvMesh& mesh,
								Foam::multiphaseCavitationMixture& mixture,
								surfaceScalarField& phi,
								Time& runTime);

	scalar calculateMeanAlphaCoNum(
									fvMesh& mesh,
									Foam::multiphaseCavitationMixture& mixture,
									surfaceScalarField& phi,
									Time& runTime);

	void printAlphaCoNumbers(scalar meanAlphaCoNum, scalar alphaCoNum);

};

}






#endif /* UTILITIES_H_ */
