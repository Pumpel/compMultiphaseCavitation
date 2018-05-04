/*
 * Utilities.C
 *
 *  Created on: May 4, 2018
 *      Author: timo
 */


#include "Utilities.H"

Utilities::Utilities()
{
}

Utilities::Utilities(const Foam::Utilities& aUtilClass)
{
}

Utilities::~Utilities()
{
}


Foam::scalar Utilities::readAdjustTimeStep(Foam::Time& runTime)
{
	Foam::scalar adjustTimeStep = runTime.controlDict().lookupOrDefault(
			"adjustTimeStep",
			false);
	return adjustTimeStep;
}

Foam::scalar Utilities::readMaxCo(Foam::Time& runTime)
{
	Foam::scalar maxCo = runTime.controlDict().lookupOrDefault("maxCo", 1.0);
	return maxCo;
}

Foam::scalar Utilities::readMaxDeltaT(Foam::Time& runTime)
{
	scalar maxDeltaT = runTime.controlDict().lookupOrDefault < scalar
			> ("maxDeltaT", GREAT);
	return maxDeltaT;
}

scalar Utilities::readMaxAcousticCo(Foam::Time& runTime)
{
	scalar maxAcousticCo = readScalar(
			runTime.controlDict().lookup("maxAcousticCo"));
	return maxAcousticCo;
}

scalar Utilities::calculateCoNumber(
									Time& runTime,
						fvMesh& mesh,
						surfaceScalarField& phi)
{
	scalar CoNum = 0.0;

	if (mesh.nInternalFaces())
	{
		scalarField sumPhi(fvc::surfaceSum(mag(phi))().primitiveField());
		CoNum = 0.5*gMax(sumPhi/mesh.V().field())*runTime.deltaTValue();
	}
	return CoNum;
}

scalar Utilities::calculateMeanCoNumber(
										Time& runTime,
						fvMesh& mesh,
						surfaceScalarField& phi) {
	scalar meanCoNum = 0.0;

	if (mesh.nInternalFaces())
	{
		scalarField sumPhi(fvc::surfaceSum(mag(phi))().primitiveField());
		meanCoNum = 0.5
				* (gSum(sumPhi) / gSum(mesh.V().field()))
				* runTime.deltaTValue();
	}
	return meanCoNum;
}


scalar Utilities::calculateAcoustincCoNumber(
												Time& runTime,
						multiphaseCavitationMixture& mixture,
						fvMesh& mesh,
						surfaceScalarField& phi)
{
	scalar acousticCoNum = 0.0;
	tmp<volScalarField> psi = mixture.mixturePsi();

	if (mesh.nInternalFaces())
	{
		scalarField sumPhi
		(
			fvc::surfaceSum(mag(phi))().primitiveField()
		);

		acousticCoNum = 0.5*gMax
		(
			fvc::surfaceSum
			(
				fvc::interpolate(scalar(1)/sqrt(psi.ref()))*mesh.magSf()
			)().primitiveField()/mesh.V().field()
		)*runTime.deltaTValue();
	}
	return acousticCoNum;
}

void Utilities::printCoNumbers(
					scalar CoNumber,
					scalar meanCoNumber,
					scalar acousticCoNumber)
{
	Info<< "phi Courant Number mean: " << meanCoNumber
			<< " max: "
			<< CoNumber
			<< " acoustic max: "
			<< acousticCoNumber
			<< endl;
}

void Utilities::setInitialExtendedDeltaT(
											bool adjustTimeStep,
								scalar CoNum,
								scalar maxCo,
								scalar maxAcousticCo,
								scalar acousticCoNum,
								Time& runTime,
								scalar maxDeltaT)
{
	if (adjustTimeStep)
	{
		if (CoNum > SMALL)
		{
			scalar maxDeltaTFact = min(
					maxCo / (CoNum + SMALL),
					maxAcousticCo / (acousticCoNum + SMALL));
			runTime.setDeltaT(
					min(maxDeltaTFact * runTime.deltaTValue(), maxDeltaT));
		}
	}
}

scalar Utilities::calculateAlphaCoNum(
										fvMesh& mesh,
							Foam::multiphaseCavitationMixture& mixture,
							surfaceScalarField& phi,
							Time& runTime)
{
	scalar alphaCoNum = 0.0;
	if (mesh.nInternalFaces())
	{
		scalarField sumPhi(
				mixture.nearInterface()().primitiveField()
						* fvc::surfaceSum(mag(phi))().primitiveField());

		alphaCoNum = 0.5
				* gMax(sumPhi / mesh.V().field())
				* runTime.deltaTValue();
	}
	return alphaCoNum;
}


scalar Utilities::calculateMeanAlphaCoNum(
											fvMesh& mesh,
								Foam::multiphaseCavitationMixture& mixture,
								surfaceScalarField& phi,
								Time& runTime)
{
	scalar meanAlphaCoNum = 0.0;
	if (mesh.nInternalFaces())
	{
		scalarField sumPhi(
				mixture.nearInterface()().primitiveField()
						* fvc::surfaceSum(mag(phi))().primitiveField());
		meanAlphaCoNum = 0.5
				* (gSum(sumPhi) / gSum(mesh.V().field()))
				* runTime.deltaTValue();
	}
	return meanAlphaCoNum;
}

void Utilities::printAlphaCoNumbers(scalar meanAlphaCoNum, scalar alphaCoNum)
{
	Info<< "Interface Courant Number mean: " << meanAlphaCoNum
			<< " max: "
			<< alphaCoNum
			<< endl;
}


