/*
 * CreateFields.H
 *
 *  Created on: May 8, 2018
 *      Author: timo
 */

#ifndef CREATESOLVERFIELDS_H_
#define CREATESOLVERFIELDS_H_

#include "fvCFD.H"
#include "multiphaseCavitationMixture.H"

namespace Foam
{
class CreateSolverFields
{
private:
	const Time& runTime_;
	const fvMesh& mesh_;
	multiphaseCavitationMixture& mixture_;

	volScalarField p_rgh_;
	volVectorField U_;
	surfaceScalarField phi_;
	volScalarField rho_;
	dimensionedScalar pMin_;

	uniformDimensionedVectorField g_;
	uniformDimensionedScalarField hRef_;
	dimensionedScalar ghRef_;
	volScalarField gh_;
	surfaceScalarField ghf_;
	// Field kinetic energy for the turbulence model
	volScalarField K_;

public:
	CreateSolverFields(	const Time& runTime,
						const fvMesh& mesh,
						multiphaseCavitationMixture& mixture);

	// Disallow copy construct
	CreateSolverFields(const Foam::CreateSolverFields& createFields) = delete;

	~CreateSolverFields();

	// Getter methods
	volScalarField& p_rgh();
	volVectorField& U();
	surfaceScalarField& phi();
	volScalarField& rho();
	dimensionedScalar& pMin();

	uniformDimensionedVectorField& g();
	uniformDimensionedScalarField& hRef();
	dimensionedScalar& ghRef();
	volScalarField& gh();
	surfaceScalarField& ghf();
	volScalarField& K();

};
}



#endif /* CREATESOLVERFIELDS_H_ */
