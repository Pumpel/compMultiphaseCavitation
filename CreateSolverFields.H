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

	volScalarField p_rgh_;
	volVectorField U_;
	surfaceScalarField phi_;

	// The multiphaseCavitationMixture is responsible for solvingthe phase
	// equations and handling all phase specific tasks
	multiphaseCavitationMixture mixture_;

	// The density field used by the solver. This is the mixture density as
	// calculated by the multiphaseCavitationMixture
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
						const fvMesh& mesh);

	// Disallow copy construct
	CreateSolverFields(const Foam::CreateSolverFields& createFields) = delete;

	~CreateSolverFields();

	// Getter methods
	volScalarField& p_rgh();
	volVectorField& U();
	surfaceScalarField& phi();
	volScalarField& rho();
	dimensionedScalar& pMin();
	multiphaseCavitationMixture& mixture();

	uniformDimensionedVectorField& g();
	uniformDimensionedScalarField& hRef();
	dimensionedScalar& ghRef();
	volScalarField& gh();
	surfaceScalarField& ghf();
	volScalarField& K();

};
}



#endif /* CREATESOLVERFIELDS_H_ */
