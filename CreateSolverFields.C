
#include "CreateSolverFields.H"

CreateSolverFields::CreateSolverFields(	const Time& runTime,
										const fvMesh& mesh) :
		runTime_(runTime),
				mesh_(mesh),
				p_rgh_(
						IOobject(
								"p_rgh",
								runTime.timeName(),
								mesh,
								IOobject::MUST_READ,
								IOobject::AUTO_WRITE),
						mesh),
				U_(
						IOobject(
								"U",
								runTime.timeName(),
								mesh,
								IOobject::MUST_READ,
								IOobject::AUTO_WRITE),
						mesh),
				phi_(
						IOobject(
								"phi",
								runTime.timeName(),
								mesh,
								IOobject::READ_IF_PRESENT,
								IOobject::AUTO_WRITE),
						fvc::flux(U_)),
				mixture_(multiphaseCavitationMixture(U_, phi_)),
				rho_(
						IOobject(
								"rho",
								runTime.timeName(),
								mesh,
								IOobject::READ_IF_PRESENT,
								IOobject::AUTO_WRITE),
						mixture_.rho()),
				pMin_(dimensionedScalar("pMin", dimPressure, mixture_)),
				g_(
						IOobject(
								"g",
								runTime.constant(),
								mesh,
								IOobject::MUST_READ,
								IOobject::NO_WRITE)),
				hRef_(
						IOobject(
								"hRef",
								runTime.constant(),
								mesh,
								IOobject::READ_IF_PRESENT,
								IOobject::NO_WRITE),
						dimensionedScalar("hRef", dimLength, 0)),
				ghRef_(
						dimensionedScalar(
								mag(g_.value()) > SMALL ?
										g_
												& (cmptMag(g_.value())
														/ mag(g_.value()))
														* hRef_ :
										dimensionedScalar(
												"ghRef",
												g_.dimensions() * dimLength,
												0))),
				gh_(volScalarField("gh", (g_ & mesh_.C()) - ghRef_)),
				ghf_(surfaceScalarField("ghf", (g_ & mesh_.Cf()) - ghRef_)),
				K_(volScalarField("K", 0.5 * magSqr(U_)))
{
}

Foam::CreateSolverFields::~CreateSolverFields()
{
}

volScalarField& Foam::CreateSolverFields::p_rgh()
{
	return p_rgh_;
}

volVectorField& Foam::CreateSolverFields::U()
{
	return U_;
}

surfaceScalarField& Foam::CreateSolverFields::phi()
{
	return phi_;
}

volScalarField& Foam::CreateSolverFields::rho()
{
	return rho_;
}

dimensionedScalar& Foam::CreateSolverFields::pMin()
{
	return pMin_;
}

multiphaseCavitationMixture& Foam::CreateSolverFields::mixture()
{
	return mixture_;
}

uniformDimensionedVectorField& Foam::CreateSolverFields::g()
{
	return g_;
}

uniformDimensionedScalarField& Foam::CreateSolverFields::hRef()
{
	return hRef_;
}

dimensionedScalar& Foam::CreateSolverFields::ghRef()
{
	return ghRef_;
}

volScalarField& Foam::CreateSolverFields::gh()
{
	return gh_;
}

surfaceScalarField& Foam::CreateSolverFields::ghf()
{
	return ghf_;
}

volScalarField& Foam::CreateSolverFields::K()
{
	return K_;
}

///*
// * CreateSolverFields.C
// *
// *  Created on: May 8, 2018
// *      Author: timo
// */
//
//
//
//
