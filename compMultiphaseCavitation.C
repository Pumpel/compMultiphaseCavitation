/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
 \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
 \\/     M anipulation  |
-------------------------------------------------------------------------------
License
 This file is part of OpenFOAM.

 OpenFOAM is free software: you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "pimpleControl.H"
#include "multiphaseCavitationMixture.H"
#include "Utilities.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{

	// Defining the utilities class object
	Foam::Utilities util = Foam::Utilities();


	// TODO Use the function calls instead of the #include files for postProcess.H
	//#include "postProcess.H"

	Foam::argList args(argc, argv);
	if (!args.checkRootCase())
	{
		Foam::FatalError.exit();
	}

	// Create Time
	Foam::Info<< "Create time\n" << Foam::endl;
	Foam::Time runTime(Foam::Time::controlDictName, args);

	// Create Mesh
	Foam::Info
			<< "Create mesh for time = "
			<< runTime.timeName()
			<< Foam::nl
			<< Foam::endl;

	Foam::fvMesh mesh
	(
			Foam::IOobject(
					Foam::fvMesh::defaultRegion,
					runTime.timeName(),
					runTime,
					Foam::IOobject::MUST_READ)
	);

	// Create Control
	#if defined(NO_CONTROL)
	#elif defined(PISO_CONTROL)
		pisoControl piso(mesh);
	#elif defined(PIMPLE_CONTROL)
		pimpleControl pimple(mesh);
	#elif defined(SIMPLE_CONTROL)
		simpleControl simple(mesh);
	#endif

	// Read the control parameters used by setDeltaT before starting the first time step
	bool adjustTimeStep 	= util.readAdjustTimeStep(runTime);
	scalar maxCo 			= util.readMaxCo(runTime);
	scalar maxDeltaT 		= util.readMaxDeltaT(runTime);
	scalar maxAcousticCo 	= util.readMaxAcousticCo(runTime);


	// Creating Fields
	Info<< "Reading field p_rgh\n" << endl;
	volScalarField p_rgh
	(
			IOobject(
					"p_rgh",
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE),
			mesh
	);

	Info<< "Reading field U\n" << endl;
	volVectorField U
	(
			IOobject(
					"U",
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE),
			mesh
	);

	// Creates and initialises the relative face-flux field phi.
	Info<< "Reading/calculating face flux field phi\n" << endl;
	surfaceScalarField phi
	(
			IOobject(
					"phi",
					runTime.timeName(),
					mesh,
					IOobject::READ_IF_PRESENT,
					IOobject::AUTO_WRITE),
			fvc::flux(U)
	);

	Info<< "Constructing multiphaseCavitationMixture\n" << endl;
	multiphaseCavitationMixture mixture(U, phi);

	volScalarField rho
	(
			IOobject(
					"rho",
					runTime.timeName(),
					mesh,
					IOobject::READ_IF_PRESENT,
					IOobject::AUTO_WRITE),
			mixture.rho()
	);

	dimensionedScalar pMin("pMin", dimPressure, mixture);

	mesh.setFluxRequired(p_rgh.name());

	// Reading the gravitational acceleration
	Info << "\nReading g" << endl;
	uniformDimensionedVectorField g(
			IOobject(
					"g",
					runTime.constant(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE));

	Info << "\nReading hRef" << endl;
	uniformDimensionedScalarField hRef(
			IOobject(
					"hRef",
					runTime.constant(),
					mesh,
					IOobject::READ_IF_PRESENT,
					IOobject::NO_WRITE),
			dimensionedScalar("hRef", dimLength, 0));

	Info << "Calculating field g.h\n" << endl;
	dimensionedScalar ghRef(
			mag(g.value()) > SMALL ?
										g & (cmptMag(g.value()) / mag(g.value())) * hRef :
										dimensionedScalar(
												"ghRef",
												g.dimensions() * dimLength,
												0));
	volScalarField gh("gh", (g & mesh.C()) - ghRef);
	surfaceScalarField ghf("ghf", (g & mesh.Cf()) - ghRef);

	// Construct compressible turbulence model
	autoPtr<compressible::turbulenceModel> turbulence
	(
			compressible::turbulenceModel::New(
					rho,
					U,
					mixture.rhoPhi(),
					mixture)
	);

	Info<< "Creating field kinetic energy K\n" << endl;
	volScalarField K("K", 0.5*magSqr(U));


	// Calculate the Courant-Numbers
	scalar CoNum 			= util.calculateCoNumber(runTime, mesh, phi);
	scalar meanCoNum 		= util.calculateMeanCoNumber(runTime, mesh, phi);
	scalar acousticCoNum 	= util.calculateAcoustincCoNumber(runTime,mixture,mesh,phi);
	util.printCoNumbers(CoNum, meanCoNum, acousticCoNum);

	// Setting the first time step deltaT. Here, the extendend version for
	// compressible flow is used.
	util.setInitialExtendedDeltaT(adjustTimeStep,
								CoNum,
								maxCo,
								maxAcousticCo,
								acousticCoNum,
								runTime,
								maxDeltaT);

	volScalarField& p = mixture.p();
	volScalarField& T = mixture.T();

	turbulence->validate();

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info << "\nStarting time loop\n" << endl;

	while (runTime.run())
	{
		// Read the control values for this time step
		bool adjustTimeStep = util.readAdjustTimeStep(
				runTime);
		scalar maxCo = util.readMaxCo(runTime);
		scalar maxDeltaT = util.readMaxDeltaT(runTime);
		scalar maxAcousticCo = util.readMaxAcousticCo(
				runTime);
		scalar maxAlphaCo(
				readScalar(runTime.controlDict().lookup("maxAlphaCo")));

		// Calculate the Courant-Numbers
		scalar CoNum = util.calculateCoNumber(
				runTime,
				mesh,
				phi);
		scalar meanCoNum = util.calculateMeanCoNumber(
				runTime,
				mesh,
				phi);
		scalar acousticCoNum = util
				.calculateAcoustincCoNumber(
				runTime,
				mixture,
				mesh,
				phi);
		util.printCoNumbers(
				CoNum,
				meanCoNum,
				acousticCoNum);

		scalar alphaCoNum = util.calculateAlphaCoNum(
				mesh,
				mixture,
				phi,
				runTime);
		scalar meanAlphaCoNum =
				util.calculateMeanAlphaCoNum(
						mesh,
						mixture,
						phi,
						runTime);
		util.printAlphaCoNumbers(
				meanAlphaCoNum,
				alphaCoNum);

		// Adjusting the time step based on the maximum Courant Numbers
		if (adjustTimeStep)
		{
			scalar maxDeltaTFact = min(
					maxCo / (CoNum + SMALL),
					maxAcousticCo / (acousticCoNum + SMALL));

			scalar deltaTFact = min(
					min(maxDeltaTFact, 1.0 + 0.1 * maxDeltaTFact),
					1.2);

			runTime.setDeltaT(
					min(deltaTFact * runTime.deltaTValue(), maxDeltaT));

			Info << "deltaT = " << runTime.deltaTValue() << endl;
		}

		runTime++;

		Info << "Time = " << runTime.timeName() << nl << endl;

		// --- Pressure-velocity PIMPLE corrector loop
		while (pimple.loop())
		{
			mixture.solve();

			solve(fvm::ddt(rho) + fvc::div(mixture.rhoPhi()));

			// The momentum equation
			fvVectorMatrix UEqn(
					fvm::ddt(rho, U) + fvm::div(mixture.rhoPhi(), U) - fvm::Sp(
							fvc::ddt(rho) + fvc::div(mixture.rhoPhi()),
							U) + turbulence->divDevRhoReff(U));

			UEqn.relax();

			if (pimple.momentumPredictor())
			{
				solve(
						UEqn == fvc::reconstruct(
								(mixture.surfaceTensionForce() - ghf * fvc::snGrad(
										rho) - fvc::snGrad(p_rgh)) * mesh.magSf()));

				K = 0.5 * magSqr(U);
			}

			// The temperature equation
			{
				fvScalarMatrix TEqn(
						fvm::ddt(rho, T) + fvm::div(mixture.rhoPhi(), T)
								- fvm::laplacian(
										mixture.alphaEff(turbulence->mut()),
										T)
								+ (fvc::div(fvc::absolute(phi, U), p)
										+ fvc::ddt(rho, K)
										+ fvc::div(mixture.rhoPhi(), K))
										* mixture
								.rCv());

				TEqn.relax();
				TEqn.solve();

				mixture.correct();
			}

			// --- Pressure corrector loop
			while (pimple.correct())
			{
				// The pressure equation
				{
					volScalarField rAU("rAU", 1.0 / UEqn.A());
					surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
					volVectorField HbyA(
							constrainHbyA(rAU * UEqn.H(), U, p_rgh));
					surfaceScalarField phiHbyA(
							"phiHbyA",
							fvc::flux(HbyA)
									+ fvc::interpolate(rho * rAU)
											* fvc::ddtCorr(
									U,
									phi));

					surfaceScalarField phig(
							(mixture.surfaceTensionForce() - ghf * fvc::snGrad(
									rho)) * rAUf * mesh.magSf());

					phiHbyA += phig;

					// Update the pressure BCs to ensure flux consistency
					constrainPressure(p_rgh, U, phiHbyA, rAUf);

					PtrList<fvScalarMatrix> p_rghEqnComps(
							mixture.phases().size());

					label phasei = 0;
					forAllConstIter
					(
							PtrDictionary<phaseModel>,
							mixture.phases(),
							phase
					){
					const rhoThermo& thermo = phase().thermo();
					const volScalarField& rho = thermo.rho()();

					p_rghEqnComps.set
					(
							phasei,
							(
									fvc::ddt(rho) +
									thermo.psi()*correction(fvm::ddt(p_rgh))
									+ fvc::div(phi, rho)
									- fvc::Sp(fvc::div(phi), rho)
							).ptr()
					);

					phasei++;
				}

				//Coefficients of the mass transfer
					const dimensionedScalar pSat = mixture.cavitationModel()
							->pSat();
					volScalarField pSatField(
							IOobject(
									"pSatField",
									mesh.time().timeName(),
									mesh,
									IOobject::NO_READ,
									IOobject::NO_WRITE),
							mesh,
							dimensionedScalar(
									"pSatField",
									pSat.dimensions(),
									0.0));
					pSatField.dimensions().reset(pSat.dimensions());
					pSatField = mixture.cavitationModel()->pSat();


					Pair<tmp<volScalarField>> vDotP = mixture.cavitationModel()
							->vDotP();
					const volScalarField& vDotcP = vDotP[0]();
					const volScalarField& vDotvP = vDotP[1]();


					// Cache p_rgh prior to solve for density update
					volScalarField p_rgh_0(p_rgh);

					while (pimple.correctNonOrthogonal())
					{

						fvScalarMatrix p_rghEqnIncomp(
								fvc::div(phiHbyA)
										- fvm::laplacian(rAUf, p_rgh)
										- (vDotvP - vDotcP) * (pSat - rho * gh)
										+ fvm::Sp(
										vDotvP - vDotcP,
										p_rgh));

						tmp<fvScalarMatrix> p_rghEqnComp;

						phasei = 0;
						forAllConstIter
						(
								PtrDictionary<phaseModel>,
								mixture.phases(),
								phase
						){
						tmp<fvScalarMatrix> hmm
						(
								(max(phase(), scalar(0))/phase().thermo().rho())
								*p_rghEqnComps[phasei]
						);

						if (phasei == 0)
						{
							p_rghEqnComp = hmm;
						}
						else
						{
							p_rghEqnComp.ref() += hmm;
						}

						phasei++;
					}

						solve(
								p_rghEqnComp + p_rghEqnIncomp,
								mesh.solver(
										p_rgh.select(pimple.finalInnerIter())));

						if (pimple.finalNonOrthogonalIter())
						{
							phasei = 0;
							forAllIter
							(
									PtrDictionary<phaseModel>,
									mixture.phases(),
									phase
							){
							phase().dgdt() =
							pos(phase())*
							(p_rghEqnComps[phasei] & p_rgh)/phase().thermo().rho();
						}

							phi = phiHbyA + p_rghEqnIncomp.flux();

							U = HbyA + rAU * fvc::reconstruct(
									(phig + p_rghEqnIncomp.flux()) / rAUf);
							U.correctBoundaryConditions();
						}



					}

					p = max(p_rgh + mixture.rho() * gh, pMin);

					mixture.correctRho(p_rgh - p_rgh_0);
					rho = mixture.rho();

					// Correct p_rgh for consistency with p and the updated densities
					p_rgh = p - rho * gh;
					p_rgh.correctBoundaryConditions();

					K = 0.5 * magSqr(U);

					Info << "max(U) " << max(mag(U)).value() << endl;
					Info << "min(p_rgh) " << min(p_rgh).value() << endl;
				}


			}

			if (pimple.turbCorr())
			{
				turbulence->correct();
			}
		}

		runTime.write();

		Info
				<< "ExecutionTime = "
				<< runTime.elapsedCpuTime()
				<< " s"
				<< "  ClockTime = "
				<< runTime.elapsedClockTime()
				<< " s"
				<< nl
				<< endl;
	}

	Info << "End\n" << endl;

	return 0;
}


// ************************************************************************* //
