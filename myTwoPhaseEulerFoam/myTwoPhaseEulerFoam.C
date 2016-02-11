/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

Application
    myTwoPhaseEulerFoam.C

Description
    currently a copy of twoPhaseEulerFoam only

    Solver for a system of 2 compressible fluid phases with one phase
    dispersed, e.g. gas bubbles in a liquid including heat-transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "twoPhaseSystem.H"
#include "PhaseCompressibleTurbulenceModel.H"
#include "pimpleControl.H"
#include "IOMRFZoneList.H"
#include "fvIOoptionList.H"
#include "fixedFluxPressureFvPatchScalarField.H"

//#include "myWaveTransmissiveFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "readGravitationalAcceleration.H"
    #include "createFields.H"
    #include "createMRFZones.H"
    #include "createFvOptions.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    #include "CourantNos.H"
    #include "setInitialDeltaT.H"
    
    //for testing only with a constant interpahse mass transfer term 	    	
    dimensionedScalar gamma_LV
 	("gamma_LV",
	 dimensionSet (1,-3,-1,0,0,0,0),
	 scalar(0.000));
	
    dimensionedScalar gamma_VL
 	("gamma_LV",
	 dimensionSet (1,-3,-1,0,0,0,0),
	 scalar(0.000));
    

    //label patchID = mesh.boundaryMesh().findPatchID("sideRight");//locate particular patch ID
    //const polyPatch& cPatch = mesh.boundaryMesh()[patchID];
    //patch().magSf()[patchID] = patch().magSf()[patchID]*scalar(0.1);				

    pimpleControl pimple(mesh);
      
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNos.H"
        #include "setDeltaT.H"

        runTime++;
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            fluid.solve();
            fluid.correct();

            volScalarField contErr1
            (
                fvc::ddt(alpha1, rho1) + fvc::div(alphaRhoPhi1)
              - (fvOptions(alpha1, rho1)&rho1)+gamma_LV-gamma_VL // add mass transfer source term here 
            );

            volScalarField contErr2
            (
                fvc::ddt(alpha2, rho2) + fvc::div(alphaRhoPhi2)
               - (fvOptions(alpha2, rho2)&rho2)-gamma_LV+gamma_VL // add mass transfer source term here 
            );


            #include "UEqns.H"
            #include "EEqns.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            #include "DDtU.H"

            if (pimple.turbCorr())
            {
                fluid.correctTurbulence();
            }
        }

        #include "write.H"

        Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
