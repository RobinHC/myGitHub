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

#include "phaseChangeTwoPhaseMixture.H"

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
    
    double celld; 
    
    double scaleFactor;

    double xDim;	

    // declare the cell centre variable 
    volScalarField flowArea // flow area of that specific cell 	
    (
     IOobject
     (

      "flowArea",
      runTime.timeName(),
      mesh

     ),

     mesh,
     dimensionedScalar("0",dimensionSet (0,2,0,0,0),0)
     
    );	
   				

    // declare the gradient of area along flow direction
    volScalarField flowAreaGrad 	
    (
     IOobject
     (

      "flowAreaGrad",
      runTime.timeName(),
      mesh

     ),

     mesh,
     dimensionedScalar("0",dimensionSet (0,1,0,0,0),0)
     
    );	
    
    // declare the area part of the source term
    volScalarField areaSource 	
    (
     IOobject
     (

      "areaSource",
      runTime.timeName(),
      mesh

     ),

     mesh,
     dimensionedScalar("0",dimensionSet (0,-1,0,0,0),0)
     
    );	
    // obtain cell length 
    const faceList & ff = mesh.faces();
    const pointField & pp = mesh.points();

    forAll( mesh.C(), celli)
    {
    	const cell & cc = mesh.cells()[celli];
	labelList pLabels(cc.labels(ff));
	pointField pLocal(pLabels.size(), vector::zero);

	forAll (pLabels, pointi)
	{
		pLocal[pointi] = pp[pLabels[pointi]];
	}

	xDim = Foam::max(pLocal & vector(1,0,0)) - Foam::min(pLocal & vector(1,0,0));

    }


    // assign flowarea value to the cell centre 
    forAll(flowArea,celli)
    {
    	flowArea[celli] = 1; // initialisation
    }

    forAll(flowArea,celli)
    {
    	if (celli >= 900 && celli < 1000)
	{
		celld = (double) celli;  
		scaleFactor = Foam::tanh((celld-899.0)/100.0);	
		flowArea[celli] = flowArea[899]-flowArea[899]*scalar(scaleFactor);	
	}
	else
	{
		flowArea[celli] = 1;
	}
    	
    }

    // assign flowAreaGrad to each cell centre
    forAll(flowAreaGrad, celli)
    {
	if (celli < 900)    
	{
		flowAreaGrad[celli] = 0; 
    	}
    	else if (celli >= 900 && celli < 1000) 
	{
		flowAreaGrad[celli] = (flowArea[celli+1] - flowArea[celli])/scalar(xDim);
	}
	else
	{	
		flowAreaGrad[celli] = 0;
	}	
    }
	
    // assign areaSource to each cell centre
    forAll(areaSource, celli)
    {
	if (celli < 900)    
	{
		areaSource[celli] = 0; 
    	}
    	else if (celli >= 900 && celli < 1000) 
	{
		areaSource[celli] = flowAreaGrad[celli]/flowArea[celli];
	}
	else
	{	
		areaSource[celli] = 0;
	}	
    }
    //label patchID = mesh.boundaryMesh().findPatchID("sideLeft");//locate particular patch ID

    //Info<< "patchID" << patchID << nl << endl; 
    
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
              - (fvOptions(alpha1, rho1)&rho1)                      //+gamma_LV-gamma_VL // fvOptions are the runtime semi-implicit source term 
              + rho1*mag(U1)*areaSource*alpha1	
	    );

            volScalarField contErr2
            (
                fvc::ddt(alpha2, rho2) + fvc::div(alphaRhoPhi2)
               - (fvOptions(alpha2, rho2)&rho2)                    //-gamma_LV+gamma_VL // 
               + rho2*mag(U2)*areaSource*alpha2	 
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
	
	//const volScalarField& test = alpha1_.db().lookupObject<volScalarField>("flowAreaGrad");


   	Info<< "flowArea=" << flowArea[999] << nl << endl; 
   	//Info<< "flowAreaGrad=" << flowAreaGrad[1000] << nl << endl; 
        //Info<< "areaSource=" << areaSource[1000] << nl << endl;
	
	Info<< "ExecutionTime = "
            << runTime.elapsedCpuTime()
            << " s\n\n" << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
