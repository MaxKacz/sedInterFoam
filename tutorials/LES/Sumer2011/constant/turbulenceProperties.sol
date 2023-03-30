/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2.4.0                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "constant";
    object      turbulenceProperties;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

simulationType  LES;

LES{
     LESModel        partDynamicLagrangian;

     turbulence      on;

     printCoeffs     on;

     delta           cubeRootVol;

     filter          simple;

     CyaMax          0.3;

     nutaMax         1e-2;

     vanDriestCoeffs
     {
         delta           cubeRootVol;
         cubeRootVolCoeffs
         {
             deltaCoeff      1;
         }

         Aplus           26;
         Cdelta          0.158;
     }
}
// ************************************************************************* //
