function printHeader()
%printHeader Summary of this function goes here
%   Detailed explanation goes here

disp('/*---------------------------------------------------------------------------*\');
disp('|           ____            |                                                 |');
disp('|          //               | Inspiration: solids4Foam / foam-extend / uFVM   |');
disp('|         //                | Version:    0.5.1                               |');
disp('|        // //\             | Web:        -/-                                 |');
disp('|       // //  \    n ano   | Thanks to:                                      |');
disp('|      // //    \   F inite |    CNPQ - Science without borders / Brazil      |');
disp('|     // //      \  V olume |    FNB - TU-Darmstadt / Germany                 |');
disp('|    // =========== M ethod |    GSC - TU-Darmstadt / Germany                 |');
disp('| __//                      | For copyright notice see file Copyright         |');
disp('\*---------------------------------------------------------------------------*/');

[ ~, hostname ] = system('hostname');
fprintf('Host     : %s', hostname);
fprintf('Date     : %s\n', date);
fprintf('Time     : %s\n', datestr(now, 'HH:MM:SS'));
fprintf('Case     : %s\n', pwd);
fprintf('\n// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n');

end

