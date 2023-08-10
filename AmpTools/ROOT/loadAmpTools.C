
{
  cout << "--------------------------------------------------" << endl;
  cout << "------      Loading AmpTools Modules        ------" << endl;
  cout << "--------------------------------------------------" << endl;

  gROOT->ProcessLine( ".include ${AMPTOOLS}" );

  cout << "Loading IUAmpTools/report.cc.." << endl;
  gROOT->LoadMacro( "IUAmpTools/report.cc+" );

  cout << "Loading IUAmpTools/ConfigurationInfo.cc.." << endl;
  gROOT->LoadMacro( "IUAmpTools/ConfigurationInfo.cc+" );

  cout << "Loading IUAmpTools/ConfigFileParser.cc.." << endl;
  gROOT->LoadMacro( "IUAmpTools/ConfigFileParser.cc+" );

  cout << "Loading IUAmpTools/NormIntInterface.cc.." << endl;
  gROOT->LoadMacro( "IUAmpTools/NormIntInterface.cc+" );

  cout << "Loading IUAmpTools/FitResults.cc.." << endl;
  gROOT->LoadMacro( "IUAmpTools/FitResults.cc+" );

  cout << "--------------------------------------------------" << endl;
  cout << "------  Finished Loading AmpTools Modules   ------" << endl;
  cout << "--------------------------------------------------" << endl;
}
