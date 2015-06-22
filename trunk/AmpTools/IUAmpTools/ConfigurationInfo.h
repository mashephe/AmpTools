#if !defined(CONFIGURATIONINFO)
#define CONFIGURATIONINFO

//******************************************************************************
// This file is part of AmpTools, a package for performing Amplitude Analysis
// 
// Copyright Trustees of Indiana University 2010, all rights reserved
// 
// This software written by Matthew Shepherd, Ryan Mitchell, and 
//                  Hrayr Matevosyan at Indiana University, Bloomington
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
// 1. Redistributions of source code must retain the above copyright
//    notice and author attribution, this list of conditions and the
//    following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright
//    notice and author attribution, this list of conditions and the
//    following disclaimer in the documentation and/or other materials
//    provided with the distribution.
// 3. Neither the name of the University nor the names of its contributors
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// 
// Creation of derivative forms of this software for commercial
// utilization may be subject to restriction; written permission may be
// obtained from the Trustees of Indiana University.
// 
// INDIANA UNIVERSITY AND THE AUTHORS MAKE NO REPRESENTATIONS OR WARRANTIES, 
// EXPRESS OR IMPLIED.  By way of example, but not limitation, INDIANA 
// UNIVERSITY MAKES NO REPRESENTATIONS OR WARRANTIES OF MERCANTABILITY OR 
// FITNESS FOR ANY PARTICULAR PURPOSE OR THAT THE USE OF THIS SOFTWARE OR 
// DOCUMENTATION WILL NOT INFRINGE ANY PATENTS, COPYRIGHTS, TRADEMARKS, 
// OR OTHER RIGHTS.  Neither Indiana University nor the authors shall be 
// held liable for any liability with respect to any claim by the user or 
// any other party arising from use of the program.
//******************************************************************************


#include <iostream>
#include <utility>
#include <string>
#include <vector>
#include <map>
#include <complex>
#include <fstream>

using namespace std;

class ConfigurationInfo;
class ReactionInfo;
class CoherentSumInfo;
class AmplitudeInfo;
class ParameterInfo;


/**
 *
 *    The ConfigurationInfo class holds four types of objects:
 *
 *      -#  REACTIONS (ReactionInfo Class)
 *            - Reactions are uniquely specified by a reaction name.
 *            - ReactionInfo must contain a vector of particle names.
 *                 (These are used by the Kinematics class, for example,
 *                  to keep track of four-vectors.)
 *            - It also holds instructions for how to read data and MC.
 *      -#  COHERENT SUMS (CoherentSumInfo Class)
 *            - Coherent sums are uniquely specified by a reaction name
 *                 and a sum name.
 *            - Every coherent sum must be associated with an existing reaction.
 *      -#  AMPLITUDES (AmplitudeInfo Class)
 *            - Amplitudes are uniquely specified by a reaction name,
 *                 a sum name, and an amplitude name
 *            - Every amplitude must be associated with an existing reaction 
 *                 and sum.
 *            - Every amplitude is a product of factors (specified by addFactor).
 *            - In every amplitude factor (vector<string>), the first element
 *                 is the amplitude class name, the other elements are passed into
 *                 the amplitude class as arguments.
 *      -#  PARAMETERS (ParameterInfo Class)
 *            - Parameters are uniquely specified by the parameter name.
 *            - Associate parameters with amplitudes using the 
 *                 AmplitudeInfo class (addParameter)
 *            - These parameters are in addition to production amplitude
 *                 parameters (which are implied for every amplitude)
 *
 *    Use the ConfigurationInfo class to create instances of the above
 *      objects with the createX methods.
 *
 *    Return lists of objects with the xList methods.
 *
 *    The purpose of the ConfigurationInfo class is to organize all information
 *    that is necessary to perform an amplitude analysis.  It can be setup
 *    directly by the user from the provided member functions or, alternatively,
 *    created by some tool.  Currently the most popular method for setting
 *    up the fit configuration is to utilize the provided ConfigFileParser
 *    to setup the ConfigurationInfo object from a user-supplied text file.
 *    ConfigurationInfo can then be used to configure the AmplitudeManager
 *    and the ParameterManager.
 *
 * \see reactionList
 * \see coherentSumList
 * \see amplitudeList
 * \see parameterList
 *
 * \ingroup IUAmpTools
 */


class ConfigurationInfo
{
public:
  
  /**
   * The constructor, which takes the name of the fit being done.
   */
  ConfigurationInfo(const string& fitName) { m_fitName = fitName; }
  ~ConfigurationInfo();
  
  /**
   * A function to retreive the fit name.
   */
  string fitName() const { return m_fitName; }


  /**
   * Name of the fit results file.
   */
  string fitOutputFileName() const { return m_fitName + ".fit"; }


  /**
   * A list of user-defined keywords.
   *
   * \see userKeywordArguments
   * \see addUserKeyword
   * \see removeUserKeyword
   */
  vector< string > userKeywords() const;

  /**
   * A list of arguments (in the form of vector<string>) associated with a 
   * given user-defined keyword.
   *
   * \see userKeywords
   * \see addUserKeyword
   * \see removeUserKeyword
   */
  vector< vector<string> > userKeywordArguments(const string& userKeyword) const;

  
  // Return lists of reactions, sums, amplitudes, or parameters.
  //   (If an argument to the call is empty, it is treated as
  //    a wildcard.  For example, reactionList("") returns all
  //    reactions.)
  
  /**
   * This returns a vector of all registered reactions (ReactionInfo) objects
   * when called without an argument.  With an argument it returns a vector
   * containg just those reactions that match the reaction name.  Note that
   * the list does not contain const pointers so the reactions themselves
   * may be modified.
   *
   * \param[in] reactionName (optional) the name of the reaction
   *
   * \see reaction
   */
  vector<ReactionInfo*>    reactionList    (const string& reactionName="") const;
  
  /**
   * This returns a vector of all coherent sums.  Optionally the user can
   * restrict the list to those registered to a particular reaction and
   * with a specific name.  Passing in an empty string to either argument
   * implies the search will not filter on that argument.
   *
   * \param[in] reactionName (optional) the name of the reaction
   * \param[in] sumName (optional) the name of the coherent sum
   *
   * \see coherentSum
   */
  vector<CoherentSumInfo*> coherentSumList (const string& reactionName="",
                                            const string& sumName="") const;
  
  /**
   * This returns a vector of all amplitudes.  Optionally the user can
   * restrict the list to those associated with a particular reaction,
   * coherent sum, or having a particular name.  Passing in an empty string 
   * to any argumentimplies the search will not filter on that argument.
   *
   * \param[in] reactionName (optional) the name of the reaction
   * \param[in] sumName (optional) the name of the sum
   * \param[in] ampName (optional) the name of the amplitude
   *
   * \see amplitude
   */
  vector<AmplitudeInfo*>   amplitudeList   (const string& reactionName="",
                                            const string& sumName="",
                                            const string& ampName="") const;
  
  /**
   * This returns a vector of all parameters.  The user can optionally
   * restrict the list to those associated with a particular reaction,
   * coherent sum, amplitude, or having a particular name.  Passing in an 
   * empty string to any argumentimplies the search will not filter on 
   * that argument.
   *
   * \param[in] reactionName (optional) the name of the reaction
   * \param[in] sumName (optional) the name of the sum
   * \param[in] ampName (optional) the name of the amplitude
   * \param[in] parName (optional) the name of the parameter
   * 
   * \see parameter
   */  
  vector<ParameterInfo*>   parameterList   (const string& reactionName="",
                                            const string& sumName="",
                                            const string& ampName="",
                                            const string& parName="") const;
  
  
  // Return a specific reaction, sum, amplitude, or parameter.
  
  /**
   * Similar to reactionList above but returns a pointer to a specific
   * reaction.  Note that the pointer is not const.  This routine can
   * be used to modify a ReactionInfo object in the configuration.
   *
   * \param[in] reactionName the name of the reaction
   *
   * \see reactionList
   */
  ReactionInfo*    reaction    (const string& reactionName) const;
  
  /**
   * Similar to coherentSumList above but returns a pointer to a specific
   * sum.  Note that the pointer is not const.  This routine can
   * be used to modify a CoherentSumInfo object in the configuration.
   *
   * \param[in] reactionName the name of the reaction
   * \param[in] sumName the name of the sum
   *
   * \see coherentSumList
   */
  CoherentSumInfo* coherentSum (const string& reactionName,
                                const string& sumName) const;
  
  /**
   * Similar to amplitudeList above but returns a pointer to a specific
   * amplitude.  Note that the pointer is not const.  This routine can
   * be used to modify a CoherentSumInfo object in the configuration.
   *
   * \param[in] reactionName the name of the reaction
   * \param[in] sumName the name of the sum
   * \param[in] ampName the name of the amplitude
   *
   * \see amplitudeList
   */
  AmplitudeInfo*   amplitude   (const string& reactionName,
                                const string& sumName,
                                const string& ampName) const;
  
  /**
   * Similar to parameterList above but returns a pointer to a specific
   * parameter.  Note that the pointer is not const.  This routine can
   * be used to modify a ParameterInfo object in the configuration.
   *
   * \param[in] parName the name of the parameter
   *
   * \see parameterList
   */
  ParameterInfo*   parameter   (const string& parName) const;
  
  
  
  // Create a new reaction, sum, amplitude or parameter.
  //   (If the object already exists, the old is replaced and
  //    its contents are cleared.)
  
  /**
   * A routine to create a new reaction.  This returns a pointer to the
   * reaction that has been created so it can be further modified.
   * If the reaction already exists, it will be overwritten with a new
   * reaction.
   *
   * \param[in] reactionName the name of the reaciton to be created
   * \param[in] particleList a vector of the strings that identify the
   * particles in the reaction
   *
   * \see removeReaction
  */
  ReactionInfo*    createReaction    (const string&         reactionName, 
                                      const vector<string>& particleList);
  
  /**
   * A routine to create a new coherent sum.  This returns a pointer to
   * the coherent sum that has been created so it can be further modified.
   * If the coherent sum already exists, it will be overwritten with a new
   * coherent sum.
   *
   * \param[in] reactionName the name of the reaction for which to create
   * the coherent sum
   * \param[in] sumName the name of the sum to create
   *
   * \see removeCoherentSum
   */
  CoherentSumInfo* createCoherentSum (const string& reactionName, 
                                      const string& sumName);
  
  /**
   * A routine to create a new amplitude.  This returns a pointer to the
   * amplitude that has been created so it can be further modified.
   * If the amplitude already exists, it will be overwritten with a new
   * amplitude.
   *
   * \param[in] reactionName the name of the reaction for which to create
   * the amplitude
   * \param[in] sumName the coherent sum that the amplitude should be
   * added to
   * \param[in] ampName the name of amplitude to be created
   *
   * \see removeAmplitude
   */
  AmplitudeInfo*   createAmplitude   (const string& reactionName, 
                                      const string& sumName, 
                                      const string& ampName);
  
  /**
   * This creates a new parameter and returns a pointer to the ParameterInfo
   * object for that parameter so it can be further modifed.  If the parameter
   * already exists, it will be overwritten with a new parameter.
   *
   * \param[in] parameterName the name of the parameter to create
   * \param[in] value the initial value of the parameter
   *
   * \see removeParameter
   */
  ParameterInfo*   createParameter   (const string& parName,
                                      double        value);
  
  
  // Remove reactions, sums, amplitudes, or parameters.
  //   (Blank arguments are treated like wildcards.
  //    This also removes associated objects -- for example,
  //    removeReaction("") removes all reactions but also all sums
  //    and amplitudes associated with those reactions.
  
  /**
   * This removes a reaction or reactions that are matched to the arguments.
   * Passing in a null string to a particular argument will match all
   * instances of that argument (like a wildcard).  The null string is
   * the default argument for all parameters.
   *
   * \param[in] reactionName the name of the reaciton
   *
   * \see createReaction
   */
  void removeReaction    (const string& reactionName="");
  
  /**
   * This removes a coherent sum or sums that are matched to the arguments.
   * Passing in a null string to a particular argument will match all
   * instances of that argument (like a wildcard).  The null string is
   * the default argument for all parameters.
   *
   * \param[in] reactionName the name of the reaciton
   * \param[in] sumName the name of the sum
   *
   * \see createCoherentSum
   */
  void removeCoherentSum (const string& reactionName="",
                          const string& sumName="");

  /**
   * This removes an amplitude or amplitudes that are matched to the arguments.
   * Passing in a null string to a particular argument will match all
   * instances of that argument (like a wildcard).  The null string is
   * the default argument for all parameters.
   *
   * \param[in] reactionName the name of the reaciton
   * \param[in] sumName the name of the sum
   * \param[in] ampName the name of the amplitude
   *
   * \see createAmplitude
   */
  void removeAmplitude   (const string& reactionName="",
                          const string& sumName="",
                          const string& ampName="");
  
  /**
   * This removes an a parameter or parameters matched to the arguments.
   * Passing in a null string to a particular argument will match all
   * instances of that argument (like a wildcard).  The null string is
   * the default argument for all parameters
   *
   * \param[in] parName the name of the parameter
   *
   * \see createParameter
   */
  void removeParameter   (const string& parName="");


  /**
   * Set the name of a fit.
   */
  void setFitName (const string& fitName) { m_fitName = fitName; }  


  /**
   * Add a user-defined keyword and its arguments.
   * A user-defined keyword can have multiple sets of arguments.
   *
   * \see userKeywords
   * \see userKeywordArguments
   * \see removeUserKeyword
   */
  void addUserKeyword(const string& uesrKeyword, const vector<string>& arguments);  


  /**
   * Remove a user-defined keyword.
   * If userKeyword is not specified, all keywords are removed.
   *
   * \see userKeywords
   * \see userKeywordArguments
   * \see addUserKeyword
   */
  void removeUserKeyword(const string& userKeyword="");
    
  /**
   * This displays a nicely formatted description of the configuration to
   * the screen (when no arguments are specified) or prints to a file named
   * by the first argument.  If the second argument is set to true the 
   * configuration info will be appeneded to the file.
   *
   * \param[in] fileName (optional) name of file to write info
   * \param[in] append (optional) set to true to append to the file, default
   * is to overwrite the file
   */
  void display(string fileName = "", bool append = false) const;
  void display(string fileName = "", bool append = false);


  /**
   * This writes the current status of this ConfigurationInfo to a file.
   * The output can be parsed by a ConfigFileParser object to recreate
   * this ConfigurationInfo object.
   *
   * \see ConfigFileParser
   */
  void write( const string& fileName ) const;
  ostream& write( ostream& output ) const;
    
  /**
   * This produces a map where the keys are the names of the amplitudes and
   * the values are vectors of the names of the amplitudes whose production
   * parameters are constrained to be the same as the key amplitude's 
   * production parameter.
   */
  map< string, vector< string > > constraintMap() const;
  
private:
  
  string                   m_fitName;
  vector<ReactionInfo*>    m_reactions;
  vector<CoherentSumInfo*> m_sums;
  vector<AmplitudeInfo*>   m_amplitudes;
  vector<ParameterInfo*>   m_parameters;
  map<string, vector< vector<string> > > m_userKeywordMap;
    
};


inline ostream& operator<<( ostream& output, const ConfigurationInfo& cfgInfo ){
  
  return cfgInfo.write( output );
}


/**
 * The reaction info holds all information that is specific to a particular
 * reaction or final state.  It is typically characterized by a set of
 * unique observable particles in the detector.  In order to create a reaction
 * it must be given a name that is unique in the fit and also a list of strings
 * that identify the particles.  
 *
 * \ingroup IUAmpTools
 */

class ReactionInfo
{
public:
  
  /**
   * The constructor.
   *
   * \param[in] reactionName a unique name to identify this reaction
   * \param[in] particleList a list of particles that define the reaction
   */
  ReactionInfo(const string&         reactionName,
               const vector<string>& particleList):
  m_reactionName(reactionName),
  m_particleList(particleList)  {clear();};
  
  
  // Every reaction has a reactionName and a particleList
  
  /**
   * This returns the unique name of the reaction.
   */
  string            reactionName()    const {return m_reactionName;}

  /**
   * This returns the list of particles that characterize the reaction.
   *
   * \see setParticleList
   */
  const vector<string>&    particleList()    const {return m_particleList;}
  
  
  // Return information about data and MC associated with this reaction
  
  /**
   * Returns the name of a DataReader and a vector of arguments (e.g. file names)
   * to that DataReader for DATA associated with this reaction.
   *
   * \see setData
   */
  const pair< string, vector<string> >& data()  const {return m_data;}

  /**
   * Returns the name of a DataReader and a vector of arguments (e.g. file names)
   * to that DataReader for BACKGROUND associated with this reaction.
   *
   * \see setBkgnd
   */
  const pair< string, vector<string> >& bkgnd()  const {return m_bkgnd;}

  /**
   * Returns the name of a DataReader and a vector of arguments (e.g. file names)
   * to that DataReader for GENERATED MC associated with this reaction.
   *
   * \see setGenMC
   */
  const pair< string, vector<string> >& genMC()  const {return m_genMC;}

  /**
   * Returns the name of a DataReader and a vector of arguments (e.g. file names)
   * to that DataReader for ACCEPTED MC associated with this reaction.
   *
   * \see setAccMC
   */
  const pair< string, vector<string> >& accMC()  const {return m_accMC;}

  /**
   * Returns the name of the file that the normalization integral cache
   * is written out to.
   *
   * \see setNormIntFile
   */
  string            normIntFile()     const {return m_normIntFile;}

  /**
   * Returns a boolean that indicates whether this normIntFile should be
   * used as input (true) or not (false).
   *
   * \see setNormIntFile
   */
  bool       normIntFileInput()     const {return m_normIntFileInput;}
  
  
  // Display or clear information for this reaction
  
  /**
   * Displays information about this reaction to the screen or writes it to
   * file if a file name is passed in.
   * 
   * \param[in] fileName (optional) name of file to write reaction info to
   * \param[in] append (optional) if true (default) information will be 
   * appended to the file, if false is passed in the file will be overwritten
   */
  void              display(string fileName = "", bool append = true);

  /**
   * Clears all internal data for this reaction
   */
  void              clear();
  
  
  // Use these methods to add new information to this reaction
  
  /**
   * Sets the vector of strings that provides the names of the final state particles.
   *
   * \param[in] particleList vector of strings of final state particles
   *
   * \see particleList
   */
  void  setParticleList (const vector<string>& particleList) {m_particleList = particleList;}

  /**
   * Sets the name of a DataReader and a vector of arguments (e.g. file names)
   * to that DataReader for DATA associated with this reaction.
   *
   * \param[in] classname the name of the DataReader that will read DATA
   * \param[in] args arguments to the DataReader
   *
   * \see data
   */
  void  setData   (const string& classname, const vector<string>& args)
                            { m_data = pair<string, vector<string> >(classname,args); }
  
  /**
   * Sets the name of a DataReader and a vector of arguments (e.g. file names)
   * to that DataReader for BACKGROUND associated with this reaction.
   *
   * \param[in] classname the name of the DataReader that will read DATA
   * \param[in] args arguments to the DataReader
   *
   * \see data
   */
  void  setBkgnd  (const string& classname, const vector<string>& args)
  { m_bkgnd = pair<string, vector<string> >(classname,args); }

  /**
   * Sets the name of a DataReader and a vector of arguments (e.g. file names)
   * to that DataReader for GENERATED MC associated with this reaction.
   *
   * \param[in] classname the name of the DataReader that will read GENERATED MC
   * \param[in] args arguments to the DataReader
   *
   * \see genMC
   */
  void  setGenMC   (const string& classname, const vector<string>& args)
                            { m_genMC = pair<string, vector<string> >(classname,args); }

  /**
   * Sets the name of a DataReader and a vector of arguments (e.g. file names)
   * to that DataReader for ACCEPTED MC associated with this reaction.
   *
   * \param[in] classname the name of the DataReader that will read ACCEPTED MC
   * \param[in] args arguments to the DataReader
   *
   * \see accMC
   */
  void  setAccMC   (const string& classname, const vector<string>& args)
                            { m_accMC = pair<string, vector<string> >(classname,args); }

  /**
   * Sets the name of the file to use as the normalization integral cache
   *
   * \param[in] normIntFile name of the normalization integral file
   * \param[in] input use this normIntFile as input
   *
   * \see normIntFile
   */
  void  setNormIntFile  (const string& normIntFile, bool input = false) 
                           { m_normIntFile = normIntFile;
                             m_normIntFileInput = input; }
  
  
private:
  
  string                         m_reactionName;
  vector<string>                 m_particleList;
  pair< string, vector<string> > m_data;
  pair< string, vector<string> > m_bkgnd;
  pair< string, vector<string> > m_genMC;
  pair< string, vector<string> > m_accMC;
  string                         m_normIntFile;
  bool                           m_normIntFileInput;
  
};


/**
 * This class holds all information for a specific coherent sum inside of
 * a final state.  Each coherent sum is defined by the reaction that it is
 * associated with and its name.  All sums for a particular reaction must
 * have unique names.  In the construction of the intensity for a particular
 * reaction, all amplitudes within a coherent sum are added coherently
 * (with interference) \f$ \left| \sum_i V_i A_i \right|^2 \f$ for amplitudes
 * and production parameters \f$ A_i \f$ and \f$ V_i \f$.  The individual
 * coherent sums within a reaction are then added incoherently, e.g. 
 * \f$ \sum_\beta \left| \sum_i V_{i,\beta} A_{i,\beta} \right|^2 \f$, where
 * \f$ \beta \f$ indexes the coherent sums and i indexes the amplitudes.
 *
 * \ingroup IUAmpTools
 */

class CoherentSumInfo
{
public:
  
  /**
   * The constructor for the class
   *
   * \param[in] reactionName the name of the reaction to associate the sum with
   * \param[in] sumName the name of the sum, must be unique to the reaction
   */
  CoherentSumInfo(const string& reactionName,
                  const string& sumName):
  m_reactionName(reactionName),
  m_sumName(sumName)  {clear();};
  
  
  // Every sum is associated with a reaction and has a sumName
  
  /**
   * This returns the name of the reaction that the sum is associated with.
   */
  string            reactionName()    const {return m_reactionName;}

  /**
   * This returns the name of the sum.
   */
  string            sumName()         const {return m_sumName;}

  /**
   * This returns the "full name" of the sum.  It is a string formed by
   * concatenating the reaction name with the sum num name using a double
   * colon:  reactionName::sumName
   */
  
  string            fullName()        const {return (m_reactionName + "::" +
                                                     m_sumName);}
  
  //  Display information for this sum
  
  /**
   * Displays information about this sum to the screen or writes it to
   * file if a file name is passed in.
   * 
   * \param[in] fileName (optional) name of file to write sum info to
   * \param[in] append (optional) if true (default) information will be 
   * appended to the file, if false is passed in the file will be overwritten
   */  
  void              display(string fileName = "", bool append = true);

  /**
   * Currently does nothing -- any future cleanup needed to return a 
   * coherent sum to its freshly-created state shoudl be done here.
   */
  void              clear() {}
  
private:
  
  string          m_reactionName;
  string          m_sumName;
  
};


/**
 * This class holds all information related to a single amplitude.  In the
 * construction of the intensity 
 * \f$ \sum_\beta \left| \sum_i V_{i,\beta} A_{i,\beta} \right|^2 \f$, where
 * \f$ \beta \f$ indexes the coherent sums and i indexes the amplitudes, this
 * class defines the \f$ A_{i,\beta} \f$.  The class is typically composed
 * of a list of user-defined factors that are multiplied together to
 * create the amplitude.
 *
 * \ingroup IUAmpTools
 */

class AmplitudeInfo
{
public:
  
  /**
   * The constructor:  each amplitude is uniquely defined by a reaction,
   * a coherent sum, and an amplitude name.
   *
   * \param[in] reactionName the name of the reaction
   * \param[in] sumName the name of the coherent sum the amplitude belongs to
   * \param[in] ampName the name of the amplitude itself
   */
  
  AmplitudeInfo(const string& reactionName,
                const string& sumName,
                const string& ampName):
  m_reactionName(reactionName),
  m_sumName(sumName),
  m_ampName(ampName)  {clear();};
  
  
  //  Every amplitude is identified by a reaction, a sum, and an ampName
  
  /**
   * This returns the name of the reaction for this amplitude.
   */
  string                     reactionName()   const {return m_reactionName;}

  /**
   * This returns the name of the sum for this amplitude.
   */
  string                     sumName()        const {return m_sumName;}
  
  /**
   * This returns the name of the amplitude.
   */
  string                     ampName()        const {return m_ampName;}
  
  /**
   * This returns the "full name" of the reaction.  It is obtained by concatenating
   * togther the reaction name, sum name, and amplitude name with a double colon,
   * e.g. reactionName::sumName::ampName
   */
  string                     fullName()       const {return (m_reactionName + "::" +
                                                             m_sumName + "::" +
                                                             m_ampName);}
  
  //  Return the amplitude factors for this reaction.
  //    Each factor has the form: <amplitude class name> (arg1) (arg2) ...
  
  /**
   * This returns information about all of the user-defined factors that make
   * up the amplitude.  It is a list of vectors -- each item in the initial
   * vector is a vector of strings.  The first string in the vector specifies
   * the name of the user-defined amplitude class that tells how to compute
   * the amplitude factor.  The remaining items in the vector of strings are
   * string arguments that can be passed into the newAmplitude routine to
   * create a new instance of the users amplitude.  The length of the returned
   * vector is the number of factors in the amplitude.
   *
   * \see Amplitude::newAmplitude
   * \see AmplitudeManager::addAmpFactor
   */
  const vector< vector<string> >&   factors()        const {return m_factors;}
  
  
  //  Return a list of different particle permuations.
  //    (For example, if there are 5 particles in this reaction, specifying
  //     a permutation of 01243 will permute the fourth and fifth particles and the
  //     amplitude becomes A(01234) + A(01243).)
  
  /**
   * This returns a vector of the different permutations for this amplitude.
   * It is used to store additional permuations besides the default permutation
   * and is used to setup the AmplitudeManager.
   *
   * \see AmplitudeManager::addAmpPermutation
   */
  
  const vector< vector<int> >&      permutations()   const {return m_permutations;}
  
  
  //  Return a list of amplitudes that are constrained to have the same
  //    production parameters as this amplitude.
  
  /**
   * This returns a vector that contains AmplitudeInfo pointers for every other
   * amplitude that is constrained to have the same production parameter
   * as the current amplitude.
   */
  const vector< AmplitudeInfo* >&   constraints()    const {return m_constraints;}


  /**
   * Check to see if this amplitude is constrained to another amplitude.
   *
   * \param[in] constraint pointer to an AmplitudeInfo object
   */
  bool hasConstraint(AmplitudeInfo* constraint) const;
  
  
  //  Return information about the production parameter for this amplitude.
  
  /**
   * This returns the current value of the production parameter.
   */
  complex< double >          value()          const {return m_value;}

  bool                       real()           const {return m_real;}
  bool                       fixed()          const {return m_fixed;}
  
  string                     scale()          const {return m_scale;} 


  //  Return a list of parameters (other than the production parameter)
  //    associated with this amplitude.
 
  /**
   * This returns a vector of pointers to all of the ParameterInfo objects
   * that are associated with this amplitude.  
   */
  const vector < ParameterInfo* >&  parameters()     const {return m_parameters;}


  //  Display or clear information for this amplitude

  /**
   * Displays information about this amplitude to the screen or writes it to
   * file if a file name is passed in.
   * 
   * \param[in] fileName (optional) name of file to write amplitude info to
   * \param[in] append (optional) if true (default) information will be 
   * appended to the file, if false is passed in the file will be overwritten
   */  
  void                       display(string fileName = "", bool append = true);

  /**
   * This clears out all the internal data for this particular AmplitudeInfo
   * object.
   */
  void                       clear();
  
  
  //  Use these methods to set up this ampltiude
  
  /**
   * This adds a new factor, i.e. a user-defined amplitude computing
   * routine, to the amplitude.
   *
   * \param[in] factor a vector of strings describing the factor.  The first
   * element of the vector is the name of the user's amplitude routine.  The
   * remaining elements are optional string arguments that are passed to
   * the newAmplitude routine of the user-defined Amplitude
   *
   * \see Amplitude::newAmplitude
   */
  void addFactor         (const vector<string>& factor)    {m_factors.push_back(factor);}
  
  /**
   * This adds a new permutation to the amplitude.  See documentation in the
   * AmplitudeManager below for info about permutations.
   *
   * \see AmplitudeManager::addAmpPermutation
   */
  void addPermutation    (const vector<int>& permutation)  {m_permutations.push_back(permutation);}
  
  /**
   * This adds the constraint that the production parameter of the current
   * amplitude must be the same as the production parameter for the amplitude
   * whose AmplitudeInfo is passed in.  If the argument's list of constraints
   * is not empty it creates constraints between those amplitudes and the
   * current amplitude also.
   *
   * \param[in] constraint a pointer to the AmplitudeInfo to which to constrain
   * the current amplitude's production parameter
   */
  void addConstraint     (AmplitudeInfo* constraint);

  /**
   * This sets the initial value of the production parameter for a particular
   * parameter.
   *
   * \param[in] value the desired initial value
   */
  void setValue          (complex<double> value)   {m_value = value;}
  
  void setReal           (bool real)               {m_real = real;}
  void setFixed          (bool fixed)              {m_fixed = fixed;}
  
  void setScale          (string scale)            {m_scale = scale;}


  /**
   * This associates some amplitude parameter described by the ParameterInfo
   * objects with the current amplitude.
   *
   * \param[in] parameter pointer the ParameterInfo object for the parameter
   */
  void addParameter      (ParameterInfo* parameter);
  
  
  //  Use these methods to remove information associated with this amplitude
  
  /**
   * This removes a particular AmplitudeInfo from the list of constraints.
   * It also removes itself from the list of constraints of each of the
   * constrained Amplitudes.
   *
   * \param[in] constraint a pointer to the AmplitudeInfo to be removed from
   * the list of constraints
   */
  void removeConstraint (AmplitudeInfo* constraint);

  /**
   * This removes an associated ParameterInfo pointer from the list of 
   * ParameterInfo objects associated with this amplitude.
   *
   * \param[in] parameter a pointer to the ParameterInfo object to removed
   */
  void removeParameter  (ParameterInfo* parameter);
  
private:
  
  string                       m_reactionName;
  string                       m_sumName;
  string                       m_ampName;
  vector< vector<string> >     m_factors;
  vector< vector<int> >        m_permutations;
  vector< AmplitudeInfo* >     m_constraints;
  complex< double >            m_value;
  bool                         m_real;
  bool                         m_fixed;
  string                       m_scale;
  vector< ParameterInfo* >     m_parameters;
  
};





class ParameterInfo
{
public:
  
  ParameterInfo(const string& parName,
                double        value):
  m_parName(parName),
  m_value(value) {clear();};
  
  
  // Every parameter has a parameter name and an initial value
  
  string parName()         const  {return m_parName;}
  double value()           const  {return m_value;}
  
  
  // If true, the parameter is fixed to its initial value.
  
  bool   fixed()           const  {return m_fixed;}
  
  
  // If the parameter is bounded, it is bounded between its
  //   lowerBound and upperBound.
  
  bool   bounded()         const  {return m_bounded;}
  double lowerBound()      const  {return m_lowerBound;}
  double upperBound()      const  {return m_upperBound;}
  
  
  // If the parameter is "gaussianBounded," an additional likelihood
  //   factor is created based on the parameter's position within a 
  //   a gaussian distribution with mean centralValue and error
  //   gaussianError.
  
  bool   gaussianBounded() const  {return m_gaussianBounded;}
  double centralValue()    const  {return m_centralValue;}
  double gaussianError()   const  {return m_gaussianError;}
  
  
  // Display or clear inforamtion for this parameter
  
  void display(string fileName = "", bool append = true);
  void clear();
  
  
  // Use these methods to set up this parameter
  
  void setValue           (double value)         {m_value = value;}
  void setFixed           (bool fixed)           {m_fixed = fixed;}
  void setBounded         (bool bounded)         {m_bounded = bounded;}
  void setLowerBound      (double lowerBound)    {m_lowerBound = lowerBound;}
  void setUpperBound      (double upperBound)    {m_upperBound = upperBound;}
  void setGaussianBounded (bool gaussianBounded) {m_gaussianBounded = gaussianBounded;}
  void setCentralValue    (double centralValue)  {m_centralValue = centralValue;}
  void setGaussianError   (double gaussianError) {m_gaussianError = gaussianError;}
  
  
private:
  
  string m_parName;
  double m_value;
  bool   m_fixed;
  bool   m_bounded;
  double m_lowerBound;
  double m_upperBound;
  bool   m_gaussianBounded;
  double m_centralValue;
  double m_gaussianError;
  
};




#endif
