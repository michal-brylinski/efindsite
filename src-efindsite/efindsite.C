/*
===============================================================================
         ___________.__            .____________.__  __          
     ____\_   _____/|__| ____    __| _/   _____/|__|/  |_  ____  
   _/ __ \|    __)  |  |/    \  / __ |\_____  \ |  \   __\/ __ \ 
   \  ___/|     \   |  |   |  \/ /_/ |/        \|  ||  | \  ___/ 
    \___  >___  /   |__|___|  /\____ /_______  /|__||__|  \___  >
        \/    \/            \/      \/       \/               \/ 

                                                  
   eFindSite - ligand-binding site prediction from meta-threading

   Computational Systems Biology Group
   Department of Biological Sciences
   Center for Computation & Technology
   Louisiana State University
   407 Choppin Hall, Baton Rouge, LA 70803, USA

   http://www.brylinski.org

   Report bugs to michal@brylinski.org

   Copyright 2013 Michal Brylinski

   This file is part of eFindSite.

   eFindSite is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   eFindSite is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with eFindSite. If not, see <http://www.gnu.org/licenses/>.

===============================================================================
*/


#include "efindsite.h"

using namespace std;

int main(int argc, char *argv[])
{
 time_t t_start, t_end, t_bench1, t_bench2;
 time(&t_start);
 
 double       cut_tmscore   = 0.40;
 double       cut_seqid     = 1.00;
 double       cut_probab    = 0.50;
 double       cut_clustdis  = 6.50;
 unsigned int cut_templates = MAXTPL;
 double       cut_binrest   = 0.2;
 int          cut_binresn   = 3;
 double       cut_clustlig  = 0.70;
 std::string  met_clustlig  = "T";
 std::string  met_clustdis  = "D";
 std::string  met_druggabl  = "R";
 double       cut_druggabl  = 0.0;
 
 cout << "------------------------------------------------------------" << endl
      << "                         efindsite" << endl
      << "                        version 1.3" << endl
      << "              ligand binding pocket prediction" << endl << endl
      << "       report bugs and issues to michal@brylinski.org" << endl
      << "------------------------------------------------------------" << endl << endl;
 
 if ( argc < 7 )
 {
  cout << " efindsite -s <target structure in PDB format>" << endl
       << "           -t <templates detected by eThread>" << endl
       << "           -e <sequence profile>" << endl
       << "           -o <output filename>" << endl << endl
       << " additional options:" << endl << endl
       << "           -l <auxiliary ligands in SDF format>" << endl
       << "           -b <sequence identity threshold (default 1.0)>" << endl
       << "           -p <eThread probability threshold (default 0.5)>" << endl
       << "           -m <TMscore threshold (default 0.4)>" << endl
       << "           -x <max number of templates (default " << MAXTPL << ")>" << endl
       << "           -r <binding residue threshold (default 0.2)>" << endl
       << "           -n <min number of binding residues (default 3)>" << endl
       << "           -g <pocket clustering method (default D)>" << endl
       << "               D - DBSCAN" << endl
       << "               L - average linkage" << endl
       << "           -d <pocket clustering cutoff (default 6.5)>" << endl
       << "           -c <fingerprint clustering method (default T)>" << endl
       << "               T - classical Tanimoto coeff" << endl
       << "               A - average Tanimoto coeff" << endl
       << "           -f <fingerprint clustering cutoff (default 0.7 for T and A)>" << endl
       << "           -u <druggability model (default R)>" << endl
       << "               R - logistic regression" << endl
       << "               D - linear discriminant analysis" << endl
       << "           -y <druggability cutoff (default 0.6)>" << endl << endl;
  
  exit(EXIT_SUCCESS);
 }
 
 string target_name;
 bool target_opt = false;
 
 string sequence_name;
 bool sequence_opt = false;
 
 string templates_name;
 bool templates_opt = false;
 
 string output_name;
 bool output_opt = false;
 
 string cmps_name;
 bool cmps_opt = false;
 
 bool drug_opt = false;
 
 for ( int i = 0; i < argc; i++ )
 {
  if ( !strcmp(argv[i],"-s") && i < argc ) { target_name    = string(argv[i+1]); target_opt    = true; }
  if ( !strcmp(argv[i],"-e") && i < argc ) { sequence_name  = string(argv[i+1]); sequence_opt  = true; }
  if ( !strcmp(argv[i],"-t") && i < argc ) { templates_name = string(argv[i+1]); templates_opt = true; }
  if ( !strcmp(argv[i],"-o") && i < argc ) { output_name    = string(argv[i+1]); output_opt    = true; }
  if ( !strcmp(argv[i],"-l") && i < argc ) { cmps_name      = string(argv[i+1]); cmps_opt      = true; }
  if ( !strcmp(argv[i],"-b") && i < argc ) { cut_seqid      = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-p") && i < argc ) { cut_probab     = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-m") && i < argc ) { cut_tmscore    = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-x") && i < argc ) { cut_templates  = atoi(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-r") && i < argc ) { cut_binrest    = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-n") && i < argc ) { cut_binresn    = atoi(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-g") && i < argc ) { met_clustdis   = string(argv[i+1]);                       }
  if ( !strcmp(argv[i],"-f") && i < argc ) { cut_clustlig   = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-c") && i < argc ) { met_clustlig   = string(argv[i+1]);                       }
  if ( !strcmp(argv[i],"-d") && i < argc ) { cut_clustdis   = atof(argv[i+1]);                         }
  if ( !strcmp(argv[i],"-u") && i < argc ) { met_druggabl   = string(argv[i+1]);                       }
  if ( !strcmp(argv[i],"-y") && i < argc ) { cut_druggabl   = atof(argv[i+1]);;  drug_opt      = true; }
 }
 
 char * path1;
 
 path1 = getenv("EF_LIB"); if ( path1==NULL ) { cout << "EF_LIB is not set" << endl; exit(EXIT_FAILURE); }
 
 path1 = getenv("EF_MOD"); if ( path1==NULL ) { cout << "EF_MOD is not set" << endl; exit(EXIT_FAILURE); }
 
 string lib_path;
 lib_path = getenv("EF_LIB");
 
 string model_path;
 model_path = getenv("EF_MOD");
 
 ifstream f01( (model_path+"/residues.model").c_str() );
 ifstream f02( (model_path+"/residues.scale").c_str() );
 
 ifstream f03( (model_path+"/ligands.model").c_str() );
 ifstream f04( (model_path+"/ligands.scale").c_str() );
 
 ifstream f05( (model_path+"/conf1.model").c_str() );
 ifstream f06( (model_path+"/conf1.scale").c_str() );
 
 ifstream f07( (model_path+"/conf2.model").c_str() );
 ifstream f08( (model_path+"/conf2.scale").c_str() );
 
 if ( !f01 || !f02 || !f03 || !f04 || !f05 || !f06 || !f07 || !f08 )
 {
  cout << "Could not find SVM models in " << model_path << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !target_opt )
 {
  cout << "Provide target structure in PDB format" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !sequence_opt )
 {
  cout << "Provide sequence profile" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !templates_opt )
 {
  cout << "Provide templates detected by eThread" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( !output_opt )
 {
  cout << "Provide output filename" << endl;
  exit(EXIT_FAILURE);
 }
 else
 {
  ofstream outprot( (output_name+".templates.pdb").c_str() );
  outprot.close();
  
  ofstream outali( (output_name+".alignments.dat").c_str() );
  outali.close();
  
  ofstream outlig( (output_name+".ligands.sdf").c_str() );
  outlig.close();
  
  ofstream outpkt1( (output_name+".pockets.pdb").c_str() );
  outpkt1.close();
  
  ofstream outpkt2( (output_name+".pockets.dat").c_str() );
  outpkt2.close();
 }
 
 if ( cut_seqid < 1.0 )
  cout << "!!! Benchmarking mode activated with max sid of " << cut_seqid << " !!!" << endl << endl;
 
 if ( cut_tmscore < 0.4 )
  cout << "!!! TMscore of " << cut_tmscore << " is below the statistical significance threshold !!!" << endl << endl;
 
 if ( cut_binresn < 1 )
 {
  cout << "!!! Min number of binding residues must be >0, setting to 1 !!!" << endl << endl;
  
  cut_binresn = 1;
 }
 
 if ( cut_probab < ( 1 / 1e6 ) )
 {
  cout << "!!! Threshold for eThread probability must be >0, setting to 0.5 !!!" << endl << endl;
  
  cut_probab = 0.5;
 }
 else if ( cut_probab > 1 )
 {
  cout << "!!! Threshold for eThread probability must be <=1, setting to 0.5 !!!" << endl << endl;
  
  cut_probab = 0.5;
 }
 
 if ( cut_binrest < ( 1 / 1e6 ) )
 {
  cout << "!!! Threshold for binding residues must be >0, setting to 0.2 !!!" << endl << endl;
  
  cut_binrest = 0.2;
 }
 else if ( cut_binrest > 1 )
 {
  cout << "!!! Threshold for binding residues must be <=1, setting to 0.2 !!!" << endl << endl;
  
  cut_binrest = 0.2;
 }
 
 if ( cut_clustlig < ( 1 / 1e6 ) )
 {
  cout << "!!! Fingerprint clustering cutoff must be >0, setting to 0.7 !!!" << endl << endl;
  
  cut_clustlig = 0.7;
 }
 else if ( cut_clustlig > 1 )
 {
  cout << "!!! Fingerprint clustering cutoff must be <=1, setting to 0.7 !!!" << endl << endl;
  
  cut_clustlig = 0.7;
 }
 
 if ( met_clustlig != "T" && met_clustlig != "A" )
 {
  cout << "!!! Fingerprint clustering method must be either T or A, setting to T !!!" << endl << endl;
  
  met_clustlig = "T";
 }
 
 if ( met_clustdis != "D" && met_clustdis != "L" )
 {
  cout << "!!! Pocket clustering method must be either D or L, setting to D !!!" << endl << endl;
  
  met_clustdis = "D";
 }
 
 if ( cut_templates > (int) MAXTPL )
 {
  cout << "!!! Max number of templates exceeded, setting to " << MAXTPL << " !!!" << endl << endl;
  
  cut_templates = MAXTPL;
 }
 
 if ( cut_templates < 1 )
 {
  cout << "!!! Max number of templates must be >0, setting to " << MAXTPL << " !!!" << endl << endl;
  
  cut_templates = MAXTPL;
 }
 
 if ( met_druggabl != "R" && met_druggabl != "D" )
 {
  cout << "!!! Druggability model must be either R or D, setting to R !!!" << endl << endl;
  
  met_druggabl = "R";
 }
 
 if ( drug_opt )
 {
  if ( cut_druggabl < ( 1 / 1e6 ) )
  {
   double cut_t = 0.0;
   
   if ( met_druggabl == "R" )
    cut_t = 0.6;
   else
    cut_t = 0.6;
   
   cout << "!!! Threshold for druggability must be >0, setting to " << setprecision(2) << cut_t << " !!!" << endl << endl;
   
   cut_druggabl = cut_t;
  }
  else if ( cut_druggabl > 1 )
  {
   double cut_t = 0.0;
   
   if ( met_druggabl == "R" )
    cut_t = 0.6;
   else
    cut_t = 0.6;
   
   cout << "!!! Threshold for druggability must be <=1, setting to " << setprecision(2) << cut_t << " !!!" << endl << endl;
   
   cut_druggabl = cut_t;
  }
 }
 else
 {
  if ( met_druggabl == "R" )
   cut_druggabl = 0.6;
  else
   cut_druggabl = 0.6;
 }
 
 /* target protein */
 
 Target * target;
 
 target = new Target( 0, 0 );
 
 if ( target->loadTarget(target_name) )
 {
  cout << "Cannot read target structure" << endl;
  exit(EXIT_FAILURE);
 }
 
 if ( target->loadSequence(sequence_name) )
 {
  cout << "Cannot read sequence profile" << endl;
  exit(EXIT_FAILURE);
 }
 
 /* auxiliary compounds */
 
 Cmps * compounds;
 
 compounds = new Cmps( 0 );
 
 if ( cmps_opt )
 {
  if ( compounds->loadCompounds(cmps_name) )
  {
   cout << "Cannot read auxiliary compounds" << endl;
   exit(EXIT_FAILURE);
  }
 }
 
 /* template list */
 
 list<string> template_list;
 
 map<string,double> template_prob1;
 map<string,double> template_prob2;
 
 getList( templates_name, template_list, template_prob1, template_prob2 );
 
 cout << "eFindSite library: " << lib_path <<  endl << endl
 
      << "Number of ligand-bound templates: " << setw(5) << template_list.size() << endl << endl;
 
 /* sequence identity */
 
 cout << "Template filtering ... " << flush;
 
 time(&t_bench1);
 
 multimap< int, Template *, greater<int> > template_set;
 
 list<string>::iterator it1;
 
 for ( it1 = template_list.begin(); it1 != template_list.end(); it1++ )
  if ( template_set.size() < cut_templates )
   if ( template_prob1.find(*it1)->second >= cut_probab )
   {
    Template * template_tmp = new Template( 0, 0, 0, template_prob1.find(*it1)->second, template_prob2.find(*it1)->second );
    
    bool load1 = template_tmp->loadTemplate( lib_path+"/data/"+(*it1).substr(1,2)+"/"+(*it1)+".gz" );
    
    if ( !load1 )
    {
     double sid1 = template_tmp->alignNW( target->getProteinSequence() );
     
     if ( sid1 <= cut_seqid )
      template_set.insert( std::pair< int,Template * >( template_tmp->getProteinResiduesTotal(), template_tmp ) );
    }
   }
 
 time(&t_bench2);
 
 cout << template_set.size() << " templates survived (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 time(&t_bench1);
 
 /* structure alignments */
 
 cout << "Calculating alignments ... " << flush;
 
 int  t_len1;
 
 char t_seq1[MAXPRO];
 
 int t_res1[MAXPRO];
 
 double t_xyz1[MAXPRO][3];
 
 t_len1 = target->getProteinResiduesTotal();
 
 char t_seqs[MAXPRO];
 
 strcpy(t_seqs, (target->getProteinSequence()).c_str());
 
 for ( int t_i = 0; t_i < target->getProteinResiduesTotal(); t_i++ )
  t_seq1[t_i] = t_seqs[t_i];
 
 for ( int t_i = 0; t_i < t_len1; t_i++ )
  t_res1[t_i] = t_i + 1;
 
 target->getProteinCoordsCA(t_xyz1);
 
 std::multimap< int, Template *, greater<int> >::iterator tpl1;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); tpl1++ )
 {
  int  t_len2;
  
  char t_seq2[MAXPRO];
  
  int t_res2[MAXPRO];
  
  double t_xyz2[MAXPRO][3];
  
  t_len2 = ((*tpl1).second)->getProteinResiduesTotal();
  
  strcpy(t_seq2, (((*tpl1).second)->getProteinSequence()).c_str());
  
  for ( int t_i = 0; t_i < t_len2; t_i++ )
   t_res2[t_i] = t_i + 1;
  
  ((*tpl1).second)->getProteinCoordsCA(t_xyz2);
  
  int    t_sco1;
  double t_sco2, t_sco3, t_sco4;
  
  int    t_alig[MAXPRO];
  double t_t[3];
  double t_u[3][3];
  
  frtmalign_( &t_sco1, &t_sco2, &t_sco3, &t_sco4, &t_len2, &t_len1, &t_seq2, &t_seq1, &t_alig, &t_res2, &t_res1, &t_xyz2, &t_xyz1, &t_t, &t_u, &t_len1 );
  
  t_sco4 = 0;
  
  for ( int t_i = 0; t_i < t_len1; t_i++ )
   if ( t_alig[t_i] > -1 )
    if ( t_seq1[t_i] == t_seq2[t_alig[t_i]] )
     t_sco4++;
  
  t_sco4 = t_sco4 / ( (double) t_sco1 );
  
  ((*tpl1).second)->setProteinLengthTM( t_sco1 );
  ((*tpl1).second)->setProteinRMSD( t_sco2 );
  ((*tpl1).second)->setProteinTMscore( t_sco3 );
  ((*tpl1).second)->setProteinSeqID2( t_sco4 );
  ((*tpl1).second)->setTMalignment(t_alig, t_len1);
  ((*tpl1).second)->setMatrix(t_t, t_u);
 }
 
 list<string> template_list_filtered;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); )
 {
  std::multimap< int, Template *, greater<int> >::iterator tpl6 = tpl1++;
  
  if ( ((*tpl6).second)->getProteinTMscore() < cut_tmscore )
   template_set.erase(tpl6);
  else
   template_list_filtered.push_back( ((*tpl6).second)->getProteinID() );
 }
 
 int ltot = 0;
 
 for ( tpl1 = template_set.begin(); tpl1 != template_set.end(); tpl1++ )
  ltot += ((*tpl1).second)->getLigandsTotal();
 
 if ( template_set.empty() )
 {
  cout << "no templates survived" << endl << endl;
  
  time(&t_end);
  
  printTime( difftime(t_end, t_start) );
  
  exit(EXIT_SUCCESS);
 }
 
 time(&t_bench2);
 
 cout << template_set.size() << "/" << ltot << " templates/ligands survived (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Detecting pockets ... " << flush;
 
 time(&t_bench1);
 
 int * clu1 = new int [ltot];
 
 int clu2 = 0;
 
 if ( met_clustdis == "D" && template_set.size() > 4 )
 {
  double *db_xyz = new double[ltot*3];
  
  int lg_xyz = 0;
  
  std::multimap< int, Template *, greater<int> >::iterator tpl2;
  
  for ( tpl2 = template_set.begin(); tpl2 != template_set.end(); tpl2++ )
   for ( int il1 = 0; il1 < ((*tpl2).second)->getLigandsTotal(); il1++ )
   {
    double t_xyz[3];
    
    ((*tpl2).second)->getLigandCenter(il1, t_xyz, true);
    
    for ( int il2 = 0; il2 < 3; il2++ )
     db_xyz[lg_xyz++] = t_xyz[il2];
   }
  
  clu2 = cluster_dbscan( db_xyz, clu1, ltot, cut_clustdis, 2, 1 );
  
  delete [] db_xyz;
 }
 
 else
 {
  double *sim1 = new double[ltot*ltot];
  
  std::multimap< int, Template *, greater<int> >::iterator tpl2;
  std::multimap< int, Template *, greater<int> >::iterator tpl3;
  
  int nl1 = 0;
  int nl2 = 0;
  
  for ( tpl2 = template_set.begin(); tpl2 != template_set.end(); tpl2++ )
   for ( int il1 = 0; il1 < ((*tpl2).second)->getLigandsTotal(); il1++ )
   {
    for ( tpl3 = template_set.begin(); tpl3 != template_set.end(); tpl3++ )
     for ( int il2 = 0; il2 < ((*tpl3).second)->getLigandsTotal(); il2++ )
     {
      double tcen1[3];
      double tcen2[3];
      
      ((*tpl2).second)->getLigandCenter(il1, tcen1, true);
      ((*tpl3).second)->getLigandCenter(il2, tcen2, true);
      
      sim1[nl1*ltot+nl2] = sqrt(pow( tcen1[0] - tcen2[0], 2.0) + pow( tcen1[1] - tcen2[1], 2.0) + pow( tcen1[2] - tcen2[2], 2.0));
      
      nl2++;
     }
    
    nl2 = 0;
    
    nl1++;
   }
  
  clu2 = cluster_avelink( sim1, clu1, nl1, cut_clustdis, "min" );
  
  delete [] sim1;
 }
 
 if ( clu2 < 1 )
 {
  cout << "no pockets found" << endl << endl;
  
  template_set.clear();
  
  time(&t_end);
  
  printTime( difftime(t_end, t_start) );
  
  exit(EXIT_SUCCESS);
 }
 
 time(&t_bench2);
 
 cout << clu2 << " pockets found (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 int nl3 = 0;
 
 std::multimap< int, Template *, greater<int> >::iterator tpl4;
 
 for ( tpl4 = template_set.begin(); tpl4 != template_set.end(); tpl4++ )
  for ( int il1 = 0; il1 < ((*tpl4).second)->getLigandsTotal(); il1++ )
   ((*tpl4).second)->setPocketNumber(il1, clu1[nl3++]);
 
 list< Pocket * > pocket_set;
 
 for ( int clu3 = 0; clu3 < clu2; clu3++ )
 {
  Pocket * pocket_tmp = new Pocket( clu3 );
  
  std::multimap< int, Template *, greater<int> >::iterator tpl5;
  
  for ( tpl5 = template_set.begin(); tpl5 != template_set.end(); tpl5++ )
   for ( int il1 = 0; il1 < ((*tpl5).second)->getLigandsTotal(); il1++ )
    if ( ((*tpl5).second)->getPocketNumber(il1) == clu3 )
     pocket_tmp->addTemplate((*tpl5).second);
  
  if ( pocket_tmp->getProteinsTotal() > 0 && pocket_tmp->getLigandsTotal() > 0 )
   pocket_set.push_back( pocket_tmp );
  else
   delete pocket_tmp;
 }
 
 delete [] clu1;
 
 cout << "Loading SVM models " << flush;
 
 ModelSVM * model_svm;
 
 model_svm = new ModelSVM( false, false, false, false, false, false, false );
 
 time(&t_bench1);
 
 model_svm->loadModel( 1, model_path+"/residues.model" ); cout << '.' << flush;
 model_svm->loadModel( 2, model_path+"/ligands.model" ); cout << '.' << flush;
 model_svm->loadModel( 3, model_path+"/conf1.model" ); cout << '.' << flush;
 model_svm->loadModel( 4, model_path+"/conf2.model" ); cout << '.' << flush;
 
 model_svm->loadScale( 1, model_path+"/residues.scale" ); cout << '.' << flush;
 model_svm->loadScale( 2, model_path+"/ligands.scale" ); cout << '.' << flush;
 model_svm->loadScale( 3, model_path+"/conf1.scale" ); cout << '.' << flush;
 model_svm->loadScale( 4, model_path+"/conf2.scale" ); cout << '.' << flush;
 
 time(&t_bench2);
 
 cout << " done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Predicting binding residues ... " << flush;
 
 time(&t_bench1);
 
 list< Pocket * > pocket_set_filtered;
 
 list< Pocket * >::iterator ipkt1;
 
 for ( ipkt1 = pocket_set.begin(); ipkt1 != pocket_set.end(); ipkt1++ )
 {
  (*ipkt1)->calculatePocketCenter();
  
  int bres1 = (*ipkt1)->calculateBindingResidues( target, model_svm, cut_binrest );
  
  if ( bres1 >= cut_binresn )
   pocket_set_filtered.push_back( *ipkt1 );
 }
 
 pocket_set.clear();
 
 time(&t_bench2);
 
 if ( pocket_set_filtered.empty() )
 {
  cout << "no pockets survived" << endl << endl;
  
  template_set.clear();
  
  pocket_set_filtered.clear();
  
  time(&t_end);
  
  printTime( difftime(t_end, t_start) );
  
  exit(EXIT_SUCCESS);
 }
 
 cout << pocket_set_filtered.size() << " pockets survived (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 double fra1 = 0.0;
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
  fra1 += (double) (*ipkt1)->getLigandsTotal();
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
 {
  double fra2 = 0.0;
  
  if ( fra1 > 0.0 )
   fra2 = ( (double) (*ipkt1)->getLigandsTotal() ) / fra1;
  
  (*ipkt1)->setPocketFraction( fra2 );
 }
 
 cout << "Constructing ligand fingerprints ... " << flush;
 
 time(&t_bench1);
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
 {
  (*ipkt1)->calculateFingerprintsSMILES( cut_clustlig, met_clustlig );
  (*ipkt1)->calculateFingerprintsMACCS( cut_clustlig, met_clustlig );
  
  if ( cmps_opt )
   (*ipkt1)->calculateCmpsScores( compounds, model_svm );
 }
 
 time(&t_bench2);
 
 cout << "done (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 cout << "Ranking pockets ... " << flush;
 
 time(&t_bench1);
 
 double rank1 = 0.0;
 std::string rank3;
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
 {
  double rank2 = (*ipkt1)->calculateConfidence( cmps_opt, model_svm );
  
  (*ipkt1)->calculateDruggability( met_druggabl, cut_binrest, cut_druggabl );
  
  if ( rank2 > rank1 )
  {
   rank1 = rank2;
   
   if ( (*ipkt1)->getDruggable() )
    rank3 = "druggable";
   else
    rank3 = "non-druggable";
  }
 }
 
 time(&t_bench2);
 
 cout << "top-ranked pocket has a confidence index of " << fixed << setprecision(1) << rank1 * 100 << "% and is " << rank3 << " (" << fixed << setprecision(0) << difftime(t_bench2, t_bench1) << " s)" << endl << endl;
 
 multimap<double,Pocket *> pocket_map_sorted;
 
 for ( ipkt1 = pocket_set_filtered.begin(); ipkt1 != pocket_set_filtered.end(); ipkt1++ )
  pocket_map_sorted.insert( pair<double,Pocket *>(-1.0*(*ipkt1)->getConfidence(),*ipkt1) );
 
 list< Pocket * > pocket_set_sorted;
 
 multimap<double,Pocket *>::iterator ipkt3;
 
 for ( ipkt3 = pocket_map_sorted.begin() ; ipkt3 != pocket_map_sorted.end(); ipkt3++ )
  pocket_set_sorted.push_back( (*ipkt3).second );
 
 pocket_set_filtered.clear();
 
 map<string,bool> chk1;
 map<string,bool> chk2;
 
 int ipkt2 = 1;
 
 for ( ipkt1 = pocket_set_sorted.begin(); ipkt1 != pocket_set_sorted.end(); ipkt1++ )
 {
  (*ipkt1)->setCenter( cut_binrest, cut_clustdis );
  
  (*ipkt1)->dumpProteinAlignments( output_name, chk1, target );
  
  (*ipkt1)->dumpPocket( output_name, target, cut_binrest, ipkt2 );
  
  (*ipkt1)->dumpLigands( output_name, chk2, ipkt2 );
  
  ipkt2++;
 }
 
 template_set.clear();
 
 pocket_set_sorted.clear();
 
 time(&t_end);
 
 printTime( difftime(t_end, t_start) );
 
 return 0;
}
