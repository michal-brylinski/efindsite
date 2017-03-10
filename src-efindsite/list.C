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


#include "list.h"

using namespace std;

void getList( string templates_name, list<string> &template_list, map<string,double> &template_prob1, map<string,double> &template_prob2 )
{
 string line1;
 
 ifstream l1_file( templates_name.c_str() );
 
 if ( !l1_file.is_open() ) { cout << "Cannot open " << templates_name << endl; exit(EXIT_FAILURE); }
 
 while (getline(l1_file,line1))
 {
  istringstream iss1(line1);
  
  int iss2 = -1;
  
  std::string template_id;
  
  double prob1[2];
  
  prob1[0] = 1.0;
  prob1[1] = 1.0;
  
  while ( !iss1.eof() )
  {
   string iss3;
   
   getline( iss1, iss3, ' ' );
   
   if ( iss2 < 0 )
    template_id = iss3;
   else
    prob1[iss2] = atof(iss3.c_str());
   
   iss2++;
  }
  
  if ( prob1[0] < 0.001 )
   prob1[0] = 0.001;
  
  if ( prob1[0] > 1.0 )
   prob1[0] = 1.0;
  
  if ( prob1[1] < 0.001 )
   prob1[1] = 0.001;
  
  if ( prob1[1] > 1.0 )
   prob1[1] = 1.0;
  
  template_prob1.insert( pair<string,double>(template_id, prob1[0]) );
  template_prob2.insert( pair<string,double>(template_id, prob1[1]) );
  
  template_list.push_back( template_id );
 }
 
 l1_file.close();
}
