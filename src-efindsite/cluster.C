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


#include "cluster.h"

using namespace std;


// ==================================================================================   average linkage clustering

int cluster_avelink( double csim1[], int csim2[], int csim3, double csim4, std::string csim5 )
{
 list< list<int> > clusters;
 
 for ( int ci1 = 0; ci1 < csim3; ci1++ )
 {
  list<int> cluster (1, ci1);
  
  clusters.push_back( cluster );
  
  csim2[ci1] = 0;
 }
 
 int cnum1 = clusters.size();
 
 bool cw1 = true;
 
 while ( cw1 && cnum1 > 1 )
 {
  list< list<int> >::iterator ci2;
  list< list<int> >::iterator ci3;
  
  int ci4, ci5;
  
  double cmin1 = 1000000.0;
  
  if ( csim5 == "max" )
   cmin1 = -1000000.0;
  
  list< list<int> >::iterator cmin2;
  list< list<int> >::iterator cmin3;
  
  for ( ci2 = clusters.begin(), ci4 = 0; ci2 != clusters.end(); ci2++, ci4++ )
  {
   for ( ci3 = clusters.begin(), ci5 = 0; ci3 != clusters.end(); ci3++, ci5++ )
   {
    if ( ci4 < cnum1 - 1 && ci5 > ci4 )
    {
     double amin1 = 0.0;
     double amin2 = 0.0;
     
     list<int>::iterator ci6;
     list<int>::iterator ci7;
     
     if ( !(*ci2).empty() && !(*ci3).empty() )
     {
      for ( ci6 = (*ci2).begin(); ci6 != (*ci2).end(); ci6++ )
       for ( ci7 = (*ci3).begin(); ci7 != (*ci3).end(); ci7++ )
       {
        if ( *ci6 < *ci7 )
         amin1 += csim1[ (*ci6) * csim3 + (*ci7) ];
        else
         amin1 += csim1[ (*ci7) * csim3 + (*ci6) ];
        
        amin2 += 1.0;
       }
       
      if ( amin2 > 0.1 )
       amin1 /= amin2;
      
      if ( ( csim5 == "min" && amin1 < cmin1 ) || ( csim5 == "max" && amin1 > cmin1 ) )
      {
       cmin1 = amin1;
       cmin2 = ci2;
       cmin3 = ci3;
      }
     }
    }
   }
  }
  
  if ( ( csim5 == "min" && cmin1 <= csim4 ) || ( csim5 == "max" && cmin1 >= csim4 ) )
  {
   (*cmin2).merge(*cmin3);
   
   (*cmin3).clear();
   
   cnum1--;
  }
  else
  {
   cw1 = false;
  }
 }
 
 multimap<int, list<int> > clusters_sorted;
 
 list< list<int> >::iterator ci8;
 
 for ( ci8 = clusters.begin(); ci8 != clusters.end(); ci8++ )
  if ( !(*ci8).empty() )
   clusters_sorted.insert( pair<int, list<int> >((*ci8).size(),*ci8) );
 
 multimap<int, list<int> >::iterator ci9;
 
 int cnum2 = cnum1;
 
 for ( ci9=clusters_sorted.begin() ; ci9 != clusters_sorted.end(); ci9++ )
 {
  cnum2--;
  
  for ( list<int>::iterator ci10 = ((*ci9).second).begin(); ci10 != ((*ci9).second).end(); ci10++)
   csim2[*ci10] = cnum2;
 }
 
 return cnum1;
}


// ==================================================================================   DBSCAN

int cluster_dbscan( double l_xyz[], int c_def[], int l_tot, double eps, int minPts, int threads )
{
 omp_set_num_threads(threads);
 
 NWUClustering::ClusteringAlgo dbs;
 
 dbs.set_dbscan_params(eps, minPts);
 
 dbs.input_data(l_tot, 3, l_xyz);
 
 dbs.build_kdtree();
 
 run_dbscan_algo_uf(dbs);
 
 int cnum1 = dbs.getClusters_uf(c_def);
 
 for ( int i = 0; i < l_tot; i++ )
  if ( c_def[i] < 0 )
   c_def[i] = cnum1++;
 
 return cnum1;
}
