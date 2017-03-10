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


#ifndef __RUNSVM_H_
#define __RUNSVM_H_

#include<cstdlib>
#include<fstream>
#include<iostream>
#include<sstream>
#include<cstdio>
#include<cctype>
#include<string>

#include "size.h"
#include "svm.h"

using namespace std;

class ModelSVM {
        
  private:
    
    bool _bres_loaded; /* binding residues flag */
    bool _alig_loaded; /* auxiliary ligands flag */
    bool _cnf1_loaded; /* confidence flag */
    bool _cnf2_loaded; /* confidence ligands flag */
    bool _vsco_loaded; /* svm scoring for vs flag */
    bool _scr1_loaded; /* screening conf 1% flag */
    bool _scr2_loaded; /* screening conf 10% flag */
    
    int _efs_attr;
    
    struct svm_model * _efs_model_bres;
    struct svm_node * _efs_node_bres;
    double _efs_scale_bres[MAXSV1][2];
    
    struct svm_model * _efs_model_alig;
    struct svm_node * _efs_node_alig;
    double _efs_scale_alig[MAXSV2][2];
    
    struct svm_model * _efs_model_cnf1;
    struct svm_node * _efs_node_cnf1;
    double _efs_scale_cnf1[MAXSV3][2];
    
    struct svm_model * _efs_model_cnf2;
    struct svm_node * _efs_node_cnf2;
    double _efs_scale_cnf2[MAXSV4][2];
    
    struct svm_model * _efs_model_vsco;
    struct svm_node * _efs_node_vsco;
    double _efs_scale_vsco[MAXSV5][2];
    
    struct svm_model * _efs_model_scr1;
    struct svm_node * _efs_node_scr1;
    double _efs_scale_scr1[MAXSV6][2];
    
    struct svm_model * _efs_model_scr2;
    struct svm_node * _efs_node_scr2;
    double _efs_scale_scr2[MAXSV6][2];
    
  public:
    
    ModelSVM( bool, bool, bool, bool, bool, bool, bool );
    
    ModelSVM( void );
    
    ~ModelSVM();
    
    void loadModel( int, std::string );
    
    void loadScale( int, std::string );
    
    double SVMpredict( int, double [] );
};

#endif
