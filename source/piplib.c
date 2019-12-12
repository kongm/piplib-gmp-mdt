/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                 piplib.c                                   *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988-2005                                        *
 *                                                                            *
 * This is free software; you can redistribute it and/or modify it under the  *
 * terms of the GNU General Public License as published by the Free Software  *
 * Foundation; either version 2 of the License, or (at your option) any later *
 * version.                                                                   *
 *                                                                            *
 * This software is distributed in the hope that it will be useful, but       *
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY *
 * or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License   *
 * for more details.                                                          *
 *                                                                            *
 * You should have received a copy of the GNU General Public License along    *
 * with software; if not, write to the Free Software Foundation, Inc.,        *
 * 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA                     *
 *                                                                            *
 * Written by Cedric Bastoul                                                  *
 *                                                                            *
 ******************************************************************************/

/* Premiere version du 30 juillet 2001. */

# include <stdlib.h>
# include <string.h>
# include <stdio.h>
# include <ctype.h>

#include <piplib-gmp/piplib-gmp.h>
#define min(x,y) ((x) < (y)? (x) : (y))


#ifndef PIPMP_LINEAR_VALUE_IS_MP
# define PIPMP_LINEAR_VALUE_IS_MP
#endif


EntierMP pipmp_UN;
EntierMP pipmp_ZERO;

long int pipmp_cross_product, pipmp_limit;
int pipmp_allocation, pipmp_comptage;
int pipmp_verbose = 0;
int pipmp_profondeur = 0;
int pipmp_compa_count;
int pipmp_deepest_cut = 0;

FILE *pipmp_dump = NULL;
char pipmp_dump_name[] = "PipXXXXXX";

/* Larger line buffer to accomodate Fr�do Vivien exemples. A version
handling arbitrary line length should be written ASAP.
*/

#define INLENGTH 2048

char pipmp_inbuff[INLENGTH];
int pipmp_inptr = 256;
int pipmp_proviso = 0;


/******************************************************************************
 *                 Fonctions d'acquisition de donn�es (ex-maind.c)            *
 ******************************************************************************/


int pipmp_dgetc(FILE *foo)
{
 char *p;
 if(pipmp_inptr >= pipmp_proviso)
   {p = fgets(pipmp_inbuff, INLENGTH, foo);
    if(p == NULL) return EOF;
    pipmp_proviso = min(INLENGTH, strlen(pipmp_inbuff));
    pipmp_inptr = 0;
    if(pipmp_verbose > 2) fprintf(pipmp_dump, "-- %s", pipmp_inbuff);
  }
 return pipmp_inbuff[pipmp_inptr++];
}



int pipmp_dscanf(FILE *foo, EntierMP  val)
{
 char * p;
 int c;

 for(;pipmp_inptr < pipmp_proviso; pipmp_inptr++)
   if(pipmp_inbuff[pipmp_inptr] != ' ' && pipmp_inbuff[pipmp_inptr] != '\n' && pipmp_inbuff[pipmp_inptr] != '\t')
				break;
 while(pipmp_inptr >= pipmp_proviso)
   {p = fgets(pipmp_inbuff, 256, foo);
    if(p == NULL) return EOF;
    pipmp_proviso = strlen(pipmp_inbuff);
    if(pipmp_verbose > 2) {
      fprintf(pipmp_dump, ".. %s", pipmp_inbuff);
      fflush(pipmp_dump);
    }
    for(pipmp_inptr = 0; pipmp_inptr < pipmp_proviso; pipmp_inptr++)
       if(pipmp_inbuff[pipmp_inptr] != ' '
       && pipmp_inbuff[pipmp_inptr] != '\n'
       && pipmp_inbuff[pipmp_inptr] != '\t') break;
  }
 if(gmp_sscanf(pipmp_inbuff+pipmp_inptr, GMP_INPUT_FORMAT, val) != 1)
 return -1;
 
 for(; pipmp_inptr < pipmp_proviso; pipmp_inptr++)
	if((c = pipmp_inbuff[pipmp_inptr]) != '-' && !isdigit(c)) break;
 return 0;
}


/******************************************************************************
 *                    Fonctions d'affichage des structures                    *
 ******************************************************************************/


/* Fonction pipmp_matrix_print :
 * Cette fonction se charge d'imprimer sur le flux 'foo' les informations
 * que contient la structure de type PipMPMatrix qu'elle recoit en parametre.
 * Premiere version : Ced. 29 juillet 2001. 
 */
void pipmp_matrix_print(FILE * foo, PipMPMatrix * Mat)
{ EntierMP * p;
  int i, j ;
  unsigned NbRows, NbColumns ;

  fprintf(foo,"%d %d\n", NbRows=Mat->NbRows, NbColumns=Mat->NbColumns) ;
  for (i=0;i<NbRows;i++) 
  { p=*(Mat->p+i) ;
    for (j=0;j<NbColumns;j++)
    #if defined(PIPMP_LINEAR_VALUE_IS_MP)
    { fprintf(foo," ") ;
      mpz_out_str(foo,10,*p++) ;
    }
    #else
    fprintf(foo," %3d", *p++) ;
    #endif
    fprintf(foo, "\n") ;
  }
} 


/* Fonction pipmp_vector_print :
 * Cette fonction se charge d'imprimer sur le flux 'foo' les informations
 * que contient la structure de type PipMPVector qu'elle recoit en parametre.
 * Premiere version : Ced. 20 juillet 2001. 
 */
void pipmp_vector_print(FILE * foo, PipMPVector * vector)
{ int i ;
  
  if (vector != NULL)
  { fprintf(foo,"#[") ;
    for (i=0;i<vector->nb_elements;i++)
    { fprintf(foo," ") ;
      #if defined(PIPMP_LINEAR_VALUE_IS_MP)
      mpz_out_str(foo,10,vector->the_vector[i]) ;
      if (mpz_cmp(vector->the_deno[i],pipmp_UN) != 0)
      #else
      fprintf(foo,FORMAT,vector->the_vector[i]) ;
      if (vector->the_deno[i] != pipmp_UN)
      #endif
      { fprintf(foo,"/") ;
        #if defined(PIPMP_LINEAR_VALUE_IS_MP)
        mpz_out_str(foo,10,vector->the_deno[i]) ;
        #else
        fprintf(foo,FORMAT,vector->the_deno[i]) ;
        #endif
      }
    }
    fprintf(foo,"]") ;
  }
}
  

/* Fonction pipmp_newparm_print :
 * Cette fonction se charge d'imprimer sur le flux 'foo' les informations
 * que contient la structure de type PipMPNewparm qu'elle recoit en parametre.
 * Le parametre indent est le nombre d'espaces blancs en debut de chaque
 * ligne avant indentation. Une valeur negative de indent signifie qu'on ne
 * desire pas d'indentation.
 * Premiere version : Ced. 18 octobre 2001. 
 */
void pipmp_newparm_print(FILE * foo, PipMPNewparm * newparm, int indent)
{ int i ;

  if (newparm != NULL)
  { do
    { for (i=0;i<indent;i++) fprintf(foo," ") ;             /* Indent. */
      fprintf(foo,"(newparm ") ;
      fprintf(foo,"%d",newparm->rank) ;
      fprintf(foo," (div ") ;
      pipmp_vector_print(foo,newparm->vector) ;
      fprintf(foo," ") ;
      #if defined(PIPMP_LINEAR_VALUE_IS_MP)
      mpz_out_str(foo,10,newparm->deno) ;
      #else
      fprintf(foo,FORMAT,newparm->deno) ;
      #endif
      fprintf(foo,"))\n") ;
    }
    while ((newparm = newparm->next) != NULL) ;
  }  
}


/* Fonction pipmp_list_print :
 * Cette fonction se charge d'imprimer sur le flux 'foo' les informations
 * que contient la structure de type PipMPList qu'elle recoit en parametre.
 * Le parametre indent est le nombre d'espaces blancs en debut de chaque
 * ligne avant indentation. Une valeur negative de indent signifie qu'on ne
 * desire pas d'indentation.
 * Premiere version : Ced. 18 octobre 2001. 
 * 16 novembre 2005 : Ced. Prise en compte du cas list->vector == NULL,
 *                         jusque l� impossible.
 */
void pipmp_list_print(FILE * foo, PipMPList * list, int indent)
{ int i ;

  if (list == NULL)
  { for (i=0;i<indent;i++) fprintf(foo," ") ;               /* Indent. */
    fprintf(foo,"()\n") ;
  }
  else
  { for (i=0;i<indent;i++) fprintf(foo," ") ;               /* Indent. */
    fprintf(foo,"(list\n") ;
    do
    { if (list->vector != NULL)
      { for (i=0;i<indent+1;i++) fprintf(foo," ") ;         /* Indent. */
        pipmp_vector_print(foo,list->vector) ;
        fprintf(foo,"\n") ;
      }
    }
    while ((list = list->next) != NULL) ;
    for (i=0;i<indent;i++) fprintf(foo," ") ;               /* Indent. */
    fprintf(foo,")\n") ;
  }
}


/* Fonction pipmp_quast_print :
 * Cette fonction se charge d'imprimer sur le flux 'foo' les informations
 * que contient la structure de type PipMPQuast qu'elle recoit en parametre.
 * Le parametre indent est le nombre d'espaces blancs en debut de chaque
 * ligne avant indentation. Une valeur negative de indent signifie qu'on ne
 * desire pas d'indentation.
 * 20 juillet 2001 : Premiere version, Ced. 
 * 18 octobre 2001 : eclatement. 
 */
void pipmp_quast_print(FILE * foo, PipMPQuast * solution, int indent)
{ int i ;
  PipMPVector * vector ;
  
  if (solution != NULL)
  { pipmp_newparm_print(foo,solution->newparm,indent) ;
    if (solution->condition == NULL)
    pipmp_list_print(foo,solution->list,indent) ;
    else
    { for (i=0;i<indent;i++) fprintf(foo," ") ;             /* Indent. */
      fprintf(foo,"(if ") ;
      pipmp_vector_print(foo,solution->condition) ;
      fprintf(foo,"\n") ;
      if (indent>=0)                                        /* Indent. */
      { pipmp_quast_print(foo,solution->next_then,indent+1) ;
        pipmp_quast_print(foo,solution->next_else,indent+1) ;
      }
      else
      { pipmp_quast_print(foo,solution->next_then,indent) ;
        pipmp_quast_print(foo,solution->next_else,indent) ;
      }
      for (i=0;i<indent;i++) fprintf(foo," ") ;             /* Indent. */
      fprintf(foo,")\n") ;
    }
  }
  else
  { for (i=0;i<indent;i++) fprintf(foo," ") ;               /* Indent. */
    fprintf(foo,"void\n") ;
  }
}    
  

/* Function pipmp_options_print:
 * This function prints the content of a PipMPOptions structure (options)
 * into a file (foo, possibly stdout).
 * March 17th 2003: first version.
 */
void * pipmp_options_print(FILE * foo, PipMPOptions * options)
{ fprintf(foo,"Option setting is:\n") ;
  fprintf(foo,"Nq          =%d\n",options->Nq) ;
  fprintf(foo,"Verbose     =%d\n",options->Verbose) ;
  fprintf(foo,"Simplify    =%d\n",options->Simplify) ;
  fprintf(foo,"Deepest_cut =%d\n",options->Deepest_cut) ;
  fprintf(foo,"Maximize    =%d\n",options->Maximize) ;
  fprintf(foo,"Urs_parms   =%d\n",options->Urs_parms);
  fprintf(foo,"Urs_unknowns=%d\n",options->Urs_unknowns);
  fprintf(foo,"\n") ;
}


/******************************************************************************
 *                       Fonctions de liberation memoire                      *
 ******************************************************************************/


/* Fonction pipmp_matrix_free :
 * Cette fonction libere la memoire reservee a la structure de type PipMPMatrix
 * que pointe son parametre.
 * Premiere version : Ced. 29 juillet 2001. 
 */
void pipmp_matrix_free(PipMPMatrix * matrix)
{ 
  #if defined(PIPMP_LINEAR_VALUE_IS_MP)
  int i, j ;
  EntierMP * p ;

  p = matrix->p_Init ;
  for (i=0;i<matrix->p_Init_size;i++) 
  mpz_clear(*p++) ;
  #endif

  if (matrix != NULL)
  { free(matrix->p_Init) ;
    free(matrix->p) ;
    free(matrix) ;
  }
}


/* Fonction pipmp_vector_free :
 * Cette fonction libere la memoire reservee a la structure de type PipMPVector
 * que pointe son parametre.
 * 20 juillet 2001 : Premiere version, Ced.
 * 18 octobre 2001 : simplification suite a l'eclatement de PipMPVector.
 * 16 novembre 2005 : Ced. Prise en compte du cas NULL.
 */
void pipmp_vector_free(PipMPVector * vector)
{ int i ;
  
  if (vector != NULL)
  { 
    #if defined(PIPMP_LINEAR_VALUE_IS_MP)
    for (i=0;i<vector->nb_elements;i++)
    { mpz_clear(vector->the_vector[i]);
      mpz_clear(vector->the_deno[i]);
    }
    #endif
  
    free(vector->the_vector) ;
    free(vector->the_deno) ;
    free(vector) ;
  }
}


/* Fonction pipmp_newparm_free :
 * Cette fonction libere la memoire reservee a la structure de type PipMPNewparm
 * que pointe son parametre. Sont liberes aussi tous les elements de la
 * liste chainee dont il pouvait etre le depart.
 * Premiere version : Ced. 18 octobre 2001. 
 */
void pipmp_newparm_free(PipMPNewparm * newparm)
{ PipMPNewparm * next ;

  while (newparm != NULL)
  { next = newparm->next ;
    #if defined(PIPMP_LINEAR_VALUE_IS_MP)
    mpz_clear(newparm->deno);
    #endif
    pipmp_vector_free(newparm->vector) ;
    free(newparm) ;
    newparm = next ;
  }
}


/* Fonction pipmp_list_free :
 * Cette fonction libere la memoire reservee a la structure de type PipMPList
 * que pointe son parametre. Sont liberes aussi tous les elements de la
 * liste chainee dont il pouvait etre le depart.
 * Premiere version : Ced. 18 octobre 2001. 
 */
void pipmp_list_free(PipMPList * list)
{ PipMPList * next ;

  while (list != NULL)
  { next = list->next ;
    pipmp_vector_free(list->vector) ;
    free(list) ;
    list = next ;
  }
}


/* Fonction pipmp_quast_free :
 * Cette fonction libere la memoire reservee a la structure de type 
 * PipMPSolution que pointe son parametre. Sont liberees aussi toutes les
 * differentes listes chainees qui pouvaient en partir.
 * 20 juillet 2001 : Premiere version, Ced.
 * 18 octobre 2001 : simplification suite a l'eclatement de PipMPVector.
 */
void pipmp_quast_free(PipMPQuast * solution)
{ if (solution != NULL)
  { if (solution->newparm != NULL)
    pipmp_newparm_free(solution->newparm) ;
  
    if (solution->list != NULL)
    pipmp_list_free(solution->list) ;

    if (solution->condition != NULL)
    { pipmp_vector_free(solution->condition) ;
      pipmp_quast_free(solution->next_then) ;
      pipmp_quast_free(solution->next_else) ;
    }
    free(solution) ;
  }
}


/* Funtion pipmp_options_free:
 * This function frees the allocated memory for a PipMPOptions structure.
 * March 15th 2003: first version.
 */
void pipmp_options_free(PipMPOptions * options)
{ free(options) ;
}


/******************************************************************************
 *                     Fonction d'initialisation des options                  *
 ******************************************************************************/


/* Funtion pipmp_options_init:
 * This function allocates the memory for a PipMPOptions structure and fill the
 * options with the default values.
 ********
 * Nq est un booleen renseignant si on cherche
 * une solution entiere (vrai=1) ou non (faux=0). Verbose est un booleen
 * permettant de rendre PipMP bavard (Verbose a vrai=1), il imprimera
 * alors la plupart de ses traitements dans le fichier dont le nom est
 * dans la variable d'environnement DEBUG, ou si DEBUG
 * n'est pas placee, dans un nouveau fichier de nom genere par mkstemp, si
 * Verbose est a faux=0, PipMP restera muet. Simplify est un booleen permettant
 * de demander a PipMP de simplifier sa solution (en eliminant les formes de type
 * 'if #[...] () ()') quand il est a vrai=1, ou non quand il est a faux=0. Max
 * n'est pas encore utilise et doit etre mis a 0. 
 ********
 * March 15th 2003: first version.
 */ 
PipMPOptions * pipmp_options_init(void)
{ PipMPOptions * options ;

  options = (PipMPOptions *)malloc(sizeof(PipMPOptions)) ;
  /* Default values of the options. */
  options->Nq          = 1 ;  /* Integer solution. */
  options->Verbose     = 0 ;  /* No comments. */
  options->Simplify    = 0 ;  /* Do not simplify solutions. */
  options->Deepest_cut = 0 ;  /* Do not use deepest cut algorithm. */
  options->Maximize    = 0 ;  /* Do not compute maximum. */
  options->Urs_parms   = 0 ;  /* All parameters are non-negative. */
  options->Urs_unknowns= 0 ;  /* All unknows are non-negative. */
  
  return options ;
}


/******************************************************************************
 *                     Fonctions d'acquisition de matrices                    *
 ******************************************************************************/


/* Fonction pipmp_matrix_alloc :
 * Fonction (tres) modifiee de Matrix_Alloc de la polylib. Elle alloue l'espace
 * memoire necessaire pour recevoir le contenu d'une matrice de NbRows lignes
 * et de NbColumns colonnes, et initialise les valeurs a 0. Elle retourne un
 * pointeur sur l'espace memoire alloue.
 * Premiere version : Ced. 18 octobre 2001. 
 */
PipMPMatrix * pipmp_matrix_alloc(unsigned NbRows, unsigned NbColumns)
{ PipMPMatrix * matrix ;
  EntierMP ** p, * q ;
  int i, j ;

  matrix = (PipMPMatrix *)malloc(sizeof(PipMPMatrix)) ;
  if (matrix == NULL) 	
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  matrix->NbRows = NbRows ;
  matrix->NbColumns = NbColumns ;
  matrix->p_Init_size = NbRows * NbColumns ;
  if (NbRows == 0) 
  { matrix->p = NULL ;
    matrix->p_Init = NULL ;
  }  
  else 
  { if (NbColumns == 0) 
    { matrix->p = NULL ;
      matrix->p_Init = NULL ;
    }
    else 
    { p = (EntierMP **)malloc(NbRows*sizeof(EntierMP *)) ;
      if (p == NULL) 
      { fprintf(stderr, "Memory Overflow.\n") ;
        exit(1) ;
      }
      q = (EntierMP *)malloc(NbRows * NbColumns * sizeof(EntierMP)) ;
      if (q == NULL) 
      { fprintf(stderr, "Memory Overflow.\n") ;
        exit(1) ;
      }
      matrix->p = p ;
      matrix->p_Init = q ;
      for (i=0;i<NbRows;i++) 
      { *p++ = q ;
	for (j=0;j<NbColumns;j++)   
        #if defined(PIPMP_LINEAR_VALUE_IS_MP)
	mpz_init_set_si(*(q+j),0) ;
	#else
	*(q+j) = 0 ;
	#endif
	q += NbColumns ;
      }
    }
  }
  return matrix ;
}


/* Fonction pipmp_matrix_read :
 * Adaptation de Matrix_Read de la polylib. Cette fonction lit les donnees
 * d'une matrice dans un fichier 'foo' et retourne un pointeur vers une
 * structure ou elle a copie les informations de cette matrice. Les donnees
 * doivent se presenter tq :
 * - des lignes de commentaires commencant par # (optionnelles),
 * - le nombre de lignes suivit du nombre de colonnes de la matrice, puis
 *   eventuellement d'un commentaire sur une meme ligne,
 * - des lignes de la matrice, chaque ligne devant etre sur sa propre ligne de
 *   texte et eventuellement suivies d'un commentaire.
 * Premiere version : Ced. 18 octobre 2001. 
 * 24 octobre 2002 : premiere version MP, attention, uniquement capable de
 *                   lire des long long pour l'instant. On utilise pas
 *                   mpz_inp_str car on lit depuis des char * et non des FILE.
 */
PipMPMatrix * pipmp_matrix_read(FILE * foo)
{ unsigned NbRows, NbColumns ;
  int i, j, n ;
  #if defined(PIPMP_LINEAR_VALUE_IS_MP)
  long long val ;
  #endif
  char *c, s[1024], str[1024] ;
  PipMPMatrix * matrix ;
  EntierMP * p ;

  while (fgets(s,1024,foo) == 0) ;
  while ((*s=='#' || *s=='\n') || (sscanf(s," %u %u",&NbRows,&NbColumns)<2))
  fgets(s, 1024, foo) ;
  
  matrix = pipmp_matrix_alloc(NbRows,NbColumns) ;

  p = matrix->p_Init ;
  for (i=0;i<matrix->NbRows;i++) 
  { do 
    { c = fgets(s,1024,foo) ;
      while (isspace(*c) && (*c != '\n'))
      c++ ;
    }
    while (c != NULL && (*c == '#' || *c == '\n'));
    
    if (c == NULL) 
    { fprintf(stderr, "Not enough rows.\n") ;
      exit(1) ;
    }
    for (j=0;j<matrix->NbColumns;j++) 
    { if (c == NULL || *c == '#' || *c == '\n')
      { fprintf(stderr, "Not enough columns.\n") ;
        exit(1) ;
      }
      /* NdCed : Dans le n ca met strlen(str). */
      if (sscanf(c,"%s%n",str,&n) == 0) 
      { fprintf(stderr, "Not enough rows.\n") ;
        exit(1) ;
      }
      #if defined(PIPMP_LINEAR_VALUE_IS_MP)
      sscanf(str,"%lld",&val) ;
      mpz_set_si(*p++,val) ;
      #else
      sscanf(str,FORMAT,p++) ;
      #endif
      c += n ;
    }
  }
  return matrix ;
}
 
 
/* initialization of pip */
static int pipmp_initialized = 0;

void pipmp_init() {
  /* Avoid initializing (and leaking) several times */
  if (!pipmp_initialized) {
    #if defined(PIPMP_LINEAR_VALUE_IS_MP)
    mpz_init_set_si(pipmp_UN, 1);
    mpz_init_set_si(pipmp_ZERO, 0);
    #else
    pipmp_UN   = VAL_UN ;
    pipmp_ZERO = VAL_ZERO ;
    #endif
    solmp_init() ;
    tabmp_init() ;
    pipmp_initialized = 1;
  }
}

void pipmp_close() {
  tabmp_close();
  solmp_close();
# if defined(PIPMP_LINEAR_VALUE_IS_MP)
  mpz_clear(pipmp_UN);
  mpz_clear(pipmp_ZERO);
# endif
  pipmp_initialized = 0;
}
 


/******************************************************************************
 *                           Fonction de resolution                           *
 ******************************************************************************/
TableauMP* pipmp_expanser();

/* Fonction pipmp_solve :
 * Cette fonction fait un appel a PipMP pour la resolution d'un probleme. Le
 * probleme est fourni dans les arguments. Deux matrices de la forme de celles
 * utilisees dans la Polylib forment les systemes d'equations/inequations :
 * un pour les inconnues, l'autre pour les parametres. Bg est le 'bignum'.
 * Le dernier argument contient les options guidant le comportement de PIP.
 * Cette fonction retourne la solution sous la forme d'un arbre de structures
 * PipMPQuast.
 * 30 juillet 2001 : Premiere version, Ced. 
 * 18 octobre 2001 : suppression de l'argument Np, le nombre de parametres. Il
 *                   est a present deduit de ineqpar. Si ineqpar vaut NULL,
 *                   c'est que Np vaut 0. S'il y a des parametres mais pas de
 *                   contraintes dessus, ineqpar sera une matrice de 0 lignes
 *                   mais du bon nombre de colonnes (Np + 2).
 * 27 f�vrier 2003 : Verbose est maintenant gradu�.
 *                  -1 -> silence absolu
 *                   0 -> silence relatif
 *                   1 -> information sur les coupures dans le cas ou on
 *                        cherche une solution enti�re.
 *                   2 -> information sur les pivots et les d�terminants
 *                   3 -> information sur les tableaux.
 *                         Chaque option inclut les pr�c�dentes. [paf]
 * 15 mars 2003    : passage a une structure d'options.
 */
PipMPQuast * pipmp_solve(inequnk, ineqpar, Bg, options)
PipMPMatrix * inequnk, * ineqpar ;
int Bg ;
PipMPOptions * options ;
{
  TableauMP * ineq, * context, * ctxt ;
  int i, Np, Nn, Nl, Nm, p, q, xq, non_vide, Shift = 0, Urs_parms = 0;
  char * g ;
  struct pipmp_high_water_mark hq ;
  EntierMP D ;
  PipMPQuast * solution ;
  int	solmp_flags = 0;
  
  pipmp_init() ;
#if defined(PIPMP_LINEAR_VALUE_IS_MP)
  //  printf ("Hi, I am MP :)\n");
#else
  //printf ("Bouh... I am long long :(\n");
#endif
  /* initialisations diverses :
   * - la valeur de Verbose et Deepest_cut sont placees dans leurs variables
   *   globales. Dans le cas ou on doit etre en mode verbose, on ouvre le
   *   fichier dans lequel ecrire les tracages. Si la variable d'environnement
   *   DEBUG est placee, on ecrira dans le nom de fichier correspondant, sinon,
   *   dans un nouveau fichier de nom genere par mkstemp,
   * - limit est mis au 0 long int (sa valeur par defaut dans PipMP original),
   * - on lance les initialisations pour tab et sol (autres mises en place
   *   de variables globales).
   */
  pipmp_verbose = options->Verbose ;
  pipmp_deepest_cut = options->Deepest_cut ;
  if (pipmp_verbose > 0)
    { g = getenv("DEBUG") ;
    if(g && *g)
      { pipmp_dump = fopen(g, "w") ;
      if(pipmp_dump == NULL)
	{ fprintf(stderr,"%s unaccessible\n",g) ;
	pipmp_verbose = 0 ;
	}
      }
    else
      { mkstemp(pipmp_dump_name) ;
      pipmp_dump = fopen(pipmp_dump_name, "w") ;
      }
    }
#if defined(PIPMP_LINEAR_VALUE_IS_MP)
  pipmp_limit = 0LL ;
#else
  pipmp_limit = pipmp_ZERO ;
#endif

  /* Si inequnk est NULL, la solution est automatiquement void (NULL). */
  if (inequnk != NULL)
    { /* Np vaut 0 si ineqpar vaut NULL, ineqpar->NbColumns - 2 sinon (on a -1
       * pour la constante et -1 pour le marqueur d'egalite/inegalite = -2).
       */
      Np = (ineqpar == NULL) ? 0 : ineqpar->NbColumns - 2 ;
      /* Calcul du nombre d'inconnues du probleme. Comme les matrices d'entree
       * sont sous la forme de la polylib.
       */
      Nn = inequnk->NbColumns - Np - 2 ;
      /* Calcul du nombre d'inequations du probleme. Le format de matrice de la
       * polylib permet les egalites, on doit donc les compter double quand il
       * y en a.
       */
      Nl = inequnk->NbRows ;
      for (i=0;i<inequnk->NbRows;i++)
#if defined(PIPMP_LINEAR_VALUE_IS_MP)
	if (mpz_sgn(**(inequnk->p + i)) == 0)
#else
	  if (**(inequnk->p + i) == 0)
#endif
	    Nl ++ ;
      
      /* On prend les 'marques' de debut de traitement. */
      pipmp_cross_product = 0 ;
      hq = tabmp_hwm() ;
      xq = p = solmp_hwm();
    
      if (options->Maximize) {
	solmp_flags |= SOL_MAX;
	Shift = 1;
      } else if (options->Urs_unknowns) {
	solmp_flags |= SOL_SHIFT;
	Shift = -1;
      } 

      if (options->Urs_parms) {
	Urs_parms = Np - (Bg >= 0);
	Np += Urs_parms;
      }

      /* Si un maximum est demande, mais sans bignum, on cr�e le bignum. */
      if (options->Maximize || options->Urs_unknowns) {
	if (Bg < 0) {
	  Bg = inequnk->NbColumns - 1 ; /* On choisit sa place. */
	  Np ++ ;                       /* On le compte comme parametre. */
	  solmp_flags |= SOL_REMOVE;      /* On le supprime apres. */
	}
      }
    
      /* On s'assure d'abord que le systeme pour le contexte n'est pas vide
       * avant de commencer le traitement. Si c'est le cas, la solution est
       * void (NULL).
       */
      if (ineqpar != NULL)
	{ /* Calcul du nombre d'inequations sur les parametres. Le format de
	   * matrice de la polylib permet les egalites, on doit donc les compter
	   * double quand il y en a.
	   */
	  Nm = ineqpar->NbRows ;
	  for (i=0;i<ineqpar->NbRows;i++)
#if defined(PIPMP_LINEAR_VALUE_IS_MP)
	    if (mpz_sgn(**(ineqpar->p + i)) == 0)
#else
	      if (**(ineqpar->p + i) == 0)
#endif
		Nm ++ ;
      
	  context = tabmp_Matrix2TableauMP(ineqpar,Nm,Np,0, Shift,Bg-Nn, Urs_parms);  
      
	  if (Nm)
	    { /* Traduction du format de matrice de la polylib vers celui de
	       * traitement de PipMP. Puis traitement proprement dit.
	       */
	      ctxt = pipmp_expanser(context, Np, Nm, Np+1, Np, 0, 0) ;
	      if (pipmp_traiter(ctxt, NULL, Pip_True, Np, 0, Nm, 0, -1) == 42)
		{
		  // Build a fake solution.
		  solution = (PipMPQuast*) malloc (sizeof(PipMPQuast));
		  solution->newparm = NULL;
		  solution->condition = NULL;
		  solution->list = (PipMPList*) malloc (sizeof(PipMPList));
		  solution->list->vector = NULL;
		  solution->list->next = NULL;
		  return solution;
		}
	      non_vide = pipmp_is_not_Nil(p) ;
	      solmp_reset(p) ;
	    }
	  else
	    non_vide = Pip_True ;
	}
      else
	{
	  Nm = 0 ;
	  PipMPMatrix* toto = pipmp_matrix_alloc(0, 2);
	  context = tabmp_Matrix2TableauMP(toto, 0, 0, 0, Shift, Bg, Urs_parms);
	  pipmp_matrix_free (toto);
	  /*     context = NULL ; /\* ATTENTION, en toute rigueur on devrait faire un */
	  /*                         * tabmp_Matrix2TableauMP d'une matrice 0*0. */
	  /* 			*\/ */
	  non_vide = Pip_True ;
	}
        
      if (pipmp_verbose > 0)
	fprintf(pipmp_dump, "%d %d %d %d %d %d\n",Nn,Np,Nl,Nm,Bg,options->Nq) ;
    
      /* S'il est possible de trouver une solution, on passe au traitement. */
      if (non_vide)
	{ ineq = tabmp_Matrix2TableauMP(inequnk,Nl,Nn,Nn, Shift,Bg, Urs_parms);
	pipmp_compa_count = 0 ;
	if (pipmp_traiter(ineq, context, options->Nq, Nn, Np, Nl, Nm, Bg) != 42)
	  {
	    if (options->Simplify)
	      solmp_simplify(xq) ;
	    q = solmp_hwm() ;
	    /* On traduit la solution du format de solution de PipMP vers un arbre
	     * de structures de type PipMPQuast.
	     */
	    solution = solmp_quast_edit(&xq, NULL, Bg-Nn-1, Urs_parms, solmp_flags);
	    solmp_reset(p) ;
	  }
	else
	  {
	    // Build a fake solution.
	    solution = (PipMPQuast*) malloc (sizeof(PipMPQuast));
	    solution->newparm = NULL;
	    solution->condition = NULL;
	    solution->list = (PipMPList*) malloc (sizeof(PipMPList));
	    solution->list->vector = NULL;
	    solution->list->next = NULL;
	    return solution;
	  }
	}
      else
	return NULL ;
      tabmp_reset(hq) ;
    }
  else
    return NULL ;
  
  return(solution) ;
}
