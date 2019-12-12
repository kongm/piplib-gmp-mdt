/******************************************************************************
 *                     PIP : Parametric Integer Programming                   *
 ******************************************************************************
 *                                   sol.h                                    *
 ******************************************************************************
 *                                                                            *
 * Copyright Paul Feautrier, 1988, 1993, 1994, 1996, 2002                     *
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
 * Written by Paul Feautrier and Cedric Bastoul                               *
 *                                                                            *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <piplib-gmp/piplib-gmp.h>


#ifndef PIPMP_LINEAR_VALUE_IS_MP
# define PIPMP_LINEAR_VALUE_IS_MP
#endif


extern long int pipmp_cross_product, pipmp_limit;
extern int pipmp_verbose;
extern FILE *pipmp_dump;

struct S
    {int flags;
     EntierMP param1, param2;
    };

#define Free 0
#define Nil  1
#define If   2
#define List 3
#define Form 4
#define New  5
#define Div  6
#define Val  7
#define Error 8

struct S * solmp_space;
static int solmp_free;

#if !defined(PIPMP_LINEAR_VALUE_IS_MP)
EntierMP mod(EntierMP, EntierMP);

EntierMP pgcd(EntierMP x, EntierMP y)
{EntierMP r;
 while(y)
     {r = mod(x, y);
      x = y;
      y = r;
     }
 return(x>= 0? x : -x);
}
#endif

void solmp_init(void)
{
 solmp_free = 0;
 solmp_space = (struct S *)malloc(SOL_SIZE*sizeof(struct S)) ;
}

void solmp_close(void)
{
 free(solmp_space) ;
}

int solmp_hwm()
{
 return(solmp_free);
}

void solmp_reset(p)
int p;
{int i;
 if(p<0 || p>=SOL_SIZE)
     {fprintf(stderr, "Syserr : solmp_reset : Memory allocation error\n");
      exit(40);
     }
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 for(i=p; i<solmp_free; i++){
   mpz_clear(solmp_space[i].param1);
   mpz_clear(solmp_space[i].param2);
 }
 #endif
 solmp_free = p;
}

struct S *solmp_alloc(void)
{struct S *r;
 r = solmp_space + solmp_free;
 r->flags = Free;
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_init_set_si(r->param1,0);
 mpz_init_set_si(r->param2,0);
 #else
 r->param1 = r->param2 = 0;
 #endif
 solmp_free++;
 if(solmp_free >= SOL_SIZE)
     {fprintf(stderr, "The solution is too complex! : sol\n");
      exit(26);
     }
     return(r);
}

void solmp_nil(void)
{
 struct S * r;
 r = solmp_alloc();
 r -> flags = Nil;
 if(pipmp_verbose > 0)
   {fprintf(pipmp_dump, "\nNil");
    fflush(pipmp_dump);
  }
}

void solmp_error(int c)
{
 struct S *r;
 r = solmp_alloc();
 r->flags = Nil;
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_set_si(r->param1, c);
 #else
 r->param1 = c;
 #endif
 if(pipmp_verbose > 0) {
     fprintf(pipmp_dump, "Erreur %d\n", c);
     fflush(pipmp_dump);
     }
}

int pipmp_is_not_Nil(p)
int p;
{
 return(solmp_space[p].flags != Nil);
}

void solmp_if(void)
{
 struct S *r;
 r = solmp_alloc();
 r -> flags = If;
 if(pipmp_verbose > 0) {
     fprintf(pipmp_dump, "\nIf ");
     fflush(pipmp_dump);
   }
}

void solmp_list(n)
int n;
{struct S * r;
 r = solmp_alloc();
 r->flags = List;
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_set_si(r->param1, n);
 #else
 r->param1 = n;
 #endif
 if(pipmp_verbose > 0) {
     fprintf(pipmp_dump, "\nList %d ", n);
     fflush(pipmp_dump);
}
}

void solmp_forme(l)
int l;
{
 struct S *r;
 r = solmp_alloc();
 r -> flags = Form;
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_set_ui(r -> param1, l);
 #else
 r -> param1 = l;
 #endif
 if(pipmp_verbose > 0) {
     fprintf(pipmp_dump, "\nForme %d ", l);
     fflush(pipmp_dump);
   }
}

void solmp_new(k)
int k;
{
 struct S *r;
 r = solmp_alloc();
 r -> flags = New;
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_set_ui(r -> param1, k);
 #else
 r -> param1 = k;
 #endif
 if(pipmp_verbose > 0) {
     fprintf(pipmp_dump, "New %d ", k);
     fflush(pipmp_dump);
   }
}

void solmp_div()
{
 struct S *r;
 r = solmp_alloc();
 r -> flags = Div;
 if(pipmp_verbose > 0) {
     fprintf(pipmp_dump, "Div ");
     fflush(pipmp_dump);
   }
}

void solmp_val(n, d)
EntierMP n, d;
{
 struct S *r;
 r = solmp_alloc();
 r -> flags = Val;
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_set(r -> param1, n);
 mpz_set(r -> param2, d);
 #else
 r -> param1 = n;
 r -> param2 = d;
 #endif
 if(pipmp_verbose > 0) {
   fprintf(pipmp_dump, "val(");
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   mpz_out_str(pipmp_dump, 10, n);
   fprintf(pipmp_dump, "/");
   mpz_out_str(pipmp_dump, 10, d);
   #else
   fprintf(pipmp_dump, FORMAT, n);
   fprintf(pipmp_dump, "/");
   fprintf(pipmp_dump, FORMAT, d);
   #endif
   fprintf(pipmp_dump, ") ");
   fflush(pipmp_dump);
  }
}


/* a` partir d'un point de la solution, sauter un objet
bien forme' ainsi qu'un e'ventuel New et pointer sur l'objet
suivant */
static int skip(int);

static
int skip_New (int i)
{
 if(solmp_space[i].flags != New) return i;
 i = skip(i+1);      /* sauter le Div */
 return i;
}
/* au lancement, i indexe une cellule qui est la te^te d'un objet.
   la valeur retourne'e est la te^te de l'objet qui suit. Les
   objets de type New sont e'limine's                                */

static
int skip (int i)
{int n, f;
 while((f = solmp_space[i].flags) == Free || f == Error) i++;
 switch (solmp_space[i].flags) {
 case Nil : case Val : i++; break;
 case New : i = skip_New(i); break;
 case If : i = skip(i+1);        /* sauter le pre'dicat */
	   i = skip(i);          /* sauter le vrai */
	   i = skip(i); break;   /* sauter le faux */
 case List : case Form :
           #if defined(PIPMP_LINEAR_VALUE_IS_MP)
           n = mpz_get_si(solmp_space[i].param1);
           #else
           n = solmp_space[i].param1;
           #endif
	   i++;
	   while(n--) i = skip(i);
	   break;
 case Div : i = skip(i+1);       /* sauter la forme */
	    i = skip(i);         /* sauter le diviseur */
	    break;
 default : fprintf(stderr,
	      "Syserr : skip : unknown %d\n", solmp_space[i].flags);
 }
 return skip_New(i);
}
/* simplification de la solution : e'limination des constructions
   (if p () ()). N'est en service qu'en pre'sence de l'option -z */

void solmp_simplify(int i)
{int j, k, l;
 if(solmp_space[i].flags == If) {
     j = skip(i+1);        /* j : debut de la partie vraie */
     k = skip(j);          /* k : debut de la partie fausse */
     solmp_simplify(k);
     solmp_simplify(j);
     if(solmp_space[j].flags == Nil && solmp_space[k].flags == Nil) {
	 solmp_space[i].flags = Nil;
	 if (k >= solmp_free - 1) 
	    solmp_reset(i+1);
	 else for(l = i+1; l<=k; l++) solmp_space[l].flags = Free;
       }
   }

}
/* e'dition de la solution */

int solmp_edit(FILE *foo, int i)
{int j, n;
 struct S *p;
 EntierMP N, D, d;
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_init(N);
 mpz_init(D);
 mpz_init(d);
 #endif
 
 p = solmp_space + i;
 for(;;) {
   if(p->flags == Free) {
     p++;
     i++;
     continue;
   }
   if(p->flags == New) {
     #if defined(PIPMP_LINEAR_VALUE_IS_MP)
     n = mpz_get_si(p->param1);
     #else
     n = p->param1;
     #endif
     fprintf(foo, "(newparm %d ", n);
     if(pipmp_verbose>0)fprintf(pipmp_dump, "(newparm %d ", n);
     i = solmp_edit(foo, ++i);
     p = solmp_space +i;
     fprintf(foo, ")\n");
     if(pipmp_verbose>0)fprintf(pipmp_dump, ")\n");
     continue;
   }
   break;
 }
 switch(p->flags){
 case Nil : fprintf(foo, "()\n");
   if(pipmp_verbose>0)fprintf(pipmp_dump, "()\n");
   i++; break;
 case Error :
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   fprintf(foo, "Error %ld\n", mpz_get_si(p->param1));
   if(pipmp_verbose>0)
   fprintf(pipmp_dump, "Error %ld\n", mpz_get_si(p->param1));
   #else
   fprintf(foo, "Error %d\n", p->param1);
   if(pipmp_verbose>0)
   fprintf(pipmp_dump, "Error %d\n", p->param1);
   #endif
   i++; break;
 case If  : fprintf(foo, "(if ");
   if(pipmp_verbose>0)fprintf(pipmp_dump, "(if ");
   i = solmp_edit(foo, ++i);
   i = solmp_edit(foo, i);
   i = solmp_edit(foo, i);
   fprintf(foo, ")\n");
   if(pipmp_verbose>0)fprintf(pipmp_dump, ")\n");
   break;
 case List: fprintf(foo, "(list ");
   if(pipmp_verbose>0)fprintf(pipmp_dump, "(list ");
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   n = mpz_get_si(p->param1);
   #else
   n = p->param1;
   #endif
   i++;
   while(n--) i = solmp_edit(foo, i);
   fprintf(foo, ")\n");
   if(pipmp_verbose>0)fprintf(pipmp_dump, ")\n");
   break;
 case Form: fprintf(foo, "#[");
   if(pipmp_verbose>0)fprintf(pipmp_dump, "#[");
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   n = mpz_get_si(p->param1);
   #else
   n = p->param1;
   #endif
   for(j = 0; j<n; j++){
     i++; p++;
     #if defined(PIPMP_LINEAR_VALUE_IS_MP)
     mpz_set(N, p->param1); mpz_set(D, p->param2);
     mpz_gcd(d, N, D);
     if(mpz_cmp(d, D) == 0){
       putc(' ', foo);
       mpz_divexact(N, N, d);
       mpz_out_str(foo, 10, N);
       if(pipmp_verbose>0){
         putc(' ', pipmp_dump);
         mpz_out_str(pipmp_dump, 10, N);
       }
     }
     else{
       mpz_divexact(N, N, d);
       mpz_divexact(D, D, d);
       putc(' ', foo);
       mpz_out_str(foo, 10, N);
       putc('/', foo);
       mpz_out_str(foo, 10, D);
       if(pipmp_verbose>0){
         putc(' ', pipmp_dump);
         mpz_out_str(pipmp_dump, 10, N);
         putc('/', pipmp_dump);
         mpz_out_str(pipmp_dump, 10, D);
       }
     }
     #else
     N = p->param1; D = p->param2;
     d = pgcd(N, D);
     if(d == D){
       putc(' ', foo);
       fprintf(foo, FORMAT, N/d);
       if(pipmp_verbose>0){
	 putc(' ', pipmp_dump);
	 fprintf(pipmp_dump, FORMAT, N/d);
       }
     }
     else{
       putc(' ', foo);
       fprintf(foo,FORMAT,N/d);
       fprintf(foo,"/");
       fprintf(foo,FORMAT, D/d);
       if(pipmp_verbose>0){
	 putc(' ', pipmp_dump);
	 fprintf(pipmp_dump,FORMAT,N/d);
	 fprintf(pipmp_dump,"/");
	 fprintf(pipmp_dump,FORMAT, D/d);
       }
     }
     #endif
   }
   fprintf(foo, "]\n");
   if(pipmp_verbose>0)fprintf(pipmp_dump, "]\n");
   i++;
   break;
 case Div : fprintf(foo, "(div ");
   if(pipmp_verbose>0)fprintf(pipmp_dump, "(div ");
   i = solmp_edit(foo, ++i);
   i = solmp_edit(foo, i);
   fprintf(foo, ")\n");
   if(pipmp_verbose>0)fprintf(pipmp_dump, ")\n");
   break;
 case Val :
   #if defined(PIPMP_LINEAR_VALUE_IS_MP)
   mpz_set(N, p->param1); mpz_set(D, p->param2);
   mpz_gcd(d, N, D);
   if(mpz_cmp(d, D) == 0){
     mpz_divexact(N, N, d);
     putc(' ', foo);
     mpz_out_str(foo, 10, N);
     if(pipmp_verbose>0){
       putc(' ', pipmp_dump);
       mpz_out_str(pipmp_dump, 10, N);
     }
   }
   else{
     mpz_divexact(N, N, d);
     mpz_divexact(D, D, d);
     putc(' ', foo);
     mpz_out_str(foo, 10, N);
     fprintf(foo, "/");
     mpz_out_str(foo, 10, D);
     if(pipmp_verbose>0){
       putc(' ', pipmp_dump);
       mpz_out_str(pipmp_dump, 10, N);
       fprintf(pipmp_dump, "/");
       mpz_out_str(pipmp_dump, 10, D);
     }
   }
   #else
   N = p->param1; D = p->param2;
   d = pgcd(N, D);
   if(d == D){putc(' ', foo);
   fprintf(foo, FORMAT, N/d);
   if(pipmp_verbose>0)
     {putc(' ', pipmp_dump);
     fprintf(pipmp_dump, FORMAT, N/d);
     }
   }
   else{putc(' ', foo);
   fprintf(foo, FORMAT, N/d);
   fprintf(foo, "/");
   fprintf(foo, FORMAT, D/d);
   if(pipmp_verbose>0)
     {putc(' ', pipmp_dump);
     fprintf(pipmp_dump, FORMAT, N/d);
     fprintf(pipmp_dump, "/");
     fprintf(pipmp_dump, FORMAT, D/d);
     }
   }
   #endif
   i++;
   break;
 default  : fprintf(foo, "Inconnu : sol\n");
   if(pipmp_verbose>0)fprintf(pipmp_dump, "Inconnu : sol\n");
 }
 #if defined(PIPMP_LINEAR_VALUE_IS_MP)
 mpz_clear(d);
 mpz_clear(D);
 mpz_clear(N);
 #endif
 return(i);
}


/* Fonction solmp_vector_edit :
 * Cette fonction a pour but de placer les informations correspondant
 * a un Vector dans la grammaire dans une structure de type PipMPVector. Elle
 * prend en parametre un pointeur vers une case memoire contenant le
 * numero de cellule du tableau solmp_space a partir de laquelle on doit
 * commencer la lecture des informations. Elle retourne un pointeur vers
 * une structure de type PipMPVector contenant les informations de ce Vector.
 * Premiere version : Ced. 20 juillet 2001. 
 */
PipMPVector * solmp_vector_edit(int *i, int Bg, int Urs_p, int flags)
{ int j, k, n, unbounded  = 0, first_urs;
  struct S *p ;
  EntierMP N, D, d ;
  PipMPVector * vector ;

  mpzvalue_init(N) ;
  mpzvalue_init(D) ;
  mpzvalue_init(d) ;
  
  vector = (PipMPVector *)malloc(sizeof(PipMPVector)) ;
  if (vector == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  p = solmp_space + (*i) ;
  n = MPZVALUE_TO_INT(p->param1);
  if (flags & SOL_REMOVE)
    --n;
  n -= Urs_p;
  first_urs = Urs_p + (Bg >= 0);
  vector->nb_elements = n ;
  vector->the_vector = (EntierMP *)malloc(sizeof(EntierMP)*n) ;
  if (vector->the_vector == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  vector->the_deno = (EntierMP *)malloc(sizeof(EntierMP)*n) ;
  if (vector->the_deno == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  
  for (j=0, k=0; k < n; j++) {
    (*i)++ ;
    p++ ;

    mpzvalue_assign(N, p->param1);
    mpzvalue_assign(D, p->param2);
    mpzvalue_gcd(d, N, D);

    if ((flags & SOL_SHIFT) && j == Bg) {
      mpzvalue_subtract(N, N, D);   /* subtract 1 */
      if (mpzvalue_notzero_p(N))
	unbounded = 1;
    }

    if ((flags & SOL_REMOVE) && j == Bg)
      continue;

    if (first_urs <= j && j < first_urs+Urs_p)
      continue;

    mpzvalue_init(vector->the_vector[k]);
    mpzvalue_divexact(vector->the_vector[k], N, d);
    if (flags & SOL_NEGATE)
      mpzvalue_oppose(vector->the_vector[k], vector->the_vector[k]);
    mpzvalue_init(vector->the_deno[k]);
    if (mpzvalue_eq(d, D))
      mpzvalue_assign(vector->the_deno[k], pipmp_UN);
    else
      mpzvalue_divexact(vector->the_deno[k], D, d);
    ++k;
  }
  if (unbounded)
    for (k=0; k < n; k++)
      mpzvalue_assign(vector->the_deno[k], pipmp_ZERO);
  (*i)++ ;

  mpzvalue_clear(d);
  mpzvalue_clear(D);
  mpzvalue_clear(N);

  return(vector) ;
}


/* Fonction solmp_newparm_edit :
 * Cette fonction a pour but de placer les informations correspondant
 * a un Newparm dans la grammaire dans une structure de type PipMPNewparm. Elle
 * prend en parametre un pointeur vers une case memoire contenant le
 * numero de cellule du tableau solmp_space a partir de laquelle on doit
 * commencer la lecture des informations. Elle retourne un pointeur vers
 * une structure de type PipMPNewparm contenant les informations de ce Newparm.
 * Premiere version : Ced. 18 octobre 2001. 
 */
PipMPNewparm * solmp_newparm_edit(int *i, int Bg, int Urs_p, int flags)
{ struct S * p ;
  PipMPNewparm * newparm, * newparm_first = NULL, * newparm_now = NULL;

  /* On place p au lieu de lecture. */
  p = solmp_space + (*i) ;

  do {
    /* On passe le New et le Div pour aller a Form et lire le VECTOR. */
    (*i) += 2 ;

    newparm = (PipMPNewparm *)malloc(sizeof(PipMPNewparm)) ;
    if (newparm == NULL)
    { fprintf(stderr, "Memory Overflow.\n") ;
      exit(1) ;
    }
    newparm->vector = solmp_vector_edit(i, Bg, Urs_p, flags);
    newparm->rank = MPZVALUE_TO_INT(p->param1);
    /* On met p a jour pour lire le denominateur (un Val de param2 pipmp_UN). */
    p = solmp_space + (*i) ;
    mpzvalue_init_set(newparm->deno, p->param1);
    if (flags & SOL_REMOVE)
      newparm->rank--;
    newparm->rank -= Urs_p;
    newparm->next = NULL ;

    if (newparm_now)
      newparm_now->next = newparm;
    else
      newparm_first = newparm;
    newparm_now = newparm ;
    if (pipmp_verbose)
    { fprintf(pipmp_dump,"\n(newparm ") ;
      fprintf(pipmp_dump,PIPMP_FORMAT,newparm->rank) ;
      fprintf(pipmp_dump," (div ") ;
      pipmp_vector_print(pipmp_dump,newparm->vector) ;
      fprintf(pipmp_dump," ") ;
      #if defined(PIPMP_LINEAR_VALUE_IS_MP)
      mpz_out_str(pipmp_dump,10,newparm->deno) ;
      #else
      fprintf(pipmp_dump,FORMAT,newparm->deno) ;
      #endif
      fprintf(pipmp_dump,"))") ;
    }
  
    /* On passe aux elements suivants. */
    (*i) ++ ;
    p = solmp_space + (*i) ;
  } while (p->flags == New);

  return newparm_first;
}


/* Fonction solmp_list_edit :
 * Cette fonction a pour but de placer les informations correspondant
 * a une List dans la grammaire dans une structure de type PipMPList. Elle
 * prend en parametre un pointeur vers une case memoire contenant le
 * numero de cellule du tableau solmp_space a partir de laquelle on doit
 * commencer la lecture des informations. Elle retourne un pointeur vers
 * une structure de type PipMPList contenant les informations de cette List.
 * Premiere version : Ced. 18 octobre 2001. 
 * 16 novembre 2005 : Ced. Prise en compte du cas 0 éléments, avant impossible.
 */
PipMPList * solmp_list_edit(int *i, int nb_elements, int Bg, int Urs_p, int flags)
{ PipMPList * list, * list_new, * list_now ;
  
  /* Pour le premier element. */
  list = (PipMPList *)malloc(sizeof(PipMPList)) ;
  if (list == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  list->next = NULL ;
  
  if (nb_elements == 0)
  { list->vector = NULL ;
    return(list) ;
  }
  
  list->vector = solmp_vector_edit(i, Bg, Urs_p, flags);

  list_now = list ;
  if (pipmp_verbose)
  { fprintf(pipmp_dump,"\n(list ") ;
    pipmp_vector_print(pipmp_dump,list->vector) ;
  }
  nb_elements-- ;

  /* Pour les elements suivants. */
  while (nb_elements--)
  { list_new = (PipMPList *)malloc(sizeof(PipMPList)) ;
    if (list_new == NULL)
    { fprintf(stderr, "Memory Overflow.\n") ;
      exit(1) ;
    }
    list_new->vector = solmp_vector_edit(i, Bg, Urs_p, flags);
    list_new->next = NULL ;
		    
    if (pipmp_verbose)
    { fprintf(pipmp_dump,"\n") ;
      pipmp_vector_print(pipmp_dump,list_new->vector) ;
    }
    list_now->next = list_new ;
    list_now = list_now->next ;
  }
  if (pipmp_verbose)
  fprintf(pipmp_dump,"\n)") ;
  
  return(list) ;
}


/* Fonction solmp_quast_edit :
 * Cette fonction a pour but de placer les informations de la solution
 * (qui sont contenues dans le tableau solmp_space) dans une structure de
 * type PipMPQuast en vue d'une utilisation directe de la solution par une
 * application exterieure. Elle prend en parametre un pointeur vers une
 * case memoire contenant le numero de cellule du tableau solmp_space
 * a partir de laquelle on doit commencer la lecture des informations. Elle
 * recoit aussi l'adresse du PipMPQuast qui l'a appelle (pour le champ father).
 * Elle retourne un pointeur vers une structure de type PipMPQuast qui
 * contient toutes les informations sur la solution (sous forme d'arbre).
 * Remarques : cette fonction lit les informations comme elles doivent
 * se presenter a la fin du traitement. Elle respecte scrupuleusement
 * la grammaire attendue et n'accepte de passer des cellules a Free
 * qu'entre une des trois grandes formes (if, list ou suite de newparm).
 * 20  juillet 2001 : Premiere version, Ced. 
 * 31  juillet 2001 : Ajout du traitement de l'option pipmp_verbose = code*2 :0( 
 * 18  octobre 2001 : Grands changements dus a l'eclatement de la structure
 *                    PipMPVector en PipMPVector, PipMPNewparm et PipMPList, et
 *                    eclatement de la fonction avec solmp_newparm_edit et
 *                    solmp_list_edit.
 * 16 novembre 2005 : (debug) Même si une liste est vide il faut la créer pour
 *                    afficher plus tard le (list), repéré par Sven Verdoolaege.
 */
PipMPQuast *solmp_quast_edit(int *i, PipMPQuast *father, int Bg, int Urs_p, int flags)
{ int nb_elements ;
  struct S * p ;
  PipMPQuast * solution ;
  PipMPList * list_new, * list_now ;
  PipMPNewparm * newparm_new, * newparm_now ;
    
  /* On place p au lieu de lecture. */
  p = solmp_space + (*i) ;
  /* En cas d'utilisation de l'option de simplification, une plage de
   * structures S peut avoir les flags a Free. On doit alors les passer.
   */
  while (p->flags == Free)
  { p ++ ;
    (*i) ++ ;
  }
  
  solution = (PipMPQuast *)malloc(sizeof(PipMPQuast)) ;
  if (solution == NULL)
  { fprintf(stderr, "Memory Overflow.\n") ;
    exit(1) ;
  }
  solution->newparm = NULL ;
  solution->list = NULL ;
  solution->condition = NULL ;
  solution->next_then = NULL ;
  solution->next_else = NULL ;
  solution->father = father ;
  
  /* On peut commencer par une chaine de nouveaux parametres... */
  if (p->flags == New)
  { solution->newparm = solmp_newparm_edit(i, Bg, Urs_p, flags & SOL_REMOVE);
    p = solmp_space + (*i) ;
  }
  
  /* ...ensuite soit par une liste (vide ou non) soit par un if. */
  (*i)++ ; /* Factorise de List, Nil et If. */
  switch (p->flags)
  { case List : 
                #if defined(PIPMP_LINEAR_VALUE_IS_MP)
                nb_elements = mpz_get_si(p->param1) ;
                #else
                nb_elements = p->param1 ;
                #endif
                solution->list = solmp_list_edit(i, nb_elements, Bg, Urs_p, flags);
		break ;
    case Nil  : if (pipmp_verbose)
		fprintf(pipmp_dump,"\n()") ;
                break ;
    case If   : solution->condition = 
			    solmp_vector_edit(i, Bg, Urs_p, flags & SOL_REMOVE);
                if (pipmp_verbose)
		{ fprintf(pipmp_dump,"\n(if ") ;
                  pipmp_vector_print(pipmp_dump,solution->condition) ;
                }
		solution->next_then = solmp_quast_edit(i, solution, Bg, Urs_p, flags);
                solution->next_else = solmp_quast_edit(i, solution, Bg, Urs_p, flags);
                if (pipmp_verbose)
		fprintf(pipmp_dump,"\n)") ;
                break ;
    default   : fprintf(stderr,"\nAie !!! Flag %d inattendu.\n",p->flags) ;
                if (pipmp_verbose)
		fprintf(pipmp_dump,"\nAie !!! Flag %d inattendu.\n",p->flags) ;
                exit(1) ;
  }
  
  return(solution) ;
}
