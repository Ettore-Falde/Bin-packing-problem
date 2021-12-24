# Libraries
library(Rglpk) 
library(slam)

# Mathprog
bins_problem <- '

/* The bin packing problem is aimed to pack I = {1,...,m}
elements of different volumes v[i] into bins each with a fixed capacity V
in a way that minimizes the number of bins used. */

  param m, integer, > 0;   /* number of elements to pack */
  set I := 1..m;           /* set of elements */
  param v{i in 1..m}, > 0; /* v[i] is volume of element i */
  param V, > 0;            /* bin capacity */

  /* worse scenario for the number of bins used */
  param k := m;
  
  
  
  /* We need to estimate an upper bound of the number of bins sufficient
   to contain all items. The number of items m can be used, however, it
   is not a good idea. To obtain a more suitable estimation an easy
   heuristic is used: we put items into a bin while it is possible, and
   if the bin is full, we use another bin. The number of bins used in
   this way gives us a more appropriate estimation. */

param z{i in I, j in 1..m} :=
/* z[i,j] = 1 if item i is in bin j, otherwise z[i,j] = 0 */

   if i = 1 and j = 1 then 1
   /* put item 1 into bin 1 */

   else if exists{jj in 1..j-1} z[i,jj] then 0
   /* if item i is already in some bin, do not put it into bin j */

   else if sum{ii in 1..i-1} v[ii] * z[ii,j] + v[i] > V then 0
   /* if item i does not fit into bin j, do not put it into bin j */

   else 1;
   /* otherwise put item i into bin j */

check{i in I}: sum{j in 1..m} z[i,j] = 1;
/* each item must be exactly in one bin */

check{j in 1..m}: sum{i in I} v[i] * z[i,j] <= V;
/* no bin must be overflowed */

param n := sum{j in 1..m} if exists{i in I} z[i,j] then 1;
/* determine the number of bins used by the heuristic; obviously it is
   an upper bound of the optimal solution */

display n;



/* set of bins */
set J := 1..n;


/* x[i,j] = 1 means elemnt i is in bin j */
var x{i in I, j in J}, binary;


/* y[j] = 1 bin j is used */
var y{j in J}, binary;



/* each element of the set of elements must be contained in exactly one bin */
s.t. one{i in I}: sum{j in J} x[i,j] = 1;


/* the content of the bin must not exceed its capacity */
s.t. lim{j in J}: sum{i in I} v[i] * x[i,j] <= V * y[j];

/* the minimization of the number of bins used is the objective function of this problem */
minimize obj: sum{j in J} y[j];

  data;

param m := 10;

param v := 1 0.61, 2 0.96, 3 0.95, 4 0.91, 5 0.13, 6 0.53, 7 0.53, 8 0.05, 9 0.65, 10 0.66;

param V := 1;

end;

'


# Solve function
solve_lp <- function(model, type = "MathProg"){
  
  writeLines(model, "model.mod")
  
  lp_model <- Rglpk_read_file("model.mod", type = type)
  
  lp_solution <- Rglpk_solve_LP(obj = lp_model$objective,
                                mat = lp_model$constraints[[1]],
                                dir = lp_model$constraints[[2]],
                                rhs = lp_model$constraints[[3]],
                                types = lp_model$types,
                                bounds = lp_model$bounds,
                                max = lp_model$maximum)
  unlink("model.mod")
  
  if(lp_solution$status == 0){
    sol <- lp_solution$solution
    names(sol) <- attr(lp_model, "objective_vars_names")
    obj = lp_solution$optimum
    return(list(sol = sol, obj = obj, success = TRUE))
  }else
    return(success = FALSE)
  
}


# Results
results <- solve_lp( bins_problem, type = "MathProg")


# Show them clearly
packages <- sprintf("package[%s]",seq(1:10))
names <- as.array(packages)
names

bins <- sprintf("bin[%s]",seq(1:results$obj))
bin_names <- as.array(bins)
bin_names

mat <- matrix(, nrow = results$obj, ncol = 10, byrow = TRUE,
              dimnames = list(bin_names, names))



for(i in 1:length(mat)){
  if(results$sol[i] == 1) {
    mat[i] <- 1
  } else {
    mat[i] <- 0
  }
}

mat <- as.table(mat)
knitr::kable(mat)