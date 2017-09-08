(** Simple 1-D harmonic oscillator, with perturbations.  We use the Lagrangian 

    L = 0.5*qdot^2 - 0.5*q^2 - eps*q^3/3*)

open Newton
open Mapping

(** States. *)
type state = {q : float; p : float}

(** [evolve dt state] returns the evolution with eps = 0 of [state] by
    [dt]. *)
let evolve dt {q = q0; p = p0} = 
  let a = q0 and 
      b = p0 in 
  let cdt = cos dt and 
      sdt = sin dt in 
  {q = a *. cdt +. b *. sdt;
   p = b *. cdt -. a *. sdt}

(** [energy state] returns the unperturbed energy of [state]. *)
let energy {q = q; p = p} = 
  0.5*.(p*.p +. q*.q)

(** [kick_pert eps dt state] implements the "kick" operator by the
    perturbing term with magnitude [eps] over a time [dt]. *)
let kick_pert eps dt {q = q0; p = p0} = 
  {q = q0;
   p = p0 -. eps*.dt*.q0*.q0}

(** [energy_pert eps state] returns the full energy of [state]
    (including perturbation). *)
let energy_pert eps ({q = q; p = p} as s) = 
  energy s +. eps*.q*.q*.q/.3.0

(** [dkd eps] returns a "drift-kick-drift" advancer of states, using
    [eps] for the magnitude of the perturbation. *)
let dkd eps = 
  mapping_advance [|0.5; 1.0; 0.5|] evolve (kick_pert eps)

(** [kdk eps] is a "kick-drift-kick" advancer. *)
let kdk eps = 
  mapping_advance [|0.5; 1.0; 0.5|] (kick_pert eps) evolve 

(** [s4b eps] is the S4B from Chambers and Murison arXiv:astro-ph/9910263 *)
let s4b eps = 
  mapping_advance [|1.0/.6.0; 0.5; 2.0/.3.0; 0.5; 1.0/.6.0|] (kick_pert eps) 
    evolve

(** [s6b eps] is the S6B from Chambers and Murison *)
let s6b eps = 
  mapping_advance 
    [|1.0/.12.0; 0.5*.(1.0-.1.0/.(sqrt 5.0)); 
      5.0/.12.0; 1.0/.(sqrt 5.0); 5.0/.12.0; 
      0.5*.(1.0-.1.0/.(sqrt 5.0)); 1.0/.12.0|]
    (kick_pert eps)
    evolve

let csc x = 1.0 /. (sin x)
let sec x = 1.0 /. (cos x)
let cot x = 1.0 /. (tan x)

(** [pow_int x n] raises [x] to the integer power [n] efficiently. *)
let rec pow_int x n = 
  if n = 0 then 
    1.0
  else if n = 1 then 
    x 
  else if n mod 2 = 0 then 
    let tmp = pow_int x (n/2) in 
    tmp *. tmp
  else
    x *. (pow_int x (n-1))

let square x = pow_int x 2
let cube x = pow_int x 3

let q1_eqn eps t0 t1 q0 p0 q1 =
  p0 +. q1*.(csc (t0 -. t1)) +. 
    ((2.0*.(square q0) +. 2.0*.q0*.q1 +. (square q1))*.eps*.
       (square (sec ((t0 -. t1)*.0.5)))*.(tan ((t0 -. t1)*.0.5)))/.6. +. 
    (q0*.(cot (t0 -. t1))*.((-3.0) +. q0*.eps*.(square (tan ((t0 -. t1)*.0.5)))))/.3.

let dq1_eqn eps t0 t1 q0 q1 = 
  (csc (t0 -. t1)) +. ((q0 +. q1)*.eps*.(square (sec ((t0 -. t1)*.0.5)))*.
      (tan ((t0 -. t1)*.0.5)))/.3.

let p1_eqn eps t0 t1 q0 q1 = 
  (6.0*.q0*.(csc (t0 -. t1)) +. 
     ((square q0) +. 2.0*.q0*.q1 +. 2.0*.(square q1))*.eps*.
     (square (sec (0.5*.(t0 -. t1))))*.(tan (0.5*.(t0 -. t1))) +. 
     2.0*.q1*.(cot (t0 -. t1))*.
     ((-3.0) +. q1*.eps*.(square (tan (0.5*.(t0 -. t1))))))/.6.

let make_variational_advancer q1_eqn dq1_eqn p1_eqn = 
  fun eps dt ({q = q0; p = p0} as s) ->
    let err = 100.0*.epsilon_float and 
        t0 = 0.0 and 
        t1 = dt in 
    let {q = q1_guess} = evolve dt s in 
    let q1 = 
      newton_solve 100 err err err (q1_eqn eps t0 t1 q0 p0)
        (dq1_eqn eps t0 t1 q0) q1_guess in 
    {q = q1;
     p = p1_eqn eps t0 t1 q0 q1}
    

(** [variational_advance eps dt state] advances [state] by [dt] using
    [eps] as the strength of the perturbation using the variational
    algorithm.  *)
let variational_advance = make_variational_advancer q1_eqn dq1_eqn p1_eqn

let q1_eqn_dkd eps t0 t1 q0 p0 q1 = 
  p0 -. q0*.(cot (t0 -. t1)) +. q1*.(csc (t0 -. t1)) +. 
   ((square (q0 +. q1))*.(t0 -. t1)*.eps*.(cube (sec ((t0 -. t1)*.0.5))))/.8.

let dq1_eqn_dkd eps t0 t1 q0 q1 = 
  (csc (t0 -. t1)) +. ((q0 +. q1)*.(t0 -. t1)*.eps*.(cube (sec ((t0 -. t1)/.2.0))))/.4.

let p1_eqn_dkd eps t0 t1 q0 q1 = 
  ~-.(q1*.(cot (t0 -. t1))) +. q0*.(csc (t0 -. t1)) +. 
   ((square (q0 +. q1))*.(t0 -. t1)*.eps*.(cube (sec ((t0 -. t1)/.2.0))))/.8.

let variational_advance_dkd = make_variational_advancer q1_eqn_dkd dq1_eqn_dkd p1_eqn_dkd

let q1_eqn_gl3 eps t0 t1 q0 p0 q1 = 
  p0 -. q0*.(cot (t0 -. t1)) +. q1*.(csc (t0 -. t1)) +. 
    ((t0 -. t1)*.eps*.(2.0*.(square q0) +. 
                         (square (q0 +. q1))*.(cube (sec ((t0 -. t1)*.0.5)))))/.12.
    
let dq1_eqn_gl3 eps t0 t1 q0 q1 = 
  (csc (t0 -. t1)) +. ((q0 +. q1)*.(t0 -. t1)*.eps*.(cube (sec ((t0 -. t1)*.0.5))))/.6.

let p1_eqn_gl3 eps t0 t1 q0 q1 = 
  ~-.(q1*.(cot (t0 -. t1))) +. q0*.(csc (t0 -. t1)) +. 
   ((t0 -. t1)*.eps*.(2.0*.(square q1) +. 
        (square (q0 +. q1))*.(cube (sec ((t0 -. t1)*.0.5)))))/.12.

let variational_advance_gl3 = make_variational_advancer q1_eqn_gl3 dq1_eqn_gl3 p1_eqn_gl3
  
