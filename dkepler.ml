open Dnewton
open Diff.DFloat

let eps = 4.0*.epsilon_float

let abs_float x = (sqrt (x*x))

(** [kepler_dE n dt ec es] computes the change in eccentric anamoly.
    [ec] and [es] are e*cos(E0) and e*sin(E0), [n] is the mean motion,
    and [dt] is the corresponding time interval.*)
let kepler_dE n dt ec es = 
  let abs_ec = abs_float ec and 
      abs_es = abs_float ((const 2.0)*es) in 
  let xmin = n*dt - abs_ec - abs_es and 
      xmax = n*dt + abs_ec + abs_es in
  bracket_solve 20 eps eps eps 
    (fun x -> x - ec*(sin x) + es*((const 1.0) - (cos x)) - n*dt)
    xmin xmax

let dot r0 v0 = 
  let sum = ref (const 0.0) in 
  for i = 0 to 2 do 
    sum := !sum + r0.(i)*v0.(i)
  done;
  !sum

let norm r0 = sqrt (dot r0 r0)

let ecc mu r0 v0 = 
  let rnrm0 = norm r0 and 
      v02 = dot v0 v0 and 
      u = dot r0 v0 in 
  let a = (const 1.0)/((const 2.0)/rnrm0 - v02/mu) in 
  if a < (const 0.0) then raise (Failure "ecc: does not work on unbound orbits");
  let n = sqrt (mu / (a*a*a)) in 
  let ec = (const 1.0) - rnrm0/a and 
      es = u / (n * a * a) in 
  sqrt (ec*ec + es*es)
  

(** [kepler_advance mu dt r0 v0] returns the [(r,v)] pair which
    corresponds to advancing the kepler orbit with initial condition
    [(r0,v0)] by [dt].  [mu] is the parameter in the Lagrangian: L =
    1/2 v^2 + [mu]/r *)
let kepler_advance mu dt r0 v0 = 
  let rnrm0 = norm r0 and 
      v02 = dot v0 v0 and 
      u = dot r0 v0 in 
  let a = (const 1.0)/((const 2.0)/rnrm0 - v02/mu) in 
  if a < (const 0.0) then raise (Failure "kepler_advance: does not work on unbound orbits.");
  let n = sqrt (mu / (a*a*a)) in 
  let ec = (const 1.0) - rnrm0/a and 
      es = u / (n * a * a) in 
  let dE = kepler_dE n dt ec es in 
  let cosdE = cos dE and 
      sindE = sin dE in 
  let aor = (const 1.0) / ((const 1.0) - ec*cosdE + es*sindE) in 
  let f = a / rnrm0 * ((cos dE) - (const 1.0)) + (const 1.0) and 
      g = dt + (const 1.0)/n*((sin dE) - dE) and 
      fdot = (const 0.0) - aor*a/rnrm0*n*sindE and 
      gdot = aor*(cosdE - (const 1.0)) + (const 1.0) in 
  let r1 = Array.make 3 (const 0.0) and 
      v1 = Array.make 3 (const 0.0) in 
  for i = 0 to 2 do 
    r1.(i) <- f*r0.(i) + g*v0.(i);
    v1.(i) <- fdot*r0.(i) + gdot*v0.(i)
  done;
  (r1,v1)

(** [kepler_dqqdqv mu dt r0 v0] returns the Jacobian matrix dY/dX,
    where Y = (q1, q2) and X = (q1, v1). *)
let kepler_dqqdqv mu dt r0 v0 = 
  jacobian 
    (fun q0v0 -> 
      let r0 = Array.sub q0v0 0 3 and 
          v0 = Array.sub q0v0 3 3 in 
      let (r1,_) = kepler_advance mu dt r0 v0 in 
      Array.append r0 r1)
    (Array.append r0 v0)
