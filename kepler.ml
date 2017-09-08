open Newton

let eps = 4.0*.epsilon_float

(** [kepler_dE n dt ec es] computes the change in eccentric anamoly,
    using up to [n] iterations of a Newton-Raphson solver.  [ec] and
    [es] are e*cos(E0) and e*sin(E0), and [dt] is the corresponding
    time interval.*)
let kepler_dE n dt ec es = 
  let x0 = n*.dt -. es in 
  newton_solve 10 eps eps eps 
    (fun x -> x -. ec*.(sin x) +. es*.(1.0 -. (cos x)) -. n*.dt)
    (fun x -> 1.0 -. ec*.(cos x) +. es*.(sin x))
    x0

let dot r0 v0 = 
  let sum = ref 0.0 in 
  for i = 0 to 2 do 
    sum := !sum +. r0.(i)*.v0.(i)
  done;
  !sum +. 0.0

let norm r0 = sqrt (dot r0 r0)

(** [kepler_advance mu dt r0 v0] returns the [(r,v)] pair which
    corresponds to advancing the kepler orbit with initial condition
    [(r0,v0)] by [dt].  [mu] is the parameter in the Lagrangian: L =
    1/2 v^2 + [mu]/r *)
let kepler_advance mu dt r0 v0 = 
  let rnrm0 = norm r0 and 
      v02 = dot v0 v0 and 
      u = dot r0 v0 in 
  let a = 1.0/.(2.0/.rnrm0 -. v02/.mu) in 
  if a < 0.0 then raise (Failure "kepler_advance: does not work on unbound orbits.");
  let n = sqrt (mu /. (a*.a*.a)) in 
  let ec = 1.0 -. rnrm0/.a and 
      es = u /. (n *. a *. a) in 
  let dE = kepler_dE n dt ec es in 
  let cosdE = cos dE and 
      sindE = sin dE in 
  let aor = 1.0 /. (1.0 -. ec*.cosdE +. es*.sindE) in 
  let f = a /. rnrm0 *. ((cos dE) -. 1.0) +. 1.0 and 
      g = dt +. 1.0/.n*.((sin dE) -. dE) and 
      fdot = ~-.aor*.a/.rnrm0*.n*.sindE and 
      gdot = aor*.(cosdE -. 1.0) +. 1.0 in 
  let r1 = Array.make 3 0.0 and 
      v1 = Array.make 3 0.0 in 
  for i = 0 to 2 do 
    r1.(i) <- f*.r0.(i) +. g*.v0.(i);
    v1.(i) <- fdot*.r0.(i) +. gdot*.v0.(i)
  done;
  (r1,v1)

let dv = epsilon_float**0.2

(** [dq1dv0 mu dt r0 v0] computes (by numerical differentiation) the
    derivative of the final position with respect to the initial
    velocity.  The resulting matrix, [m], has [m.(i).(j) =
    dq.(i)/.dv0.(j)]. *)
let dq1dv0 mu dt r0 v0 = 
  let my_v0 = Array.copy v0 and 
      dq1dv0 = Array.init 3 (fun _ -> Array.make 3 0.0) in 
  for i = 0 to 2 do
     for j = 0 to 2 do 
       let {Gsl_fun.res = dq1dv0ij; err = err} = 
         Gsl_diff.central 
           (fun x -> 
             my_v0.(j) <- x;
             let (q1, _) = kepler_advance mu dt r0 my_v0 in 
             my_v0.(j) <- v0.(j);
             q1.(i))
           dv
            in 
       dq1dv0.(i).(j) <- dq1dv0ij
     done;
  done;
  dq1dv0
         
(** [dv0dq1 mu dt r0 v0] is the inverse matrix of [dq1dv0 mu dt r0
    v0]. *)
let dv0dq1 mu dt r0 v0 = 
  let dq1dv0 = dq1dv0 mu dt r0 v0 in 
  let dv0dq1_gsl = Gsl_linalg.invert_LU ~protect:true (`AA dq1dv0) in 
  Gsl_vectmat.to_arrays dv0dq1_gsl
