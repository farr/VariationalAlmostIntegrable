(** Three-body keplerian system. *)

open Diff.DFloat
open Dkepler

(** Bodies. *)
type body = {m : diff; q : diff array; v : diff array}

(** A system.  We assume that the first mass is much greater than the
    others. *)
type system = body * body * body

(** The masses are the original body masses; [qs.(0)] is the center of
  mass, the remaining [qs] are the jacobi coordinates.  [vs.(0)] is
  the velocity of the center of mass, and the remaining [vs] are the
  velocity of the Jacobi coordinates. *)
type jacobi_system = 
    {ms : diff * diff * diff;
     qs : diff array * diff array * diff array;
     vs : diff array * diff array * diff array}

let jup_sat_sys = 
  let vconv = 365.0/.(2.0*.3.1415926) in 
  ({m = C 1.0;
    q = Array.make 3 (C 0.0);
    v = Array.make 3 (C 0.0)},
   {m = C 0.000954638698; 
    q = [|C 1.941455636102168; 
          C (-4.777567725893647); 
          C (-2.361189258559663e-02)|];
    v = [|C (vconv*.6.903769195824462E-03);
          C (vconv*.3.201282408702222E-03);
          C (vconv*.(-1.677858439891209E-04))|]},
   {m = C 0.000285838546;
    q = [|C (-8.861625899291923E+00);
          C 2.913731954830950E+00;
          C 3.019596218199777E-01|];
    v = [|C (vconv*.(-2.035550738813113E-03));
          C (vconv*.(-5.310432810319549E-03));
          C (vconv*.1.738187681708300E-04)|]})

let ecc_jup_sat_sys = 
  let vconv = 365.0/.(2.0*.3.1415926) in 
  ({m = C 1.0;
    q = Array.make 3 (C 0.0);
    v = Array.make 3 (C 0.0)},
   {m = C 0.000954638698; 
    q = [|C 1.941455636102168; 
          C (-4.777567725893647); 
          C (-2.361189258559663e-02)|];
    v = [|C (0.775*.vconv*.6.903769195824462E-03);
          C (0.775*.vconv*.3.201282408702222E-03);
          C (0.775*.vconv*.(-1.677858439891209E-04))|]},
   {m = C 0.000285838546;
    q = [|C (-8.861625899291923E+00);
          C 2.913731954830950E+00;
          C 3.019596218199777E-01|];
    v = [|C (vconv*.(-2.035550738813113E-03));
          C (vconv*.(-5.310432810319549E-03));
          C (vconv*.1.738187681708300E-04)|]})

(** [jacobi_transform sys] returns the Jacobi coordinates for
    [sys]. *)
let jacobi_transform 
    ({m = m1; q = q1; v = v1},
     {m = m2; q = q2; v = v2},
     {m = m3; q = q3; v = v3}) = 
  let mtot = m1+m2+m3 and 
      m12 = m1+m2 in 
  let r0 = 
    Array.init 3 
      (fun i -> 
        (m1*q1.(i) + m2*q2.(i) + m3*q3.(i))/mtot) and 
      r1 = 
    Array.init 3
      (fun i -> q2.(i) - q1.(i)) and 
      r2 = 
    Array.init 3 
      (fun i -> 
        q3.(i) - (m1*q1.(i)+m2*q2.(i))/m12) and 
      v0 = 
    Array.init 3 
      (fun i -> 
        (m1*v1.(i) + m2*v2.(i) + m3*v3.(i))/mtot) and 
      v1 = 
    Array.init 3
      (fun i -> v2.(i) - v1.(i)) and 
      v2 = 
    Array.init 3 
      (fun i -> 
        v3.(i) - (m1*v1.(i)+m2*v2.(i))/m12) in 
  {ms = (m1,m2,m3);
   qs = (r0,r1,r2);
   vs = (v0,v1,v2)}

let jjup_sat_sys = jacobi_transform jup_sat_sys
let jecc_jup_sat_sys = jacobi_transform ecc_jup_sat_sys

let inv_jacobi_transform
    {ms = (m1,m2,m3);
     qs = (q0,q1,q2);
     vs = (v0,v1,v2)} = 
  let m12 = m1+m2 and 
      mtot = m1+m2+m3 in 
  let r3 = 
    Array.init 3 (fun i -> q0.(i) + m12*q2.(i)/mtot) and 
      r1 = 
    Array.init 3 (fun i -> q0.(i) - m2*q1.(i)/m12 - m3*q2.(i)/mtot) and 
      r2 = 
    Array.init 3 (fun i -> q0.(i) + m1*q1.(i)/m12 - m3*q2.(i)/mtot) and 
      v3 = 
    Array.init 3 (fun i -> v0.(i) + m12*v2.(i)/mtot) and 
      v1 =
    Array.init 3 (fun i -> v0.(i) - m2*v1.(i)/m12 - m3*v2.(i)/mtot) and 
      v2 = 
    Array.init 3 (fun i -> v0.(i) + m1*v1.(i)/m12 - m3*v2.(i)/mtot) in 
  ({m = m1; q = r1; v = v1},
   {m = m2; q = r2; v = v2},
   {m = m3; q = r3; v = v3})
  

let dot v1 v2 = 
  let sum = ref zero in 
  for i = 0 to 2 do 
    sum := !sum + v1.(i)*v2.(i)
  done;
  !sum

let norm v = sqrt (dot v v)

let distance_squared v1 v2 = 
  let sum = ref zero in 
  for i = 0 to 2 do 
    let tmp = v1.(i) - v2.(i) in 
    sum := !sum + tmp*tmp
  done;
  !sum

let distance v1 v2 = sqrt (distance_squared v1 v2)

let lagrangian 
    ({m = m1; q = q1; v = v1},
     {m = m2; q = q2; v = v2},
     {m = m3; q = q3; v = v3}) = 
  (const 0.5)*(m1*(dot v1 v1) +
         m2*(dot v2 v2) + 
         m3*(dot v3 v3))  +
    m1*m2/(distance q1 q2) + 
    m1*m3/(distance q1 q3) + 
    m2*m3/(distance q2 q3) 

let jacobi_perturbing_lagrangian 
    ({ms = (m1,m2,m3);
      qs = (q0,q1,q2);
      vs = (v0,v1,v2)}) = 
  let m12 = m1+m2 in 
  let q1pq2 = 
    Array.init 3 (fun i -> q1.(i)+q2.(i)) in 
  let q22 = dot q2 q2 and 
      q1pq22 = dot q1pq2 q1pq2 and 
      q2dq1pq2 = dot q2 q1pq2 and 
      m2om1 = m2 / m1 in 
  let q = (const 2.0)*m2om1*q2dq1pq2/q22 + m2om1*m2om1*q1pq22/q22 in 
  let q2coeff = one + m2om1 in 
  let denom13 = 
    Array.init 3 (fun i -> q1.(i) - q2coeff*q2.(i)) in 
  m2om1*m12*m3/(norm denom13) -
    m3*m12/(norm q2)*(q/((sqrt (one+q))+one+q))

let jacobi_lagrangian 
    ({ms = (m1,m2,m3);
      qs = (q0,q1,q2);
      vs = (v0,v1,v2)} as sys) = 
  let mtot = m1+m2+m3 and 
      m12 = m1+m2 in 
  (const 0.5)*(mtot*(dot v0 v0) +
                 m1*m2/m12*(dot v1 v1) + 
                 m12*m3/mtot*(dot v2 v2)) + 
     m1*m2/(norm q1) + 
    m12*m3/(norm q2) + 
    (jacobi_perturbing_lagrangian sys) 

let jacobi_hamiltonian 
    ({ms = (m1,m2,m3);
      qs = (q0,q1,q2);
      vs = (v0,v1,v2)} as sys) = 
  let mtot = m1+m2+m3 and 
      m12 = m1+m2 in 
  (const 0.5)*(mtot*(dot v0 v0) +
                 m1*m2/m12*(dot v1 v1) + 
                 m12*m3/mtot*(dot v2 v2)) - 
    m1*m2/(norm q1) - 
    m12*m3/(norm q2) - 
    (jacobi_perturbing_lagrangian sys) 

let jacobi_momenta {ms = (m1,m2,m3); vs = (v0,v1,v2)} = 
  let mtot = m1+m2+m3 and 
      m12 = m1+m2 in 
  (Array.map (( * ) mtot) v0,
   Array.map (( * ) (m1*m2/m12)) v1,
   Array.map (( * ) (m12*m3/mtot)) v2)

let jacobi_velocities_of_momenta {ms = (m1,m2,m3)} (p0, p1, p2) = 
  let mtot = m1+m2+m3 and 
      m12 = m1+m2 in 
  (Array.map (fun p0x -> p0x/mtot) p0,
   Array.map (fun p1x -> p1x*m12/(m1*m2)) p1,
   Array.map (fun p2x -> p2x*mtot/(m12*m3)) p2)

let jacobi_kepler_advance 
    dt {ms = (m1,m2,m3) as ms; qs = (q0,q1,q2); vs = (v0,v1,v2)} = 
  let mtot = m1+m2+m3 and 
      m12 = m1+m2 in 
  let q01 = Array.init 3 (fun i -> q0.(i) + v0.(i)*dt) and 
      (q11, v11) = kepler_advance m12 dt q1 v1 and 
      (q21, v21) = kepler_advance mtot dt q2 v2 in 
  {ms = ms; 
   qs = (q01, q11, q21);
   vs = (v0, v11, v21)}

let jacobi_kick dt ({ms = (m1,m2,m3); qs = (q0,q1,q2); vs = (v0,v1,v2)} as sys) = 
  let m12 = m1+m2 and 
      mtot = m1+m2+m3 in 
  let f1 = 
    gradient 
      (fun q1 -> jacobi_perturbing_lagrangian {sys with qs = (q0,q1,q2)})
      q1 and 
      f2 = 
    gradient
      (fun q2 -> jacobi_perturbing_lagrangian {sys with qs = (q0,q1,q2)})
      q2 in 
  let v1 = 
    Array.init 3 (fun i -> v1.(i) + dt*f1.(i)*m12/(m1*m2)) and 
      v2 = 
    Array.init 3 (fun i -> v2.(i) + dt*f2.(i)*mtot/(m12*m3)) in 
  {sys with vs = (v0,v1,v2)}
    
let jacobi_dkd dt sys = 
  jacobi_kepler_advance (dt*(C 0.5))
    (jacobi_kick dt 
       (jacobi_kepler_advance (dt*(C 0.5)) sys))

let jacobi_s4b dt sys = 
  jacobi_kick (dt*(C (1.0/.6.0)))
    (jacobi_kepler_advance (dt*(C 0.5))
       (jacobi_kick (dt*(C (2.0/.3.0)))
          (jacobi_kepler_advance (dt*(C 0.5))
             (jacobi_kick (dt*(C (1.0/.6.0)))
                sys))))

let jacobi_s6b dt sys = 
  let b1 = dt * (C (1.0 /. 12.0)) and 
      b2 = dt * (C (5.0 /. 12.0)) and 
      a1 = dt * (C (0.5*.(1.0 -. 1.0/.(Pervasives.sqrt 5.0)))) and 
      a2 = dt * (C (1.0/.(Pervasives.sqrt 5.0))) in 
  jacobi_kick b1
    (jacobi_kepler_advance a1
       (jacobi_kick b2
          (jacobi_kepler_advance a2
             (jacobi_kick b2
                (jacobi_kepler_advance a1
                   (jacobi_kick b1 sys))))))

let integrand i ({qs = (q0,q1,q2); vs = (v0,v1,v2)} as sys) t = 
  match 
    partial i 
      (fun qv -> 
        let q1 = Array.sub qv 0 3 and 
            q2 = Array.sub qv 3 3 and 
            v1 = Array.sub qv 6 3 and 
            v2 = Array.sub qv 9 3 in 
        jacobi_perturbing_lagrangian 
          (jacobi_kepler_advance (C t) {sys with qs = (q0,q1,q2); vs = (v0,v1,v2)}))
      (Array.concat [q1; q2; v1; v2]) with 
  | C x -> x
  | _ -> raise (Failure "integrand: didn't evaluate to numerical value")

let integration_eps = 16.0*.epsilon_float

let dlddqv dt ({qs = (q0,q1,q2); vs = (v0,v1,v2)} as sys) = 
  let ws = Gsl_integration.make_ws 10000 in 
  let dt_float = match dt with 
  | C x -> x
  | _ -> raise (Failure "dlddqv: cannot deal with non-constant dt") in 
  Array.init 12
    (fun i -> 
      let {Gsl_fun.res = res} = 
        Gsl_integration.qag
          (integrand i sys)
          0.0
          dt_float
          integration_eps
          integration_eps
          Gsl_integration.GAUSS61
          ws in 
      (C res))

let matrix_solve m b = 
  let m_float = 
    Array.map 
      (fun row -> 
        Array.map 
          (function 
            | C x -> x 
            | _  -> raise (Failure "matrix_solve: cannot deal with non-constant matrix"))
          row)
      m and 
      b_float = 
    Array.map 
      (function 
        | C x -> x
        | _ -> raise (Failure "matrix_solve: cannot deal with non-constant vector"))
      b in 
  let x_float = 
    Gsl_linalg.solve_LU ~protect:true (`AA m_float) (`A b_float) in 
  Array.map (fun elt -> C elt) x_float

let jac1 dt sys = 
  let (m1,m2,_) = sys.ms and 
      (_,q1,_) = sys.qs and 
      (_,v1,_) = sys.vs in 
  kepler_dqqdqv (m1+m2) dt q1 v1

let jac2 dt sys = 
  let (m1,m2,m3) = sys.ms and 
      (_,_,q2) = sys.qs and 
      (_,_,v2) = sys.vs in 
  kepler_dqqdqv (m1+m2+m3) dt q2 v2

let transpose m = 
  let ni = Array.length m and 
      nj = Array.length m.(0) in 
  Array.init nj
    (fun i -> 
      Array.init ni 
        (fun j -> 
          m.(j).(i)))

let dlddqs dt sys = 
  let dlddqv = dlddqv dt sys in 
  let dlddqv1 = Array.append (Array.sub dlddqv 0 3) (Array.sub dlddqv 6 3) and 
      dlddqv2 = Array.append (Array.sub dlddqv 3 3) (Array.sub dlddqv 9 3) in 
  let dlddqq1 = matrix_solve (transpose (jac1 dt sys)) dlddqv1 and 
      dlddqq2 = matrix_solve (transpose (jac2 dt sys)) dlddqv2 in 
  Array.concat
    [Array.sub dlddqq1 0 3;
     Array.sub dlddqq2 0 3;
     Array.sub dlddqq1 3 3;
     Array.sub dlddqq2 3 3]

let variational_evolve dt sys = 
  let eps = C (4.0*.integration_eps) in 
  let (p0,p1,p2) = jacobi_momenta sys in 
  let rec loop i sys = 
    if Pervasives.(>) i 10 then 
      raise (Failure "variational_evolve: loop not converging")
    else
      let dlddqs = dlddqs dt sys in 
      let dp1 = Array.sub dlddqs 0 3 and 
          dp2 = Array.sub dlddqs 3 3 in 
      let p1_new = Array.init 3 (fun i -> p1.(i) + dp1.(i)) and 
          p2_new = Array.init 3 (fun i -> p2.(i) + dp2.(i)) in 
      let (_,p1_old,p2_old) = jacobi_momenta sys in 
      let dp = (distance p1_new p1_old) + (distance p2_new p2_old) in 
      let dp = dp / ((norm p1_new) + (norm p2_new)) in 
      let (v0_new,v1_new,v2_new) = jacobi_velocities_of_momenta sys (p0,p1_new,p2_new) in 
      let new_sys = {sys with vs = (v0_new,v1_new,v2_new)} in 
      if dp < eps then 
        new_sys
      else
        loop (Pervasives.(+) i 1) new_sys in 
  let new_sys = loop 0 sys in 
  let evolved_sys = jacobi_kepler_advance dt new_sys in 
  let (p0,p1,p2) = jacobi_momenta evolved_sys and 
      dlddqs = dlddqs dt new_sys in 
  let dp1 = Array.sub dlddqs 6 3 and 
      dp2 = Array.sub dlddqs 9 3 in 
  let p1_new = Array.init 3 (fun i -> p1.(i) + dp1.(i)) and 
      p2_new = Array.init 3 (fun i -> p2.(i) + dp2.(i)) in 
  let (v0,v1,v2) = jacobi_velocities_of_momenta new_sys (p0,p1_new,p2_new) in 
  {evolved_sys with vs = (v0,v1,v2)} 
