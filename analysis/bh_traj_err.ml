open Dthree_body
open Diff.DFloat

let t_tot = 1000.0

let rec pow_int_int n m = 
  if Pervasives.(=) m 0 then 
    1
  else if Pervasives.(=) m 1 then 
    n
  else if Pervasives.(=) (m mod 2) 0 then  
    let nm2 = pow_int_int n (Pervasives.(/) m 2) in 
    Pervasives.( * ) nm2 nm2
  else
    Pervasives.( * ) n (pow_int_int n (Pervasives.(-) m 1))

let ns = List.map (fun x -> pow_int_int 2 x) [0; 1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12; 13; 14; 15]

let dts = List.map (fun n -> t_tot /. (float_of_int n)) ns

let jsys = 
  let (s,j,sat) = jup_sat_sys in 
  let bh_sys = (s, {j with m = const 1e-6}, {sat with m = const 2e-7}) in 
  jacobi_transform bh_sys

let _ = 
  let {ms = (m1,m2,_);
       qs = (_,qjup,_);
       vs = (_,vjup,_)} = jsys in 
  Printf.fprintf stdout "Jup e = ";
  print_diff stdout (Dkepler.ecc (m1+m2) qjup vjup);
  Printf.fprintf stdout "\n%!"

let exact_final_state = 
  let dt = List.nth dts 12 and (* This is the min energy error step *)
      n = List.nth ns 12 in 
  let s = ref jsys in 
  for i = 1 to n do 
    s := jacobi_s6b (C dt) !s
  done;
  !s
  
let state_distance 
    {qs = (q11, q12, q13);
     vs = (v11, v12, v13)}
    {qs = (q21, q22, q23);
     vs = (v21, v22, v23)} =
  sqrt (((distance_squared q11 q21) +
           (distance_squared q12 q22) + 
           (distance_squared q13 q23) + 
           (distance_squared v11 v21) + 
           (distance_squared v12 v22) + 
           (distance_squared v13 v23)) /
        ((dot q11 q11) + 
           (dot q12 q12) + 
           (dot q13 q13) + 
           (dot v11 v11) + 
           (dot v12 v12) +
           (dot v13 v13)))

let test_integrator evolve outfile = 
  let output = open_out outfile in
  List.iter2 
    (fun n dt -> 
      let s = ref jsys in 
      for i = 1 to n do
        s := evolve (C dt) !s
      done;
      let dqp = state_distance exact_final_state !s in 
      match dqp with 
      | C dqp -> Printf.fprintf output "%g %g\n%!" dt dqp
      | _ -> raise (Failure "state distance did not produce numerical output"))
    ns
    dts;
  close_out output

let _ = test_integrator jacobi_dkd "output/bh-traj-err/dkd.dat"
let _ = test_integrator jacobi_s4b "output/bh-traj-err/s4b.dat"
let _ = test_integrator jacobi_s6b "output/bh-traj-err/s6b.dat"
let _ = test_integrator variational_evolve "output/bh-traj-err/var.dat"
