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

let e0 = jacobi_hamiltonian jsys

let _ = 
  let {ms = (m1,m2,_);
       qs = (_,qjup,_);
       vs = (_,vjup,_)} = jsys in 
  Printf.fprintf stdout "Jup e = ";
  print_diff stdout (Dkepler.ecc (m1+m2) qjup vjup);
  Printf.fprintf stdout "\n%!"

let test_integrator evolve outfile = 
  let output = open_out outfile in
  List.iter2 
    (fun n dt -> 
      let dE_max = ref neg_infinity and 
          s = ref jsys in 
      for i = 1 to n do 
        s := evolve (C dt) !s;
        let e = jacobi_hamiltonian !s in 
        let e0 = match e0 with 
        | C x -> x
        | _ -> raise (Failure "test_integrator: got bad e0 value") and 
            e = match e with 
            | C x -> x
            | _ -> raise (Failure "test_integrator: got bad e value") in 
        let dE = abs_float ((e -. e0)/.e0) in 
        if Pervasives.(>) dE !dE_max then dE_max := dE
      done;
      Printf.fprintf output "%g %g\n%!" dt !dE_max)
    ns
    dts;
  close_out output

let _ = test_integrator jacobi_dkd "output/bh-dE/dkd.dat"
let _ = test_integrator jacobi_s4b "output/bh-dE/s4b.dat"
let _ = test_integrator jacobi_s6b "output/bh-dE/s6b.dat"
let _ = test_integrator variational_evolve "output/bh-dE/var.dat"
