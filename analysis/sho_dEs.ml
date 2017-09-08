open Sho

let t_tot = 1000.0 

let rec pow_int_int n m = 
  if m = 0 then 
    1
  else if m = 1 then 
    n
  else if m mod 2 = 0 then  
    let nm2 = pow_int_int n (m/2) in 
    nm2*nm2
  else
    n*(pow_int_int n (m-1))

let ns = List.map (fun x -> pow_int_int 2 x) [0; 1; 2; 3; 4; 5; 6; 7; 8; 9; 10; 11; 12; 13; 14; 15]

let dts = List.map (fun n -> t_tot /. (float_of_int n)) ns

let eps = 1e-5

let s0 = {q = 1.0; p = 0.0}

let e0 = energy_pert eps s0

let test_integrator evolve outfile = 
  let output = open_out outfile in 
  List.iter2 
    (fun n dt -> 
      let dE_max = ref neg_infinity and 
          s = ref s0 in 
      for i = 1 to n do 
        s := evolve eps dt !s;
        let e = energy_pert eps !s in 
        let dE = abs_float ((e -. e0)/.e0) in 
        if dE > !dE_max then dE_max := dE
      done;
      Printf.fprintf output "%g %g\n" dt !dE_max)
    ns
    dts;
  close_out output

let _ = test_integrator dkd "output/small-eps-dE/dkd.dat"
let _ = test_integrator s4b "output/small-eps-dE/s4b.dat"
let _ = test_integrator s6b "output/small-eps-dE/s6b.dat"
(* let _ = test_integrator variational_advance_dkd "output/small-eps-dE/var-dkd.dat" *)
(* let _ = test_integrator variational_advance_gl3 "output/small-eps-dE/var-gl3.dat" *)
let _ = test_integrator variational_advance "output/small-eps-dE/var.dat"
