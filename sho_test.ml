open Sho
open OUnit

let random_state () = 
  {q = Random.float 1.0 -. 0.5;
   p = Random.float 1.0 -. 0.5}

let test_evolve () = 
  for i = 1 to 100 do 
    let s1 = random_state () in 
    let e = energy s1 and 
        s2 = evolve 10.0 s1 in 
    let e2 = energy s2 in 
    assert_bool "evolve breaks energy conservation" 
      ((abs_float (e2 -. e))/.e < 100.0*.epsilon_float)
  done

let test_evolver ev dt eg ratio = 
  for i = 1 to 100 do 
    let s0 = random_state () in 
    let e0 = eg s0 and 
        s1 = ev dt s0 and 
        s2 = ev (0.5*.dt) (ev (0.5*.dt) s0) in 
    let e1 = eg s1 and 
        e2 = eg s2 in 
    let de1 = (abs_float (e1 -. e0))/.e0 and 
        de2 = (abs_float (e2 -. e0))/.e0 in 
    let r = de1 /. de2 in 
    assert_bool (Printf.sprintf "energy error ratio (%g = %g/%g) outside limits"
                   r de1 de2)
      (r > ratio/.(sqrt 2.0) && r < ratio*.(sqrt 2.0))
  done

let test_dkd () = 
  let eps = 1e-4 in 
  test_evolver (dkd eps) 0.01 (energy_pert eps) 4.0

let test_kdk () = 
  let eps = 1e-4 in 
  test_evolver (kdk eps) 0.01 (energy_pert eps) 4.0

let test_s4b () = 
  let eps = 1e-6 in 
  test_evolver (s4b eps) 0.5 (energy_pert eps) 16.0

let test_s6b () = 
  let eps = 1e-6 in 
  test_evolver (s6b eps) 1.0 (energy_pert eps) 64.0

let test_variational () = 
  for i = 1 to 100 do 
    let eps = Random.float 1e-4 in 
    let s = random_state () in 
    let s2 = 
      let s = ref s in 
      for i = 1 to 100 do 
        s := variational_advance eps 0.01 !s
      done;
      !s in 
    let e = energy_pert eps s and 
        e2 = energy_pert eps s2 in 
    let de = abs_float ((e2 -. e)/.e) in 
    assert_bool (Printf.sprintf "energy error too big: %g" de)
      (de < 1e-8)
  done
      

let tests = "sho.ml tests" >:::
  ["evolve test" >:: test_evolve;
   "dkd test" >:: test_dkd;
   "kdk test" >:: test_kdk;
   "s4b test" >:: test_s4b;
   "s6b test" >:: test_s6b;
   "variational test" >:: test_variational]
